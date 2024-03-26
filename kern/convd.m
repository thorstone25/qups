function [C, lags] = convd(x, y, dim, shape, kwargs)
% CONVD - GPU-enabled Convolution in one dimension
%
% C = CONVD(A, B) computes the convolution of A with B in dimension 1. A 
% and B may differ in size in dimension 1 only. All other dimensions must
% be of equal size.
% 
% C = CONVD(A, B, dim) executes in dimension "dim" instead of dimension 1
% 
% C = CONVD(A, B, dim, shape) selects the shape of the output. Must be one 
% of {'full'*, 'same', 'valid'}. The default is 'full'.
% 
% C = CONVD(..., 'gpu', true) selects whether to use a gpu. A ptx-file will 
% be used if compiled. The default is true if x or y is a gpuArray.
%
% C = CONVD(..., 'ocl', true) selects whether to use an OpenCL device if
% OpenCL support is available. If the currently selected device does not 
% support the precision of the data, this argument is ignored. The default
% is true if OpenCL support is available via Matlab-OpenCL.
%
% C = CONVD(..., 'parenv', clu) or C = CONVD(..., 'parenv', pool) or 
% performs the convolution in parallel on the parcluster clu or the parpool
% pool when operating with native MATLAB code. If clu is 0, no parallel
% environment is used. The default is the current parallel pool returned
% by gcp.
%
% [C, lags] = CONVD(...) returns the lags of y in the same dimension as the
% computation.
% 
% Example:
% % correlate an integer sequence
% y = convd([1 -2 3 -4 5], [5 -4 3 -2 1])
% isalmostn(y, xcorr([1 -2 3 -4 5]))
% 
% % Create two 16 x 4 exponential sequences
% n = (1:16)-1;
% m = (1:4)'-1;
% xa = single(0.84).^(n+m);
% xb = single(0.92).^(n+m);
% 
% % Compute iteratively
% z0 = 0*xa; % initialize
% for i = 1:4
%    z0(i,:) = conv(xa(i,:), xb(i,:), 'same');
% end
% 
% % Compute with convd
% z = convd(xa,xb,2,'same');
% isalmostn(z, z0)
% 
% See also CONV CONV2 CONVN

arguments
    x {mustBeFloat}
    y {mustBeFloat}
    dim (1,1) {mustBePositive, mustBeInteger} = findSingletonDim(x, y)
    shape (1,1) string {mustBeMember(shape, ["full", "same", "valid"])} = 'full'
    kwargs.gpu (1,1) logical = isa(x, 'gpuArray') || isa(y, 'gpuArray')
    kwargs.ocl (1,1) logical = exist('oclDeviceCount', 'file') && oclDeviceCount() && ~isempty(oclDevice())
    kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate') % parallel environment
    kwargs.lowmem (1,1) logical = 16 * max(numel(x), numel(y)) > 16 * 2^30; % default true if complex double storage of either argument exceeds 16GB
end

% permute so that dim becomes the first dimension
% flip 2nd arg on CPU so that stride is positive
D = max(ndims(x),ndims(y)); % number of dimensions in the data
idim = setdiff(1:max(D, dim), dim); % inverse (not) dim

% check data sizes
sz_x = size(x, 1:D);
sz_y = size(y, 1:D);
assert(numel(sz_x) == numel(sz_y) && all(sz_x(idim) == sz_y(idim)),...
    "Incompatible sizes [" + join(string((sz_x))+",") + "], and [" + join(string((sz_y))+",") + "].");

% get computation/output data type and precision
complex_type = ~isreal(x) || ~isreal(y);
gpu_type    = isa(x, 'gpuArray') || isa(y, 'gpuArray');
dtype = getDtype(x, y); % get the data type using casting rules

% get function output prototype
To = zeros(0);
dfun = str2func(dtype); % casting function / constructor
To = dfun(To);
if gpu_type, To = gpuArray(To); end
if complex_type, To = complex(To); end

% get data/kernel sizing info (M >= N)
M = size(x,dim); 
N = size(y,dim);

% get the proper lags for this computation
switch shape
    case 'full'
        lags = colon(-(N - 1), M - 1);
    case 'same'
        lags = colon(0,        M - 1) - floor((N-1)/2);
    case 'valid'
        lags = colon(0,        M - N);
end

% get kernel sizing info
if isempty(lags), l0 = 0; else, l0 = -lags(1); end
L = numel(lags); % number of lags
S = prod(sz_x(dim+1:D)); % number of slices


% row strides
C = prod(esize(x, 1:dim-1));
[xstr, ystr, zstr] = deal(C);
sizes = [xstr, M, ystr, N, zstr, L, C];

% whether/how to operate with CUDA/OpenCL
use_gdev = kwargs.gpu && exist('convd.ptx', 'file') && exist('convd.cu', 'file');
use_odev = kwargs.ocl && exist('oclDeviceCount', 'file') && oclDeviceCount() && ~isempty(oclDevice()) && exist('convd.cl', 'file');
if use_odev % validate the device supports the type
    dev = oclDevice(); % get device
    switch dtype
        case 'double', use_odev = dev.SupportsDouble;
        case 'halfT' , use_odev = dev.SupportsHalf;
    end
end

% dispatch based on device (native vs. kernel)
if use_gdev || use_odev
    % if complex, eshiftdimnsure both all arguments complex
    if complex_type && isreal(x), x = complex(x); end
    if complex_type && isreal(y), y = complex(y); end

    % ensure all have dtype type
    x = dfun(x); 
    y = dfun(y);

    % allocate the output implicitly
    switch shape
        case 'valid', z = sub(x, 1:L, dim); % implicit pre-allocation - z will be smaller than x (shared copy ?)
        case 'same',  z = x; % implicit pre-allocation - z will be the same size as x - (shared copy)
        case 'full'
            z_ex_sz = size(x);
            z_ex_sz(dim) = L - size(x, dim); % expansion size
            x = cat(dim, x, zeros(z_ex_sz, 'like', x)); % implict pre-allocation - maintain valid region by zero-padding
            z = x;
            x = sub(z, 1:M, dim); % (shared-copy ?)
    end
    
    % get the kernel
    if use_gdev
        % get the kernel suffix
        suffix = '';
        if complex_type, suffix = [suffix 'c']; end
        switch dtype
            case 'single', suffix = [suffix 'f'];
            case 'halfT',  suffix = [suffix 'h'];
        end

        % specify the kernel
        kern = parallel.gpu.CUDAKernel( ...
            'convd.ptx', 'convd.cu', ['conv' suffix] ...
            );

        % define the constant data parameters: care must be taken to match the
        % datatype here: for NVIDIA gpus, size_t <-> uint64
        kern.setConstantMemory('L0', int32(l0));
    
    else
        % reference kernel
        kern = oclKernel('convd.cl');

        % define typing
        switch dtype
            case 'single', prc = 32; typ = 'single';
            case 'double', prc = 64; typ = 'double';
            case 'halfT' , prc = 16; typ = 'half'  ;
        end
        kern.defineTypes({typ}); % {T}
        kern.macros = [kern.macros, "QUPS_PRECISION="+prc];
        if complex_type, kern.macros = [kern.macros, "QUPS_COMPLEX"]; end

        % define the constant data parameters
        kern.macros = [kern.macros, "QUPS_L0="+int32(l0)];
    end
    
    % setup the execution size
    if C == 1, blk = [1 min(L, kern.MaxThreadsPerBlock) 1];
    else,      blk = [min(C, kern.MaxThreadsPerBlock), 1, 1];
    end
    kern.ThreadBlockSize = blk;
    kern.GridSize = ceil([C L S] ./ blk);
    
    % run the kernel
    z = kern.feval(x,y,z,sizes); 
    
else
    % vectorized MATLAB on CPU - perform convolution manually for vectors
    % parfor handles implicit allocation and avoids data copies

    % move to GPU if requested
    if kwargs.gpu, [x, y] = deal(gpuArray(x), gpuArray(y)); end

    % Use the given parallel environment if not on a GPU
    clu = kwargs.parenv;
    if isempty(clu) || isa(x, 'gpuArray') || isa(y, 'gpuArray'), clu = 0; end

    % get vectorized sizing
    sz = size(x);
    idim = 1:dim-1;
    jdim = dim+1:ndims(x);
    I = prod(sz(idim));
    J = prod(sz(jdim));

    % vectorize
    x = reshape(x, [I M J]);
    y = reshape(y, [I N J]);
    zproto = x(1)*y(1);
    z = zeros([I L J], 'like', zproto);

    if kwargs.lowmem || (J < I && J <= 64)  % J small - par over I
        for (j = 1:J)
            parfor(i = 1:I, clu)
                zj(i,:) = convn(x(i,:,j), y(i,:,j), shape);
            end
            z(:,:,j) = zj;
        end
    else % J large - par over J
        parfor (j = 1:J, clu)
            zj = zeros([I L], 'like', zproto);
            for(i = 1:I)
                zj(i,:) = convn(x(i,:,j), y(i,:,j), shape);
            end
            z(:,:,j) = zj;
        end
    end

    % return to original sizing
    sz(dim) = L;
    z = reshape(z, sz);
end

% cast to output type / dimensions - try to use implicit in-place assignment
% z = cast(z, 'like', To); 
C = z; % shared-copy
lags = swapdim(lags, 2, dim); % put lags in same dimensions as operation

function d = findSingletonDim(A, B)
dA = find(size(A) ~= 1,1,'first');
dB = find(size(B) ~= 1,1,'first');
d = min([dA, dB]);
if isempty(d), d = 1; end

function dtype = getDtype(A, B)

single_type = isa(gather(A(1:0)), 'single') || isa(gather(B(1:0)), 'single');
half_type = isa(A, 'halfT') || isa(B, 'halfT');

if half_type,       dtype = 'halfT';
elseif single_type, dtype = 'single'; 
else,               dtype = 'double';
end

function sz = esize(x, dim), if isempty(dim), sz = []; else, sz = size(x, dim); end