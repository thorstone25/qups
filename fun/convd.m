function [C, lags] = convd(A, B, dim, shape, varargin)
% CONVD - GPU-enabled Convolution in one dimension
%
% C = CONVD(A, B) computes the convolution of A with B in dimension 1. A 
%   and B may differ in size in dimension 1 only. All other dimensions must
%   be of equal size.
% C = CONVD(A, B, dim) executes in dimension "dim" instead of dimension 1
% C = CONVD(..., shape) selects the shape of the returned convolution. Can 
%   be one of {'full'*|'same'|'valid'}. The default is 'full'.
% C = CONVD(..., 'device', dev) selects which gpu device to use. 
%   dev = -1 specifies the device returned by gpuDevice() 
%   dev = n where 1 <= n <= gpuDeviceCount selects device n.
%   dev = 0 specifies no device and operates in native MATLAB code (default) 
%
% [C, lags] = CONVD(...) returns the lags of y in the same dimension as the
% computation.
% 
% To select default behaviour, pass an empty argument.
%
% See also CONV CONV2 CONVN

% TODO: switch to kwargs structure

% parse the inputs and set defaults
if nargin < 3 || isempty(dim), dim = 1; end
if nargin < 4 || isempty(shape), shape = 'full'; end
device = 0;

% get optional inputs
for i = 1:2:numel(varargin)
    switch varargin{i}
        case 'device'
            device = varargin{i+1};
    end
end

% permute so that dim becomes the first dimension
% flip 2nd arg on CPU so that stride is positive
ord = [dim, setdiff(1:max(ndims(A), dim), dim)];
x =      permute(A, ord); % shared-copy?
y = flip(permute(B, ord),1); % shared-copy?

% check data sizes
sz_x = size(x);
sz_y = size(y);
assert(numel(sz_x) == numel(sz_y) && all(sz_x(2:end) == sz_y(2:end)),...
    'Incompatible sizes.');

% get computation/output data type and precision
complex_type = ~isreal(A) || ~isreal(B);
gpu_type    = isa(A, 'gpuArray') || isa(B, 'gpuArray');
single_type = isa(A, 'single'  ) || isa(B, 'single'   ) ...
    || isa(A, 'gpuArray') && strcmp(classUnderlying(A), 'single') ...
    || isa(B, 'gpuArray') && strcmp(classUnderlying(B), 'single');

% get function output prototype
To = zeros(0);
if single_type, To = single(To); end
if gpu_type, To = gpuArray(To); end
if complex_type, To = complex(To); end

% get data/kernel sizing info (M >= N)
M = size(x,1); 
N = size(y,1);

% get the proper lags for this computation
switch shape
    case 'full'
        lags = colon(-(N - 1), M - 1).';
    case 'same'
        lags = colon(0,        M - 1).' - floor(N/2);
    case 'valid'
        lags = colon(0,        M - N).';
end

% get kernel sizing info
if isempty(lags), l0 = 0; else, l0 = -lags(1); end
L = numel(lags); % number of lags
S = prod(sz_x(2:end)); % number of slices
sz = [L, sz_x(2:end)]; % output is lags by size of A

% operate on GPU
if device
    if device > 0
        g = gpuDevice(device); % access a specific gpu device
    else
        g = gpuDevice; % access the current gpu device
    end
    
    % get the kernel suffix
    suffix = '';
    if complex_type, suffix = [suffix 'c']; end
    if single_type,  suffix = [suffix 'f']; end
    
    % the kernel file
    src.folder = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
    src.name = 'conv_cuda';
    
    % specify the kernel
    kern = parallel.gpu.CUDAKernel(...
        [src.name '.ptx'],... % must be on the path
        fullfile(src.folder, [src.name '.cu']), ... % must be in source
        ['conv' suffix]);
    
    % setup the execution size
    lags_per_thread = g.MaxThreadsPerBlock;
    kern.ThreadBlockSize = [lags_per_thread, 1, 1];
    kern.GridSize = [ceil(L / lags_per_thread), S, 1];
    
    % define the constant data parameters: care must be taken to match the
    % datatype here: for NVIDIA gpus, size_t <-> uint64
    kern.setConstantMemory('M', uint64(M), 'N', uint64(N), 'L', uint64(L), 'L0', int32(l0));
    
    % if complex, ensure both input arguments complex
    if complex_type && isreal(x), x = complex(x); end
    if complex_type && isreal(y), y = complex(y); end
        
    % allocate the output
    switch shape
        case 'full',  z = zeros(sz, 'like', gpuArray(To)); % explicit pre-allocation - z larger than x
        case 'same',  z = x; % implicit pre-allocation - z will be the same size as x
        case 'valid', z = sub(x, 1:L, 1); % implicit pre-allocation - z will be smaller than x
    end
        
    % queue the kernel
    z = kern.feval(x, y, z);
    
else
    % vectorized MATLAB on CPU - perform convolution manually for vectors
    % parfor handles implicit allocation and avoids data copies

    % TODO: use the current pool if not on a GPU
    
    % just for-loop it
    parfor (s = 1:S, 0), z(:,s) = conv(x(:,s), y(:,s), shape); end
end

% cast to output type / dimensions - try to use implicit in-place assignment
z = ipermute(reshape(z, sz), ord);
z = cast(z, 'like', To); 
C = z; % shared-copy
lags = ipermute(lags, ord);

