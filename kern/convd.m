function [C, lags] = convd(A, B, dim, shape, varargin)
% CONVD - GPU-enabled Convolution in one dimension
%
% C = CONVD(A, B) computes the convolution of A with B in dimension 1. A 
%   and B may differ in size in dimension 1 only. All other dimensions must
%   be of equal size.
% 
% C = CONVD(A, B, dim) executes in dimension "dim" instead of dimension 1
% 
% C = CONVD(A, B, dim, shape) selects the shape of the output. Must be one 
%   of {'full'*|'same'|'valid'}. The default is 'full'.
% 
% C = CONVD(..., 'device', dev) selects which gpu device to use. 
%   dev = -1 specifies the current device (returned by gpuDevice())
%   dev = n where 1 <= n <= gpuDeviceCount selects device n.
%   dev = 0 specifies no device and operates in native MATLAB code (default) 
%
% C = CONVD(..., 'parcluster', clu) performs the convolution in parallel on
%   the parcluster clu. The default is the current parallel pool. If the
%   data is on a GPU, the cluster object is not used.
%
% [C, lags] = CONVD(...) returns the lags of y in the same dimension as the
% computation.
% 
% To select default behaviour, pass an empty argument.
%
% See also CONV CONV2 CONVN

% TODO: update this for tall variables?
% TODO: extend for half type

% parse the inputs and set defaults
if nargin < 3 || isempty(dim), dim = 1; end
if nargin < 4 || isempty(shape), shape = 'full'; end

% defaults
kwargs.device = 0; 
kwargs.parcluster = gcp('nocreate'); 

% get optional inputs
for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

% permute so that dim becomes the first dimension
% flip 2nd arg on CPU so that stride is positive
D = max(ndims(A),ndims(B)); % number of dimensions in the data
ord = [dim, setdiff(1:max(D, dim), dim)];
x = permute(A, ord); % shared-copy?
y = permute(B, ord); % shared-copy?

% check data sizes
sz_x = size(x, 1:D);
sz_y = size(y, 1:D);
assert(numel(sz_x) == numel(sz_y) && all(sz_x(2:end) == sz_y(2:end)),...
    "Incompatible sizes [" + join(string((sz_x))+",") + "], and [" + join(string((sz_y))+",") + "].");

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

% specify the kernel file
src.folder = fullfile(fileparts(mfilename('fullpath')), '..', 'src');
src.name = 'conv_cuda';

% initialize the GPU if requested
if kwargs.device > 0
    g = gpuDevice(kwargs.device); % access a specific gpu device
elseif kwargs.device < 0
    g = gpuDevice; % access the current gpu device
end

% whether/how to operate on GPU
if ~isempty(kwargs.device) && kwargs.device && exist([src.name '.ptx'], 'file')
    % get the kernel suffix
    suffix = '';
    if complex_type, suffix = [suffix 'c']; end
    if single_type,  suffix = [suffix 'f']; end
    
    % specify the kernel
    kern = parallel.gpu.CUDAKernel(...
        [src.name '.ptx'],... % must be on the path
        fullfile(src.folder, [src.name '.cu']), ... % must be in source
        ['conv' suffix]); % function name in the kernel
    
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
        
    % run the kernel
    z = kern.feval(x, flip(y,1), z); % y is flipped for the kernel
    
else
    % vectorized MATLAB on CPU - perform convolution manually for vectors
    % parfor handles implicit allocation and avoids data copies

    % move to GPU if requested
    if kwargs.device, [x, y] = deal(gpuArray(x), gpuArray(y)); end

    % Use the current pool if not on a GPU
    clu = kwargs.parcluster;
    if isempty(clu) || isa(x, 'gpuArray') || isa(y, 'gpuArray'), clu = 0; end
    
    % for-loop it in, parallel if a pool exists
    parfor (s = 1:S, clu), z(:,s) = conv(x(:,s), y(:,s), shape); end
end

% cast to output type / dimensions - try to use implicit in-place assignment
z = ipermute(reshape(z, sz), ord);
z = cast(z, 'like', To); 
C = z; % shared-copy
lags = ipermute(lags, ord); % put lags in same dimensions as operation

