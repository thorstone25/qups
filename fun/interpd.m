function y = interpd(x, t, dim, interp)
% INTERPD Interpolation in one dimension
%
% y = interpd(x, t, dim, interp)
%
% x is (T x N')
% t is (I x N' x M')
%
% Sampling is performed in dimension matching I/T, broadcasted in N' (must
% match) and the data is replicated across M' i.e. x(t, n', m') = x(t, n')
%
% 

%% validate dimensions
if nargin < 4, interp = 'linear'; end
if nargin < 3, dim = 1; end
assert(isreal(t)); % sampling at a real data type

maxdims = max(ndims(x), ndims(t));
ddms = setdiff(1:maxdims, dim); % data dimensions
xsz = size(x, ddms);
tsz = size(t, ddms);
assert(all(xsz == tsz | xsz == 1), 'Data sizing is not broadcastable.'); % data sizing must be broadcastable

bdms = ddms(xsz == tsz); % matching dimensions
rdms = ddms(xsz ~= tsz & xsz == 1); % outer dimensions

% set the data order: sampling dimensions, broadcasting dimensions,
% replicating dimensions
ord = [dim, bdms, rdms]; % this should be all of the dimensions

% get the data sizing
T = size(x, dim);
I = size(t, dim);
if isempty(bdms), N = 1; else, N = prod(size(t, bdms)); end
if isempty(rdms), M = 1; else, M = prod(size(t, rdms)); end

% move the input data to the proper dimensions for the GPU kernel
x = permute(x, ord);

% grab the kernel reference
issingle = @(x) strcmp(class(x), 'single') || any(ismember(class(x), ['tall', 'gpuArray'])) && strcmp(classUnderlying(x), 'single'); %#ok<STISA> 
if issingle(x) || issingle(t)
    [x, t] = deal(single(x), single(t));
    suffix = "f"; 
else
    suffix = "";
end

k = parallel.gpu.CUDAKernel('bf.ptx', 'interpd.cu', 'interpd' + suffix); % TODO: fix path
k.setConstantMemory('I', uint64(I), 'T', uint64(T), 'N', uint64(N), 'M', uint64(M));
k.ThreadBlockSize = k.MaxThreadsPerBlock; % why not?
k.GridSize = [ceil(I ./ k.ThreadBlockSize(1)), N, 1]; % no need to go further
switch interp
    case "nearest", flagnum = 0;
    case "linear",  flagnum = 1;
    case "cubic",  flagnum = 2;
    case "lanczos3",flagnum = 3;
end

% sample
y = nan(size(t), 'like', x);
y = k.feval(y, x, t, flagnum); 

% place back in prior dimensions
y = ipermute(y, ord);



