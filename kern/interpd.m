function [y, k, ord] = interpd(x, t, dim, interp, extrapval)
% INTERPD GPU-enabled interpolation in one dimension
%
% y = INTERPD(x, t) interpolates the data x at the indices t. It uses the
% matching and non-matching dimensions of the data to implicitly broadcast
% operations across the data. In other words, the data is formatted into
% the following:
%
%   x is (T x N' x      F')
%   t is (I x N' x M'     )
%   y is (I x N' x M' x F')
%
% Sampling is performed in the dimension matching I/T, element-wise in N'
% (where matched) and broadcasted across M' and F' such that x(k,n',m',f')
% = x(k,n',f') and t(i,n',m',f') = t(i,n',m').
%
% x is 0-base indexed implicitly i.e. x(1) in MATLAB <-> t == 0
%
% y = INTERPD(x, t, dim) samples the data in dimension dim instead of
% dimension 1
%
% y = INTERPD(x, t, dim, interp) specifies the interpolation method. It
% must be one of {"nearest", "linear"*,"cubic","lanczos3"} or any option 
% supported by interp1.
%
% y = INTERPD(x, t, dim, interp, extrapval) uses extrapval as the
% extrapolation value when outside of the domain of t.
%
% See also INTERPN INTERP1
%

%% validate dimensions
if nargin < 5 || isempty(extrapval), extrapval = nan; end
if nargin < 4 || isempty(interp),    interp = 'linear'; end
if nargin < 3 || isempty(dim),       dim = 1; end
assert(isreal(t)); % sampling at a real data type
assert(isscalar(extrapval), 'Only a scalar value accepted for extrapolation.');

k = parallel.gpu.CUDAKernel.empty; % default CUDAKernel output argument

maxdims = max(ndims(x), ndims(t)); % maximum dimension
ddms = setdiff(1:maxdims, dim); % data dimensions
xsz = size(x, ddms);
tsz = size(t, ddms);
mdms  = ddms(xsz == tsz); % matching dimensions
rdmsx = ddms(xsz ~= tsz & tsz == 1); % outer dim for x
rdmst = ddms(xsz ~= tsz & xsz == 1); % outer dim for t
rdms = union(rdmsx, rdmst); % all outer dimensions

% set the data order: sampling dimensions, broadcasting dimensions,
% replicating dimensions (in t, then in x)
ord = [dim, mdms, rdmst, rdmsx]; % this should be all of the dimensions
assert(isempty(setxor(ord, 1:maxdims)), 'Unable to identify all dimensions: this is a bug.');

% get the data sizing
T = size(x, dim);
I = size(t, dim);
N = prod(esize(t, mdms));
M = prod(esize(t, rdms));
F = prod(esize(x, rdms));

% move the input data to the proper dimensions for the GPU kernel
if ~isequal(ord, 1:length(ord)) % avoid data copy if already ordered
x = permute(x, ord);
t = permute(t, ord);
end

% function to determine type
isftype = @(x,T) strcmp(class(x), T) || any(arrayfun(@(c)isa(x,c),["tall", "gpuArray"])) && strcmp(classUnderlying(x), T);

if exist('oclDeviceCount','file') && oclDeviceCount(), ocl_dev = oclDevice(); else, ocl_dev = []; end % get oclDevice if supported

use_gdev = exist('interpd.ptx', 'file') ...
        && ( isa(x, 'gpuArray') || isa(t, 'gpuArray') || isa(x, 'halfT') && x.gtype) ...
        && (ismember(interp, ["nearest", "linear", "cubic", "lanczos3"])); ... 

use_odev = exist('interpd.cl', 'file') ...
        && (ismember(interp, ["nearest", "linear", "cubic", "lanczos3"])) ... 
        && ~isempty(ocl_dev) && ...
        (  (isftype(x,'double') && ocl_dev.SupportsDouble) ...
        || (isftype(x,'single')                          ) ... always supported
        || (isftype(x,'halfT' ) && ocl_dev.SupportsHalf  ));

% use ptx on gpu if available or use native MATLAB
if use_gdev || use_odev
    
    % determine the data type
    if     isftype(x, 'double')
        suffix = "" ; prc = 64; [x,t] = dealfun(@double, x, t);
    elseif isftype(x, 'single')
        suffix = "f"; prc = 32; [x,t] = dealfun(@single, x, t);
    elseif isftype(x, 'halfT' )
        suffix = "h"; prc = 16; [x,t] = dealfun(@(x)gpuArray(halfT(x)), x, t); % custom type
    else
        suffix = "f"; prc = 32;
        warning("Datatype " + class(x) + " not recognized as a GPU compatible type.");
    end

    % translate the interp flag
    switch interp
        case "nearest", flagnum = 0;
        case "linear",  flagnum = 1;
        case "cubic",   flagnum = 2;
        case "lanczos3",flagnum = 3;
        otherwise, error('Interp option not recognized: ' + string(interp));
    end

    % cache data prototype
    x0 = zeros(0, 'like', x);

    % condition inputs/outputs
    osz = [I, max(size(t,2:maxdims), size(x,2:maxdims))];
    x = complex(x); % enforce complex type
    switch suffix
        case "h" % halfT type
            y = complex(gpuArray(halfT(repelem(extrapval,osz))));
            [y_,x_,t_] = dealfun(@(x)x.val, y,x,t);
        otherwise % others
            y = complex(repmat(cast(extrapval, 'like', x), osz));
            [y_,x_,t_] = deal(y,x,t);
    end
    
    if use_gdev
        % grab the kernel reference
        k = parallel.gpu.CUDAKernel('interpd.ptx', 'interpd.cu', 'interpd' + suffix);
        d = gpuDevice();
        
        % set constants
        k.setConstantMemory('QUPS_I', uint64(I), 'QUPS_T', uint64(T), 'QUPS_N', uint64(N), 'QUPS_M', uint64(M), 'QUPS_F', uint64(F));
    
    elseif use_odev
        % get the kernel reference
        k = oclKernel(which('interpd.cl'), 'interpd');
        d = k.Device;

        % set the data types
        switch prc, case 16, tp = 'uint16'; case 32, tp = 'single'; case 64, tp = 'double'; end
        k.defineTypes({tp,tp}); % all aliases are this type

        % configure the kernel sizing and options
        k.macros = "QUPS_" + ["I", "T", "S", "N", "F"] + "=" + uint64([I T M N F]); % size constants
        k.macros = [k.macros, "QUPS_INTERPD_" + ["NO_V","FLAG"] ...
            + "=" + ["0."+suffix, flagnum]]; % input constants
        k.macros(end+1) = "QUPS_PRECISION="+prc;
        % k.opts = ["-cl-fp32-correctly-rounded-divide-sqrt", "-cl-mad-enable"];
    end
    
    % kernel sizing
    K = d.MaxGridSize(2);
    L = min(k.MaxThreadsPerBlock, ceil(N / K));
    k.ThreadBlockSize = [min(I*M,floor(k.MaxThreadsPerBlock / L)), L, 1]; % I and M are indexed together
    k.GridSize = ceil([I*M, min(N,K*L), N/(K*L)] ./ k.ThreadBlockSize);
    
    % sample
    y_ = k.feval(y_, x_, t_, flagnum); % compute

    % restore type
    if isreal(x0), y_ = real(y_); end
    switch suffix, case "h", y.val = y; otherwise, y = y_; end
else
    % get new dimension mapping
    [~, tmp] = cellfun(@(x) ismember(x, ord), {dim, mdms, rdmst, rdmsx}, 'UniformOutput',false);
    [Tdim, Ndim, Mdim, Fdim] = deal(tmp{:});

    % promote half types
    if isftype(x, 'halfT'), [xc, tc] = dealfun(@single, x, t); % promote until MATLAB's half type is native
    else [xc, tc] = deal(x,t); % keep original otherwise
    end
    
    % compute using interp1, performing for each matching dimension
    pdims = [Tdim, Mdim, Fdim]; % pack all except matching dims
    [xc, tc] = deal(num2cell(xc, pdims), num2cell(tc, pdims));
    parfor(i = 1:numel(xc), 0), y{i} = interp1(xc{i},1+tc{i},interp,extrapval); end

    % fold the non-singleton dimensions of x back down into the singleton
    % dimensions of the output
    if isempty(Mdim), Dt = 1; else, Dt = max(Mdim(size(t,Mdim) ~= 1)); end % max non-singleton dim in t
    % if isempty(Fdim), Dx = 1; else, Dx = max(Fdim(size(x,Fdim) ~= 1)); end % max non-singleton dim in x
    nsing = 1+find(size(t,2:maxdims) == 1); % where t singular in dims of x
    % swap out entries of t that are singular corresponding to where x is non-singular
    y = cellfun(@(y) {swapdim(y, nsing, Dt+nsing-1)}, y);
    y = cat(maxdims+1, y{:}); % unpack Ndims in upper dimensions
    if ~isempty(Ndim), lsz = size(t,Ndim); else, lsz = []; end % forward empty
    y = reshape(y, [size(y,1:maxdims), lsz]); % restore data size in upper dimension
    y = swapdim(y,Ndim,maxdims+(1:numel(Ndim))); % fold upper dimensions back into original dimensions
    
    % return to original type (if it was promoted)
    if isftype(x, 'halfT'), y = halfT(y); end     
    
end

% place back in prior dimensions
if ~isequal(ord, 1:length(ord)) % avoid data copy if already ordered
    y = ipermute(y, ord);
end

function sz = esize(x, sz), if ~isempty(sz), sz = size(x, sz); end
