function y = wsinterpd(x, t, dim, w, sdim, interp, extrapval, varargin)
% WSINTERPD GPU-enabled interpolation in one dimension
%
% y = WSINTERPD(x, t) interpolates the data x at the indices t. It uses the
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
% y = WSINTERPD(x, t, dim) samples the data in dimension dim instead of
% dimension 1
%
% y = WSINTERPD(x, t, dim, w) applies the weighting array via point-wise 
% multiplication after sampling the the data. The dimensions must be 
% compatible with the sampling array t in the sampling dimension dim. 
% The default is 1.
%
% y = WSINTERPD(x, t, dim, w, sdim) sums the data in the dimension(s) sdim
% after the weighting matrix has been applied. The default is [].
%
% y = WSINTERPD(x, t, dim, w, sdim, interp) specifies the interpolation 
% method. It must be one of {"nearest", "linear"*,"cubic","lanczos3"} or 
% any option supported by interp1.
%
% y = WSINTERPD(x, t, dim, w, sdim, interp, extrapval) uses extrapval as the
% extrapolation value when outside of the domain of t.
%
% See also INTERPN INTERP1 INTERPD
%

%% validate dimensions
if nargin < 7 || isempty(extrapval), extrapval = nan; end
if nargin < 6 || isempty(interp),    interp = 'linear'; end
if nargin < 5 || isempty(sdim),      sdim = []; end 
if nargin < 4 || isempty(w        ), w = 1; end
if nargin < 3 || isempty(dim),       dim = 1; end
assert(isreal(t), 'Sample indices must be real.'); % sampling at a real data type
assert(isscalar(extrapval), 'Only a scalar value accepted for extrapolation.');

maxdims = max([ndims(x), ndims(t), ndims(w)]); % maximum dimension
ddms = setdiff(1:maxdims, dim); % data dimensions (all except time/image)
xsz = size(x, ddms);
tsz = size(t, ddms);
mdms  = ddms(xsz == tsz); % matching dimensions
rdmsx = ddms(xsz ~= tsz & tsz == 1); % outer dim for x
rdmst = ddms(xsz ~= tsz & xsz == 1); % outer dim for t
rdms = union(rdmsx, rdmst); % all outer dimensions

% remove summation in singleton dimensions
if ~isempty(sdim), sdim(size(x,sdim) == 1 & size(t,sdim) == 1) = []; end

% set the data order: sampling dimensions, broadcasting dimensions,
% replicating dimensions (in t, then in x)
ord = [dim, mdms, rdmst, rdmsx]; % this should be all of the dimensions
assert(isempty(setxor(ord, 1:maxdims)), 'Unable to identify all dimensions: this is a bug.');
assert(all(...
    size(w,1:maxdims) == 1 | ...
    size(w,1:maxdims) == max(size(x,1:maxdims), size(t,1:maxdims)) ...
    ), 'The weighting vector w must have dimensions compatible with the data and may not broadcast.' ...
    );

% get the data sizing
T = size(x, dim);
I = size(t, dim);
N = prod(esize(t, mdms));
M = prod(esize(t, rdms));
F = prod(esize(x, rdms));

% move the input data to the proper dimensions for the GPU kernel
x = permute(x, ord);
t = permute(t, ord);
w = permute(w, ord);

S = maxdims; % number of (full) dimensions
dsizes = [I, max(size(t,2:S), size(x,2:S))]; % data sizes in new dimensions
sdimo = arrayfun(@(d)find(d==ord), sdim); % new dimensions of summation

% get stride for weighting 
wstride = size(w,1:S);
wstride = cumprod([1, wstride(1:end-1)]) .* (wstride ~= 1);

% get stride for output
osz = dsizes; osz(sdimo) = 1; % these dimensions are summed
ystride = cumprod([1, osz(1:end-1)]) .* (osz ~= 1);

% function to determine type
isftype = @(x,T) strcmp(class(x), T) || any(arrayfun(@(c)isa(x,c),["tall", "gpuArray"])) && strcmp(classUnderlying(x), T);

% use ptx on gpu if available or use native MATLAB
if exist('interpd.ptx', 'file') ...
        && ( isa(x, 'gpuArray') || isa(t, 'gpuArray') || isa(x, 'halfT') && x.gtype) ...
        && (ismember(interp, ["nearest", "linear", "cubic", "lanczos3"]))
    
    % determine the data type
    if     isftype(x, 'double')
        suffix = "" ; [x,t,w] = dealfun(@double, x, t, w);
    elseif isftype(x, 'single')
        suffix = "f"; [x,t,w] = dealfun(@single, x, t, w);
    elseif isftype(x, 'halfT'  )
        suffix = "h"; [x,t,w] = dealfun(@(x)gpuArray(halfT(x)), x, t, w); % custom type
    else
        warning("Datatype " + class(x) + " not recognized as a GPU compatible type.");
        suffix = "f" ;
    end

    % grab the kernel reference
    k = parallel.gpu.CUDAKernel('interpd.ptx', 'interpd.cu', 'wsinterpd' + suffix); 
    k.setConstantMemory( ...
        'QUPS_I', uint64(I), 'QUPS_T', uint64(T), 'QUPS_S', uint64(S), ...
        'QUPS_N', uint64(N), 'QUPS_M', uint64(M), 'QUPS_F', uint64(F) ...
        );
    k.ThreadBlockSize = min(k.MaxThreadsPerBlock,I*M); % I and M are indexed together
    k.GridSize = [ceil(I*M ./ k.ThreadBlockSize(1)), N, 1];

    % translate the interp flag
    switch interp
        case "nearest", flagnum = 0;
        case "linear",  flagnum = 1;
        case "cubic",   flagnum = 2;
        case "lanczos3",flagnum = 3;
        otherwise, error('Interp option not recognized: ' + string(interp));
    end

    % sample
    [x,w] = deal(complex(x), complex(w)); % enforce complex type
    switch suffix
        case "h", 
            y = complex(gpuArray(halfT(zeros(osz))));
            [w_,x_,t_,y_] = dealfun(@(x)x.val, w,x,t,y); % extract data

        otherwise, 
            [w_,x_,t_] = deal(w,x,t); % copy data
            y_ = cast(zeros(osz), 'like', x_);
    end
     % zeros: uint16(0) == storedInteger(half(0)), so this is okay
    y_ = k.feval(y_, w_, x_, t_, dsizes, wstride, ystride, flagnum); % compute

    % for halfT, store the data back in the output
    switch suffix, case "h", y.val = y_; otherwise, y = y_; end
else
    % get new dimension mapping
    [~, tmp] = cellfun(@(x) ismember(x, ord), {dim, mdms, rdmst, rdmsx}, 'UniformOutput',false);
    [Tdim, Ndim, Mdim, Fdim] = deal(tmp{:});

    % promote half types
    if isftype(x, 'halfT'), [xc, tc] = dealfun(@single, gather(x), t); % promote until MATLAB's half type is native
    else [xc, tc] = deal(x,t); % keep original otherwise
    end
    
    % get the dimensions after using interp1
    if isempty(Mdim), Dt = 1; else, Dt = max(Mdim(size(t,Mdim) ~= 1)); end % max non-singleton dim in t
    % if isempty(Fdim), Dx = 1; else, Dx = max(Fdim(size(x,Fdim) ~= 1)); end % max non-singleton dim in x

    % identify outer dimensions for summing
    dsume = 1+ndims(xc)+ndims(tc); % add a scalar dimension - avoid sum over empty dims
    dsume = [dsume,      sdimo(~ismember(sdimo, Fdim))]; % add dimensions that aren't shifted high (external)
    dsumi = [dsume, Dt-1+sdimo( ismember(sdimo, Fdim))]; % add dimensions that are    shifted high (internal only)
    w = swapdim(w, Fdim, Dt-1+Fdim); % move weights into shifted upper dimensions

    % on CPU, force w to be full for matching dimensions by adding zeros
    % TODO: remove this by making a hybrid approach: instead of saving to
    % cells, expand the NDarray and add atomically
    for d = Ndim, w = w + shiftdim(zeros([size(x,d),1]),1-d); end
    
    % compute using interp1, performing for each matching dimension,
    % weighting and summing over outer dimensions as requested
    pdims = unique([Tdim, Mdim, Fdim, Dt-1+Fdim]); % pack all except matching dims
    [xc, tc, wc] = deal(num2cell(xc, pdims), num2cell(tc, pdims), num2cell(w, pdims));
    parfor(i = 1:numel(xc), 0), y{i} = sum(wc{i} .* interp1(xc{i},1+tc{i},interp,extrapval), dsumi, 'omitnan'); end

    % fold the non-singleton dimensions of x back down into the singleton
    % dimensions of the output
    % swap out entries of t that are singular corresponding to where x is non-singular
    nsing = 1+find(size(t,2:maxdims) == 1); % where t singular in dims of x
    y = cellfun(@(y) {swapdim(y, nsing, Dt+nsing-1)}, y);
    y = cat(maxdims+1, y{:}); % unpack Ndims in upper dimensions
    if ~isempty(Ndim), lsz = size(t,Ndim); else, lsz = []; end % forward empty
    y = reshape(y, [size(y,1:maxdims), lsz]); % restore data size in upper dimension
    y = swapdim(y,Ndim,maxdims+(1:numel(Ndim))); % fold upper dimensions back into original dimensions
    y = sum(y, dsume, 'omitnan'); % sum across matching dimensions too
    
    % demote for half types
    if isftype(x, 'halfT'), y = halfT(y); end % convert half types back
end

% place back in prior dimensions
y = ipermute(y, ord);


function d = esize(varargin), if isempty(varargin{2}), d = []; else, d = size(varargin{:}); end
% ESIZE - Empty-forwarding size
%
% d = ESIZE(x) or d = ESIZE(x,dim) returns [] if x is empty or returns
% size(x,dim) otherwise.
