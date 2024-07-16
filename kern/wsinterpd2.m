function y = wsinterpd2(x, t1, t2, dim, w, sdim, interp, extrapval, omega)
% WSINTERPD2 GPU-enabled weight-and-sum interpolation in one dimension with
% 2 delay inputs
%
% y = WSINTERPD2(x, t1, t2) interpolates the data x at the indices 
% t = t1 + t2. It uses the matching and non-matching dimensions of the data
% to implicitly broadcast operations across the data. In other words, the
% data is formatted into the following:
%
%   x is (T x N' x      F')
%   t is (I x N' x M'     )
%   y is (I x N' x M' x F')
%
% Defining t1 and t2 separately allows for a lower memory footprint.
% 
% Sampling is performed for all I, element-wise in N' (where matched) and 
% broadcasted across M' and F' such that x(k,n',m',f') = x(k,n',f') and
% t(i,n',m',f') = t(i,n',m') 
%
% x is 0-base indexed implicitly i.e. x(1) in MATLAB <-> t == 0
%
% y = WSINTERPD2(x, t1, t2, dim) samples the data in dimension dim instead of
% dimension 1
%
% y = WSINTERPD2(x, t1, t2, dim, w) applies the weighting array via point-wise 
% multiplication after sampling the the data. The dimensions must be 
% compatible with the sampling array t in the sampling dimension dim. 
% The default is 1.
%
% y = WSINTERPD2(x, t1, t2, dim, w, sdim) sums the data in the dimension(s) sdim
% after the weighting matrix has been applied. The default is [].
%
% y = WSINTERPD2(x, t1, t2, dim, w, sdim, interp) specifies the interpolation 
% method. It must be one of {"nearest", "linear"*,"cubic","lanczos3"} or 
% any option supported by interp1.
%
% y = WSINTERPD2(x, t1, t2, dim, w, sdim, interp, extrapval) uses extrapval as the
% extrapolation value when outside of the domain of t.
%
% y = WSINTERPD2(x, t1, t2, dim, w, sdim, interp, extrapval, omega) applies
% an exponential weighting `exp(omega .* (t1 + t2))` after sampling. The
% default is 0.
%
% This can be useful when omega is a complex value representing a phasor
% rotation.
%
% See also WSINTERPD INTERPN INTERP1 INTERPD
%

%% validate dimensions
persistent k prc0; % cache prior kernel calls
if nargin < 9 || isempty(omega),     omega = 0; end
if nargin < 8 || isempty(extrapval), extrapval = nan; end
if nargin < 7 || isempty(interp),    interp = 'linear'; end
if nargin < 6 || isempty(sdim),      sdim = []; end 
if nargin < 5 || isempty(w        ), w = 1; end
if nargin < 4 || isempty(dim),       dim = 1; end
assert(isreal(t1) && isreal(t2), 'Sample indices must be real.'); % sampling at a real data type
assert(isscalar(extrapval), 'Only a scalar value accepted for extrapolation.');

% first: move dim to 1st dimension, as it's always required
x  = swapdim(x ,dim,1);
t1 = swapdim(t1,dim,1);
t2 = swapdim(t2,dim,1);
w  = swapdim(w ,dim,1);
sdim_ = sdim; % store
[sdim(sdim_ == 1), sdim(sdim_ == dim)] = deal(dim, 1); % swap
[ogdim, dim] = deal(dim, 1); % store original

maxdims = max([ndims(x), ndims(t1), ndims(t2), ndims(w)]); % maximum dimension
ddms = setdiff(1:maxdims, dim); % data dimensions (all except time)
xsz = size(x, ddms);
tsz = max(size(t1, ddms), size(t2,ddms)); % defined by broadcasted addition of t1 and t2
mdms  = ddms(xsz == tsz); % matching dimensions
rdmsx = ddms(xsz ~= tsz & tsz == 1); % outer dim for x
rdmst = ddms(xsz ~= tsz & xsz == 1); % outer dim for t
% rdms = union(rdmsx, rdmst); % all outer dimensions

% remove summation in singleton dimensions
if ~isempty(sdim), sdim(all([size(x,sdim); size(t1,sdim); size(t2,sdim)] == 1,1)) = []; end

% set the data order: sampling dimensions, broadcasting dimensions,
% replicating dimensions (in t, then in x)
alldims = [dim, mdms, rdmst, rdmsx]; % this should be all of the dimensions
assert(isempty(setxor(alldims, 1:maxdims)), 'Unable to identify all dimensions: this is a bug.');
assert(all(...
    size(w,ddms) == 1 | ...
    size(w,ddms) == max(xsz, tsz) ...
    ) & any(size(w,dim) == [1, size(t1,dim), size(t2,dim)]), 'The weighting vector w must have dimensions compatible with the data and may not broadcast.' ...
    );

S = maxdims; % number of (full) dimensions
dsizes = [max(size(t1,1),size(t2,1)), max([size(t1,2:S); size(t2,2:S); size(x,2:S)],[],1)]; % data sizes in new dimensions
% sdimo = arrayfun(@(d)find(d==ord), sdim); % new dimensions of summation

% function to determine type
isftype = @(x,T) strcmp(class(x), T) || any(arrayfun(@(c)isa(x,c),["tall", "gpuArray"])) && strcmp(classUnderlying(x), T);

use_gdev = exist('interpd.ptx', 'file') ...
        && ( isa(x, 'gpuArray') || isa(t1, 'gpuArray') || isa(t2, 'gpuArray') || isa(x, 'halfT') && x.gtype) ...
        && (ismember(interp, ["nearest", "linear", "cubic", "lanczos3"])) ... 
        && real(omega) == 0;

use_odev = exist('interpd.cl', 'file') ...
        && (ismember(interp, ["nearest", "linear", "cubic", "lanczos3"])) ... 
        && real(omega) == 0;

if use_odev && ~use_gdev && exist('oclDeviceCount','file') && oclDeviceCount(), ocl_dev = oclDevice(); else, ocl_dev = []; end % get oclDevice if supported

use_odev = use_odev && ~use_gdev ... 
        && ~isempty(ocl_dev) && ...
        (  (isftype(x,'double') && ocl_dev.SupportsDouble) ...
        || (isftype(x,'single')                          ) ... always supported
        || (isftype(x,'halfT' ) && ocl_dev.SupportsHalf  ));

% use ptx on gpu if available or use native MATLAB
if use_gdev || use_odev

    % get stride for weighting
    wstride = sz2stride(size(w,1:S));

    % get stride for data
    xstride = sz2stride([1 size(x,2:S)]); % treat as singular in dim 1

    % get stride for time
    t1stride = sz2stride(size(t1,1:S));
    t2stride = sz2stride(size(t2,1:S));

    % get stride for output
    osz = dsizes; osz(sdim) = 1; % these dimensions are summed
    ystride = sz2stride(osz);

    % determine the data type
    if     isftype(x, 'double')
        suffix = "" ; prc = 64; [x,t1,t2,w] = dealfun(@double, x, t1, t2, w);
    elseif isftype(x, 'single')
        suffix = "f"; prc = 32; [x,t1,t2,w] = dealfun(@single, x, t1, t2, w);
    elseif isftype(x, 'halfT'  )
        suffix = "h"; prc = 16; [x,t1,t2,w] = dealfun(@(x)gpuArray(halfT(x)), x, t1, t2, w); % custom type
    else
        suffix = "f" ; prc = 32;
        warning("Datatype " + class(x) + " not recognized as a GPU compatible type.");
    end

    % get the data sizing
    T = size(x, dim);
    I = prod(max(esize(t1, [dim, rdmst]), esize(t2, [dim, rdmst])));
    N = prod(esize(x, mdms));
    F = prod(esize(x, rdmsx));

    % translate the interp flag
    switch interp
        case "nearest", flagnum = 0;
        case "linear",  flagnum = 1;
        case "cubic",   flagnum = 2;
        case "lanczos3",flagnum = 3;
        otherwise, error('Interp option not recognized: ' + string(interp));
    end

    % enforce complex type for gpu data
    if use_odev || isa(x, 'gpuArray') || isa(x, 'halfT'), x = complex(x); end
    if use_odev || isa(w, 'gpuArray') || isa(w, 'halfT'), w = complex(w); end
    switch suffix
        case "h", 
            y = complex(gpuArray(halfT(zeros(osz))));
            [w_,x_,t1_,t2_,y_] = dealfun(@(x)x.val, w,x,t1,t2,y); % extract data

        otherwise, 
            [w_,x_,t1_,t2_] = deal(w,x,t1,t2); % copy data
            y_ = zeros(osz, 'like', x_); % pre-allocate output
    end
    % zeros: uint16(0) == storedInteger(half(0)), so this is okay

    % use a cached kernel (reduces likelihood of overloading OpenCL)?
    use_cached = isscalar(k) && isvalid(k) && existsOnGPU(k) && (prc == prc0) ...
            && (( use_gdev && isa(k,'parallel.gpu.CUDAKernel')) ...
            || (  use_odev && isa(k,'oclKernel') && k.Device.Index == ocl_dev.Index) ...
            );
    % if use_cached, warning("Using cached " + class(k)), disp(k), end
    if use_gdev
        % grab the kernel reference
        if ~use_cached
        k = parallel.gpu.CUDAKernel('interpd.ptx', 'interpd.cu', 'wsinterpd2' + suffix);
        end
        k.setConstantMemory( ...
            'QUPS_I', uint64(I), 'QUPS_T', uint64(T), 'QUPS_S', uint64(S), ...
            'QUPS_N', uint64(N), 'QUPS_F', uint64(F) ...
            );
    elseif use_odev
        % get the kernel reference
        if ~use_cached
        k = oclKernel(which('interpd.cl'), 'wsinterpd2');

        % set the data types
        switch prc, case 16, t = 'uint16'; case 32, t = 'single'; case 64, t = 'double'; end
        k.defineTypes({t,t}); % all aliases are this type
        end

        % configure the kernel sizing and options
        k.macros = "QUPS_" + ["I", "T", "S", "N", "F"] + "=" + uint64([I T S N F]); % size constants
        k.macros = [k.macros, "QUPS_INTERPD_" + ["NO_V","FLAG","OMEGA"] ...
            + "=" + ["0."+suffix, flagnum, imag(omega)]]; % input constants
        k.macros(end+1) = "QUPS_PRECISION="+prc;
        % k.opts = ["-cl-fp32-correctly-rounded-divide-sqrt", "-cl-mad-enable"];

        % rebuild for this data size
        k.build();
    end

    % constant arguments
    if use_gdev, d = gpuDevice(); cargs = {flagnum, imag(omega)}; else, d = k.Device; cargs = {}; end 

    % kernel sizing
    K = d.MaxGridSize(2);
    L = max(1,min(k.MaxThreadsPerBlock, ceil(N / K)));
    k.ThreadBlockSize = [min(I,floor(k.MaxThreadsPerBlock / L)), L, 1]; % local group size
    k.GridSize = max(1,ceil([I, min(N,K*L), N/(K*L)] ./ k.ThreadBlockSize));
    [dsz, prc0] = deal([I T S N F], prc); % record the data size / precision

    % index label flags
    iflags = zeros([1 maxdims], 'uint8'); % increase index i
    iflags(mdms) = 1; % increase index n
    iflags(rdmsx) = 2; % increase index f
    
    % strides
    strides = cat(1,wstride,ystride,t1stride,t2stride,xstride);

    % compute
    y_ = k.feval(y_, w_, x_, t1_, t2_, dsizes, iflags, strides, cargs{:}); 

    % for halfT, store the data back in the output
    switch suffix, case "h", y.val = y_; otherwise, y = y_; end
else

    % promote half types
    if isftype(x, 'halfT'), [xc, t1c, t2c] = dealfun(@single, gather(x), t1, t2); % promote until MATLAB's half type is native
    else, [xc, t1c, t2c] = deal(x,t1,t2); % keep original otherwise
    end
    
    % identify cumulative permutation of the non-matching dimensions
    Dt = max(rdmst(find(esize(t1, rdmst) ~= 1 | esize(t2, rdmst) ~= 1,1,'last')), 1); if isempty(Dt), Dt = 1; end
    Dx = max(rdmsx(find(esize(x, rdmsx) ~= 1,1,'last')), 1); if isempty(Dx), Dx = 1; end

    % identify outer dimensions for summing
    dsume = 1+ndims(xc)+max(ndims(t1c), ndims(t2c)); % add a scalar dimension - avoid sum over empty dims
    dsume = [dsume,      sdim(~ismember(sdim, rdmsx))]; % add dimensions that aren't shifted high (external)
    dsumi = [dsume, Dt-1+sdim( ismember(sdim, rdmsx))]; % add dimensions that are    shifted high (internal only)

    % compute using interp1, performing for each matching dimension,
    % weighting and summing over outer dimensions as requested
    pdims = setdiff(1:maxdims, mdms); % pack all except matching dims
    [xc, t1c, t2c, wc] = deal(num2cell(xc, pdims), num2cell(t1c, pdims), num2cell(t2c, pdims), num2cell(w, pdims));

    % get the sub indexing
    assert(all(numel(xc) >= cellfun(@numel, {t1c, t2c, wc})), "Internal error: data size should be larger than temporal and weighting sampling");
    [xsz, t1sz, t2sz, wcsz] = dealfun(@(x)size(x,1:maxdims), xc, t1c, t2c, wc);
    [t1st, t2st, wcst] = dealfun(@sz2stride, t1sz, t2sz, wcsz);
    % wc = repmat(wc, size(xc,1:maxdims) ./ size(wc,1:maxdims)); % ensure full size for indexing TODO: remove this limitation

    % determine the parallel environment
    if any(cellfun(@(x)isa(x,'gpuArray'), {x,t1,t2,w})), parenv = 0; % no parenv for gpuArrays
    elseif isempty(gcp("nocreate")),                     parenv = 0; % no parpool if none active
    elseif isa(gcp("nocreate"), "parallel.ProcessPool"), parenv = 0; % no parpool - memory overhead outweighs benefits
    else,                                                parenv = gcp; % ideally a ThreadPool
    end

    nsing = 1+find((size(t1,2:maxdims) == 1) & (size(t2,2:maxdims) == 1)); % where t singular in dims of x

    parfor(i = 1:numel(xc), parenv)
    % for i = numel(xc):-1:1
        % get 0-based sub indices
        inds = cell(1,maxdims);
        [inds{:}] = ind2sub(xsz, i);

        % get 1-based linear index for each variable
        [jt1, jt2, jwc] = dealfun(@(st) 1 + sum(([inds{:}] - 1) .* st), t1st, t2st, wcst);
        
        % index, sample, weight, sum
        tau = (t1c{jt1} + t2c{jt2}); %#ok<PFBNS> % time
        osz = [esize(tau,1:Dt), esize(xc{i},2:Dx),1,1]; % proper output sizing
        amp = wc{jwc}; %#ok<PFBNS> % weighting
        amp = swapdim(amp, nsing, Dt-1+nsing); % move singleton dimensions into upper dimensions
        yi = sum(exp(omega .* tau) .* amp .* reshape(interp1(xc{i},1+tau,interp,extrapval), osz), dsumi, 'omitnan'); 
        
        % store output
        yi = swapdim(yi, nsing, Dt-1+nsing); % move singleton dimensions into upper dimensions
        y{i} = yi;
    end

    % fold the non-singleton dimensions of x back down into the singleton
    % dimensions of the output
    % swap out entries of t that are singular corresponding to where x is non-singular
    y = cat(maxdims+1, y{:}); % unpack matching dims in upper dimensions
    lsz = max(esize(t1,mdms), esize(t2,mdms)); % get original shape of matching dims
    y = reshape(y, [size(y,1:maxdims), lsz]); % restore data size in upper dimension
    y = swapdim(y,mdms,maxdims+(1:numel(mdms))); % fold upper dimensions back into original dimensions
    y = sum(y, dsume, 'omitnan'); % sum across matching dimensions too
    
    % demote for half types
    if isftype(x, 'halfT'), y = halfT(y); end % convert half types back
end

% restore dimension dim
y = swapdim(y,ogdim,1);


function d = esize(varargin), if nargin >= 2 && isempty(varargin{2}), d = []; else, d = size(varargin{:}); end
% ESIZE - Empty-forwarding size
%
% d = ESIZE(x) or d = ESIZE(x,dim) returns [] if dim is empty or returns
% size(x,dim) otherwise.

function stride = sz2stride(sz), stride = cumprod([1, sz(1:end-1)]) .* (sz ~= 1);
