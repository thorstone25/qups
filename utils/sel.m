function [y, ind] = sel(x, ind, dim, kwargs)
% SEL Subscripting expression to index in one dimension only.
%
% y = SEL(x, ind) selects x at indices ind in the first non-singleton 
% dimension of ind. This allows selection expressions in one dimension to
% be broadcast to all other dimensions.
%
% If any value of ind is NaN, the corresponding value of y will be missing.
%
% y = SEL(x, ind, dim) selects x in dimension dim instead of the first
% non-singleton dimension
%
% [y, ind] = SEL(...) returns the linear indices
% 
% Example:
% x = rand(2,3,4,5);
% [y, ix] = max(x, [], 3);
% [y(1), ix(1)] = deal(NaN); % remove
% ys = sel(x, ix, 3);
% isequaln(y, ys)
% 
% See also SUB
arguments
    x
    ind {mustBeIndex}
    dim (1,1) double {mustBeInteger, mustBePositive} = find(size(ind,1:max(ndims(x), ndims(ind)+1))==1, 1, 'first')
    kwargs.legacy (1,1) logical = false
end

% create index axes for each dimension from ind
sz = size(x); % data size
sz(dim) = size(ind,dim); % use dimensions of the indexing
ONE = ones('like', ind); % match type
if isgpuarray(x), ONE = gpuArray(ONE); end % cast to GPU
i = arrayfun(@(I,d) {reshape((ONE):I, [ones(1,d-1) I 1])}, sz, 1:numel(sz)); % indices per dim - use int32 to save data?

if kwargs.legacy
% expand grid to all dimensions in x (except dim)
[i{:}] = ndgrid(i{:}); % implies [a, b, ...] = ndgrid(a, b, ...);

% broadcast ind to all dimensions in x (except dim)
D = max([ndims(x), ndims(ind)+1, dim]); % maximum dimension
rsz = size(x,1:D) ./ max(1,size(ind,1:D)); % broadcasting size
rsz(dim) = 1; % don't repeat in selected dimension
ind = repmat(ind, rsz); % make full

% get linear indices
ind =  sub2ind(size(x), i{1:dim-1}, ind, i{dim+1:end}); % linear index

else
% get linear indices
ind = vsub2ind(size(x), dim, i{1:dim-1}, ind, i{dim+1:end}); % linear index

end

% select data at ind in dimension dim
val = ~isnan(ind) & logical(ind); % invalid indices 
y(val) = x(ind(val)); % index the valid indices
if any(~val,'all'), y(~val) = missing; end % propagate missing
y = reshape(y, size(ind));

    
function mustBeIndex(ind)
arguments
    ind {mustBeNumeric}
end
j = ~isnan(ind);
mustBeNonnegative(ind(j));
mustBeInteger(ind(j));
    
% local sub2ind - overload types, broadcast, and avoid error-checking
function [ndx, v] = vsub2ind(siz, dim, varargin)

% sizing
K = numel(varargin); % number of indices / dimensions
siz(end+1:K) = 1; % enforce at least K dims 
k = [1 cumprod(siz)]; % stride

% explicit broadcast (implicit pre-allocation)
rsz = siz(1:K) ./ max(1,size(varargin{dim},1:K)); % replication size
rsz(dim) = 1; % exempt indexing dimension
varargin{dim} = repmat(varargin{dim}, rsz); % make full size

% accumulate for each dimension
[ndx, v] = deal(1, true); % init
for i = [dim, 1:dim-1, dim+1:K] % start with input dimension
    ndx = ndx + (varargin{i} - 1) .* k(i); % index
    v   = v   & (1 <= varargin{i} & varargin{i} <= siz(i)); % validity
end
ndx(~v) = 0; % set invalid

