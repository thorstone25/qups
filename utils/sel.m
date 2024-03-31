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
% Example:
% x = rand(2,3,4);
% [y, ix] = max(x, [], 3);
% ys = sel(x, ix, 3);
% all(y == ys, 'all')
% 
% See also SUB
function y = sel(x, ind, dim)
arguments
    x
    ind {mustBeIndex}
    dim (1,1) double {mustBeInteger, mustBePositive} = find(size(ind,1:max(ndims(x), ndims(ind)+1))==1, 1, 'first')
end

% create index axes for each dimension from ind
i = arrayfun(@(I) {int32(1):I}, size(x)); % use int32 to save data?
i{dim} = int32(1):size(ind, dim); % use dimensions of the indexing

% expand grid to all dimensions in x (except dim)
[i{:}] = ndgrid(i{:}); % implies [a, b, ...] = ndgrid(a, b, ...);

% broadcast ind to all dimensions in x (except dim)
D = max([ndims(x), ndims(ind)+1, dim]); % maximum dimension
rsz = size(x,1:D) ./ size(ind,1:D); % broadcasting size
rsz(dim) = 1; % don't repeat in selected dimension
ind = repmat(ind, rsz); % make full

% select data at ind in dimension dim
lind = sub2ind(size(x), i{1:dim-1}, ind, i{dim+1:end}); % linear index
val = ~isnan(lind); % invalid indices 
y(val) = x(lind(val)); % index good indices
if any(~val,'all'), y(~val) = missing; end % propagate missing
y = reshape(y, size(lind));


function mustBeIndex(ind)
arguments
    ind {mustBeInteger, mustBeNonnegative}
end

if isnumeric(ind)
    mustBePositive(ind);
end
    
