function z = tenmul(y, x, i, m)
% TENMUL - Tensor multiplication
%
% z = TENMUL(y, x, i, m) multiplies the tensors x and y using dimensions i
% as the inner dimensions (inner-product), m as the matching dimensions
% to be point-wise multiplied(hadamard), and the rest of the dimensions are
% assumed to be the outer dimensions (outer-product). 
% 
% The inner-product dimension is multiplied and summed via matrix
% multiplication. The output is singleton in this dimension. All outer 
% dimensions must be singelton in either x or y.
%
% Example:
%
% x = rand([4, 3, 1, 1, 5]); 
% y = rand([4, 1, 2, 1, 5]); 
% z = tenmul(y, x, 1, 5); % matrix multiply in dim 1
% assert(isequal(size(z), [1, 3, 2, 1, 5]))
%
% See also PAGEMTIMES PERMUTE SWAPDIM

% TODO: make sure i, m are row vectors within dims of x, y
% TODO: handle empty cases better i.e. i or m is empty.

% all dims
d = max(ndims(x), ndims(y));

% hack: if i or m is empty, force it to a singleton dimension
if nargin < 4 || isempty(m), m = d + 1; d = d + 1; end % expand up one      dimension
if nargin < 3 || isempty(i), i = m + 1; d = d + 1; end % expand up one more dimension

% get outer dims
o = setdiff(1:max(d), [i, m]);

% get whether outer in x or y
ox = o(size(x, o) ~= 1);
oy = o(size(y, o) ~= 1);

% move singleton dimensions of o to m
m = union(m, setdiff(o, union(ox, oy)));
o = union(ox, oy);

% check the sizing
% inner dimensions must match
assert(all(size(x,i) == size(y,i)), "Inner dimensions do not match (" + sprintf('%i,', size(x,i)) + "), (" + sprintf('%i,', size(y,i)) + ")")
% matching dimensions must match
assert(all(size(x,m) == size(y,m)), "Outer dimensions do not match (" + sprintf('%i,', size(x,m)) + "), (" + sprintf('%i,', size(y,m)) + ")")
% outer dimensions must be singleton in x or y
assert(isempty(intersect(ox, oy)),  "Matching dimensions are not broadcasting (" + sprintf('%i,', intersect(ox, oy)) + ")");
% should be no dimensions unaccounted for
assert(isempty(setdiff(o, union(ox, oy))),  "Unexpected dimensions: found dimension(s) (" + sprintf('%i,',setdiff(o, union(ox, oy))) + ")");

% get the exact sizing
[szI, szOX, szOY, szM] = deal(size(x, i), size(x, ox), size(y, oy), size(x, m));
[I, OX, OY, M] = dealfun(@prod, szI, szOX, szOY, szM);

% reorder dimensions:
% inner, {outer-x | outer-y}, matching
xr = permute(x, [i, o, m]);
yr = permute(y, [i, o, m]);

% resize to vectorize in these dimension
xr = reshape(xr, [I, OX, M]);
yr = reshape(yr, [I, OY, M]);

% matrix multiply the inner dimension for each matching dimension
% result has the size of the outer product. 
% Y is transposed without conjugate
% TODO: allow conjugate transpose for more numerical accuracy
z = pagemtimes(yr, 'transpose', xr, 'none'); % (OY x OX x M x 1s)

% resize to non-vectorized dimensions
z = reshape(z, [szOY, szOX, szM]); % (OY x OX x M [x 1s])

% reorder to original broadcasting dims, placing collapsed (singleton) dims 
% at the end
z = ipermute(z, [oy, ox, m, i]); 

function varargout = dealfun(fun, varargin), varargout = cellfun(fun, varargin, 'UniformOutput', false);