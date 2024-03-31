% SUB Subscripting expression to slice in one dimension only.
%
% y = SUB(x, ind) returns x(ind,:,:,...,:) where x is any subscriptable
% object.
%
% y = SUB(x, ind, dim) slices x in dimension dim. The default is the first
% non-singleton dimension.
%
% y = SUB(x, vecind, vecdim) slices x indices in each specified dimension 
% in vecdim i.e. returns x(:,...,:,ind,:,...,:,ind,:,...,:). vecind must be
% a cell array of indices.
%
% i = SUB(..., true) returns the indexing expression used to index the
% object. This is useful to write statements using subsasgn and subsref.
% 
% Example:
% x = rand(2,3,4,5,4,3,2);
% i = sub(x,2,1,true), % the indexing expression
% assert(isequal(x(2,:,:,:,:,:,:), sub(x,2)));
% assert(isequal(x(2,:,:,:,:,:,:), subsref(x, i)));
% assert(isequal(x(:,2:3,:,4,:,:,[true,false]), ...
%            sub(x,{ 2:3,  4,    [true,false]}, [2 4 7])))
% 
% 
% See also SEL SUBSREF SUBSASGN SUBSTRUCT

function y = sub(x, ind, dim, expr)
if nargin < 3, dim = max([1, find(size(x) ~= 1, 1, 'first')]); end
if nargin < 4, expr = false; end

% ensure indices placed in a cell
if ~iscell(ind), ind = {ind}; end 

% check vecind/vecdim correspond
assert(numel(ind) == numel(dim), "Must have as many dims as indices.") 

% short-circuit: empty inds/dims -> identity
if isempty(ind) && isempty(dim), dim = 0; end

% place ':' in all dimensions
subs = cellstr(repmat(":", [1, max(gather(ndims(x)), max(dim))])); 

% replace ':' with index in selected dimensions
for i = 1:numel(ind), subs{dim(i)} = ind{i}; end

% get object indexing expression
iy = substruct('()', subs);

% return the indexed object, or the expression
if expr, y = iy; else, y = subsref(x, iy); end
