% SUB Subscripting expression to slice in one dimension only.
%
% y = SUB(x, ind) returns x(ind,:,:,...,:) where x is any subscriptable
% object.
%
% y = SUB(x, ind, dim) slices x in dimension dim instead of dimension 1
%
% y = SUB(x, vecind, vecdim) slices x indices in each specified dimension 
% in vecdim i.e. returns x(:,...,:,ind,:,...,:,ind,:,...,:). 
% 
% Example:
% x = rand(2,3,4);
% all(x(:,2:3,:) == sub(x,2:3,2), 'all') 
% 

function y = sub(x, ind, dim)
if nargin < 3, dim = 1; end
subs = cellstr(repmat(":", [1, max(ndims(x), max(dim))]));
if ~iscell(ind), ind = {ind}; end
assert(numel(ind) == numel(dim), "Must have as many dims as indices.")
for i = 1:numel(ind)
    subs{dim(i)} = ind{i};
end
y = subsref(x, substruct('()', subs));
