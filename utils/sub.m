% SUB Subscripting expression to slice in one dimension only.
%
% y = SUB(x, ind) returns x(ind,:,:,...,:) where x is any subscriptable
% object.
%
% y = SUB(x, ind, dim) slices x in dimension dim instead of dimension 1
%
% y = SUB(x, vecind, vecdim) slices x indices in each specified dimension 
% in vecdim i.e. returns x(:,...,:,ind,:,...,:,ind,:,...,:). vecind must be
% a cell array of indices. 
% 
% Example:
% x = rand(2,3,4,5,4,3,2);
% assert(isequal(x(2,:,:,:,:,:,:), sub(x,2,1))); 
% assert(isequal(x(:,2:3,4,:,:,:,:), sub(x,{2:3,4},[2 3])))
%
% 

function y = sub(x, ind, dim)
% default to dim 1. TODO: default to 1st non-singleton
if nargin < 3, dim = 1; end 

% ensure indices placed in a cell
if ~iscell(ind), ind = {ind}; end 

% check vecind/vecdim correspond
assert(numel(ind) == numel(dim), "Must have as many dims as indices.") 

% place ':' in all dimensions
subs = cellstr(repmat(":", [1, max(gather(ndims(x)), max(dim))])); 

% replace ':' with index in selected dimensions
for i = 1:numel(dim), subs{dim(i)} = ind{i}; end

% index the object
y = subsref(x, substruct('()', subs));
