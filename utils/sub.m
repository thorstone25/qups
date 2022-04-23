% SUB Subscripting expression to slice in one dimension only.
%
% y = SUB(x, ind) returns x(ind,:,:,...,:) where x is any subscriptable
% object.
%
% y = SUB(x, ind, dim) slices x in dimension dim instead of dimension 1
%
% y = SUB(x, ind, vecdim) slices x with the same indices in each dimension 
% in vecdim i.e. returns x(:,...,:,ind,:,...,:,ind,:,...,:)
%
% Example:
% x = rand(2,3,4);
% all(x(:,2:3,:) == sub(x,2:3,2), 'all') 
% 
% 
function y = sub(x, ind, dim)
arguments
    x
    ind (1,:) {mustBeIndex}
    dim (1,:) double {mustBeInteger, mustBePositive} = 1
end
subs = cellstr(repmat(":", [1, max(ndims(x), max(dim))]));
[subs{dim}] = deal(ind);
y = subsref(x, substruct('()', subs));


function mustBeIndex(ind)
arguments
    ind (1,:) {mustBeInteger, mustBeNonnegative}
end

if isnumeric(ind)
    mustBePositive(ind);
end
    
