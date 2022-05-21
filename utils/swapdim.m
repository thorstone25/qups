function [x, ord] = swapdim(x, i, o)
% SWAPDIM Swap dimensions
%
% y = SWAPDIM(x, i, o) swaps dimensions i and o in the array x. 
%
% Example:
%   x = rand([1,2,3,4,5]);
%   y = swapdim(x,[2 3],[5 4]); % swap 2<->5, 3<->4
%   assert(isequal(size(y), [1,5,4,3,2]))
%
% See also PERMUTE SUB

% make sure we swap the right dimensions
assert(length(i) == length(o), 'Dimensions to swap must be the same length.');
ord = 1:max([ndims(x), i, o]); % all dimensions we have to worry about
ord(i) = o; % swap dimension 
ord(o) = i; 
x = permute(x, ord); 