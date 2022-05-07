function [x, ord] = swapdim(x, i, o)
% SWAPDIM Swap dimensions
%
% y = SWAPDIM(x, i, o) swaps dimensions i and o in the array x. 
%
% Example:
%   x = rand([1,2,3,4,5]);
%   y = swapdim(x,2,4);
%   assert(all(size(y) == [1,4,3,2,5]))

assert(length(i) == length(o), 'Dimensions to swap must be the same length.');
ord = 1:max([ndims(x), i, o]);
ord(i) = o;
ord(o) = i;
x = permute(x, ord);