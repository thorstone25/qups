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
ord(i) = o; % swap dimensions
ord(o) = i; 

% if the data is already ordered, do nothing
if isequal(ord, 1:length(ord)), return, end

% if at most 1 dim between i and o is non-scalar ...
iosz = size(x, [(i:o), (o:i)]); % sizing from i<->o inclusive
if sum(iosz ~= 1) <= 1
    % we can reshape to save a data copy operation
    sz = size(x, 1:length(ord));
    [sz(i), sz(o)] = deal(sz(o), sz(i)); % swap sizing
    x = reshape(x, sz); % set the new size
else % implement a generalized transpose
    x = permute(x, ord);
end