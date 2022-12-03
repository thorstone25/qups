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
arguments
    x
    i (1,:) {mustBeInteger}
    o (1,:) {mustBeInteger}
end

% make sure we swap the right dimensions
assert(length(i) == length(o), 'Dimensions to swap must be the same length.');
L = max([ndims(x), i, o]);
ord = 1:L; % all dimensions we have to worry about
ord(i) = o; % swap dimensions
ord(o) = i; 

% if the data is already ordered, do nothing
if isequal(ord, 1:numel(ord)), return, end

% check the sizing to see if we can reshape
e = [1 : min([i,o]) - 1, max([i,o]) + 1 : L]; % external indices - not involved
m = setdiff(1:L, [i, o, e]); % middle indices - involved, not moving
msz = prod(esize(x,m)); % middle indices total size
isz = prod(esize(x,i)); % input  indices total size
osz = prod(esize(x,o)); % output indices total size

% if at most 1 dim between i, m, and o is non-scalar 
% and the order of dims in i matches the order in o
if sum([isz msz osz] ~= 1) <= 1 && ...
    isequal(argsort(i), argsort(o))
    
    % we can reshape to save a data copy operation
    sz = size(x, 1:numel(ord));
    [sz(i), sz(o)] = deal(sz(o), sz(i)); % swap sizing
    x = reshape(x, sz); % set the new size
else % implement a generalized transpose
    x = permute(x, ord);
end

function n = esize(x, i), if isempty(i), n = []; else; n = size(x, i); end

function o = argsort(x), [~, o] = sort(x); 