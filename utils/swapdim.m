function [x, ord] = swapdim(x, i, o)
% SWAPDIM Swap dimensions
%
% y = SWAPDIM(x, i, o) swaps dimensions i and o in the array x. 
%
% y = SWAPDIM(x, i) swaps dimensions i with dimensions 1:numel(i). This is
% equivalent to calling y = permute(x, i); 
% 
% [y, ord] = SWAPDIM(...) also returns the permution order ord.
% 
% SWAPDIM optimistically calls reshape when feasible to avoid an implicit
% copy and falls back to permute when required.
% 
% Example:
%   x = rand([1,2,3,4,5]);
%   [y, ord] = swapdim(x,[2 3],[5 4]); % swap 2<->5, 3<->4
%   assert(isequal(size(y), [1,5,4,3,2]))
%   assert(isequal(y, permute(x, ord)));
%
% See also SUB PERMUTE RESHAPE
arguments
    x
    i (1,:) {mustBeInteger}
    o (1,:) {mustBeInteger} = 1:numel(i)
end

% make sure we swap the right dimensions
assert(length(i) == length(o), 'Dimensions to swap must be the same length.');
L = max([ndims(x), i, o]);
ord = 1:L; % all dimensions we have to worry about
ord0 = ord;

% get permutation
if isempty(intersect(i, o)) % can be handled separately
    ord([i o]) = [o i]; % swap dimensions
else % must fill in missing indices
    l = min(min(i),min(o)) : max(max(i),max(o)); % all indices within swap
    i = [i, setdiff(l, i)]; % expanded input indices
    o = [o, setdiff(l, o)]; % expanded output indices
    ord(o) = i; % full permutation ordering
end

% if the data is already ordered, do nothing
if isequal(ord, ord0), return, end

% get the non-singleton output dimensions, ordered
ons = ord(size(x, ord) > 1); 

% if the output (and input) dimensions are ordered
if issorted(ons, 'ascend')  % we can reshape to save a data copy operation
    x = reshape(x, size(x, ord)); % set the new size
else % otherwise implement a generalized transpose
    x = permute(x, ord);
end

