function I = argmax(varargin), [~, I] = max(varargin{:}); end
% ARGMAX - Index of the maximum elements of an array
%
% I = ARGMAX(X) returns the indices corresponding to the maximum
% values. The values in I index into the dimension of X that is being
% operated on. If X contains more than one element with the maximum
% value, then the index of the first one is returned.
%
% I = ARGMAX(X,[],'all') also returns the linear index into X that 
% corresponds to the maximum value over all elements in X.
%
% I = ARGMAX(X,[],DIM) operates along the dimension DIM.
%
% ARGMAX calls max to perform the calculation. Other inputs accepted by max
% are valid here.
%
% Example:
% x = rand([1 10]);
% [y, i] = max(x);
% j = argmax(x);
% 
% assert(isequal(i, j));
% assert(isequal(y, x(j)));
% 
% See also MAX ARGMIN
