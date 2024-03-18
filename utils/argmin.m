function I = argmin(varargin), [~, I] = min(varargin{:}); end
% ARGMIN - Index of the minimum elements of an array
%
% I = ARGMIN(X) returns the indices corresponding to the minimum
% values. The values in I index into the dimension of X that is being
% operated on. If X contains more than one element with the minimum
% value, then the index of the first one is returned.
%
% I = ARGMIN(X,[],'all') also returns the linear index into X that 
% corresponds to the minimum value over all elements in X.
%
% I = ARGMIN(X,[],DIM) operates along the dimension DIM.
%
% ARGMIN calls min to perform the calculation. Other inputs accepted by min
% are valid here.
%
% Example:
% x = rand([1 10]);
% [y, i] = min(x);
% j = argmin(x);
% 
% assert(isequal(i, j));
% assert(isequal(y, x(j)));
% 
% See also MIN ARGMAX