function [tf, tol] = isalmostn(a, b, tol)
% ISALMOSTN - Are two arrays almost equal (omitnan)
% 
% tf = ISALMOSTN(a, b) computes wheteher two arrays are almost equal
%
% tf = ISALMOSTN(a, b, tol) uses the tolerance tol for comparing the
% difference rather than a hueristic based on the maximum value.
%
% [tf, tol] = ISALMOSTN(...) returns the tolerance used for the difference
%
% See also ISEQUAL ISEQUALN

if nargin < 3, tol = 1e2*eps(max(max(a,[],'all'), max(b,[],'all'))); end

D = max(ndims(a), ndims(b)); % max dimension

if(~all(size(a,1:D) == size(b,1:D)) )
    tf = false;
    return;
end

% find if any invalid
ntf = xor(isnan(a), isnan(b)); % one, but not the other, is NaN
vtf = abs(a - b) >= tol; % outside of tolerance
ttf = ntf | vtf;
tf = ~any(ttf(:));

