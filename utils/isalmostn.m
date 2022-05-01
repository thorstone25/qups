function [tf, tol] = isalmostn(a, b, tol)
% ISALMOSTN Are two arrays almost equal (omitnan)
% 
% tf = ISALMOSTN(a, b) computes wheteher two arrays are almost equal
%
% tf = ISALMOSTN(a, b, tol) uses the tolerance tol for comparing the
% difference rather than a hueristic based on the maximum value
%
% [tf, tol] = ISALMOSTN(...) returns the tolerance used for the difference

if nargin < 3, tol = 1e2*eps(max(max(a(:)), max(b(:)))); end

if(~all(size(a) == size(b)) )
    tf = false;
    return;
end

ntf = isnan(a) & isnan(b);
vtf = (a - b) < tol;
ttf = ntf | vtf;
tf = all(ttf(:));

