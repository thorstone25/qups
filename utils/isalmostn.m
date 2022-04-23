function [tf, tol] = isalmostn(a, b, tol)

if nargin < 3, tol = 1e2*eps(max(max(a(:)), max(b(:)))); end

if(~all(size(a) == size(b)) )
    tf = false;
    return;
end

ntf = isnan(a) & isnan(b);
vtf = (a - b) < tol;
ttf = ntf | vtf;
tf = all(ttf(:));

