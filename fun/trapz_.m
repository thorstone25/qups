function z = trapz_(y, dim)
% TRAPZ_ - Alternate trapezoidal numerical integration
%
% z = TRAPZ_(y, dim) implements z = TRAPZ(y, dim) with fewer type-checking
% guardrails by using only sum, sub, and size.
%
% See also TRAPZ SUB SUM SIZE

    % x = sum(x, dim) - sum(sub(x, [1, size(x,dim)], dim), dim) / 2;
    z = (sum(sub(y, 1:size(y,dim)-1, dim), dim) + sum(sub(y, 2:size(y,dim), dim), dim)) / 2;
    % x = sum(movsum(x / 2, [0 1], dim, 'Endpoints', 'discard'), dim);
end