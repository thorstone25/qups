function x = trapz_(x, dim)
    % x = sum(x, dim) - sum(sub(x, [1, size(x,dim)], dim), dim) / 2;
    x = (sum(sub(x, 1:size(x,dim)-1, dim), dim) + sum(sub(x, 2:size(x,dim), dim), dim)) / 2;
    % x = sum(movsum(x / 2, [0 1], dim, 'Endpoints', 'discard'), dim);
end