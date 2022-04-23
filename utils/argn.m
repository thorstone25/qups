function out = argn(n, fun, varargin)
% draw n(2) args, output n(1)th arg from output

if isscalar(n)
    m = n;
elseif length(n) == 2
    [n,m] = deal(n(1),n(2));
end
out = cell([m,1]);
[out{:}] = fun(varargin{:});
out = out{n};
end