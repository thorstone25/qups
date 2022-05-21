function out = argn(n, fun, varargin)
% ARGN - Return the nth output from a function
%
% out = ARGN(n, fun, x) returns the nth output of calling the function fun 
% on x. It is the same as
%
% [~,...,~,out] = fun(x)
%
% where out is the nth output.
%
% out = ARGN([n,m], fun, x) draws m outputs from fun and returns the nth
% output.
%
% out = ARGN(..., fun, x, y, z, ...) applies fun to all following inputs.
%
% Example:
%
% % Make a 4D array
% A = rand([2,3,4,5]);
%
% % ask for output 3/3 
% % this vectorizes dimension 4 into dimension 3
% [~,~,o1] = size(A); 
% o2 = argn(3, @size, A); % use argn to ask for output 3/3
% assert(all([o1, o2] == 4*5))
%
% i1 = size(A, 3); % use the dim argument to get dim 3
% [~,~,i2,~] = size(A); % ask for output 3/4 - this does NOT vectorize dimension 4
% i3 = argn([3 4], @size, A); % use argn to ask for output 3/4
% assert(all([i1 i2 i3] == 4))
%
% See also DEAL DEALFUN

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