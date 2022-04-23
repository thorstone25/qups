function varargout = dealfun(fun, varargin)
%DEAL Deal inputs to outputs after applying fun to each input.
%    [A,B,C,...] = DEAL(fun,X,Y,Z,...) simply matches up the input and
%       output lists, but applies fun to each argument. It is the same as 
%       A=fun(X), B=fun(Y), C=fun(Z), ...
%    [A,B,C,...] = DEAL(X) copies the single input to all
%       the requested outputs.  It is the same as A=fun(X), B=fun(X), 
%       C=fun(X), ...
%
%    See also DEAL.

if (nargin-1)==1,
  varargout = cellfun(fun, varargin(ones(1,nargout)), 'UniformOutput', false);
else
  if nargout ~= (nargin-1)
    error(message('MATLAB:dealfun:narginNargoutMismatch'))
  end
  varargout = cellfun(fun, varargin, 'UniformOutput', false);
end