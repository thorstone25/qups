function varargout = dealfun(fun, varargin)
% DEALFUN - Apply fun to inputs and deal return values to outputs
%
% [A,B,C,...] = DEALFUN(fun,X,Y,Z,...) simply matches up the input and
% output lists, but applies fun to each argument. It is the same as 
%       A=fun(X), B=fun(Y), C=fun(Z), ...
% [A,B,C,...] = DEALFUN(fun,X) copies the single input to all
% the requested outputs.  It is the same as 
%       A=fun(X), B=fun(X), C=fun(X), ...
%
% See also DEAL.


if (nargin-1)==1,
  varargout = cellfun(fun, varargin(ones(1,nargout)), 'UniformOutput', false);
else
  if ~nargout || (nargout ~= (nargin-1))
    error('QUPS:dealfun:narginNargoutMismatch',  ['Error using dealfun',newline,'The number of outputs should match the number of inputs.'])
  end
  varargout = cellfun(fun, varargin, 'UniformOutput', false);
end