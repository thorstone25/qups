function x = nan2zero(x), x(isnan(x)) = 0;
% NAN2ZERO - Convert NaN values to a 0
%
% y = NAN2ZERO(x) returns an output where all nan-values in x are replaced
% with zeros.
%
% See also ISNAN