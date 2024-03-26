function x = nan2zero(x)
% NAN2ZERO - Convert NaN values to a 0
%
% y = NAN2ZERO(x) returns an output where all nan-values in x are replaced
% with zeros.
%
% Example:
% nan2zero([1 2 nan 4 5])
% 
% See also ISNAN
x(isnan(x)) = 0;