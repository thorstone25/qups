function x = rsqrt(x)
% RSQRT - Reciprical square root
%
% y = rsqrt(x) returns the reciprocal square root of x by computing
% 1 ./ sqrt(x), while handling MATLAB casting issues for gpuArrays or
% integer types. 
% 
% gpuArray types are cast to complex if x contains negative values. 
% Integer types are upcast to double.
%
% Example:
% a = rsqrt(int16(64));
% 
% b = single(1 : 10) .^ 2;
% assert(isalmostn(rsqrt(b), 1 ./ (1 : 10)))
%  
% c = rsqrt(-4);
% 
% 
% See also: SQRT

arguments, x {mustBeNumeric}, end
if isinteger(x),    x = cast(x,'double'); end % cannot rsqrt integers
if any(x < 0,'all'),x = 1 ./ sqrt(complex(x)); % must be done manually for gpuArrays
else,               x = 1 ./ sqrt(        x );
end
end