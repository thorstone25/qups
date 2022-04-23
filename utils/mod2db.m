function y = mod2db(x), y = mag2db(abs(x));
% MOD2DB - modulus (magnitude) in dB
%
% y = mod2db(x) returns the modulus of x in dB. Underneath it calls the
% mag2db function on the absolute value of x.
%
% See also MAG2DB, ABS
 