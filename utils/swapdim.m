function [x, ord] = swapdim(x, i, o)
ord = 1:max([ndims(x), i, o]);
ord(i) = o;
ord(o) = i;
x = permute(x, ord);