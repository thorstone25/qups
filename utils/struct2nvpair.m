% STRUCT2NVPAIR - Convert a struct to a cell array of Name/Value pairs.
%
% nv = STRUCT2NVPAIR(s) converts the struct s to cell array of Name/Value
% pairs. This is convenient for using structs to pass input arguments to
% functions requiring Name/Value pairs.
%
% Example:
% s = struct('a', 1, 'b', 2, 'c', 3);
% nv = struct2nvpair(s); % make into Name/Value
% s2 = struct(nv{:}); % recreate the struct using Name/Value pairs
% assert(isequal(s, s2));
% 
% See also STRUCT STRUCT2CELL FIELDNAMES
function nv = struct2nvpair(s), arguments, s (1,1) struct, end
nv = cat(1, shiftdim(fieldnames(s),-1), shiftdim(struct2cell(s),-1));
end