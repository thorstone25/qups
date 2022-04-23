function nv = struct2nvpair(s)
    nv = reshape(cat(1, shiftdim(fieldnames(s),-1), shiftdim(struct2cell(s),-1)), [], 1);
end