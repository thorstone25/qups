function reformat1liners(fld, kwargs)
% reformat1liners - Reformat all 1-line functions for all m-files
%
% reformat1liners(fld) reformats all m-files in the folder fld or any of
% it's subfolders so that one-liners with a following line of help text
% follow a stricter formatting standard by inserting newlines where
% appropriate. It is reformatted to a function declaration and signature,
% followed by a comment block of help text, followed by a matching end,
% each on a separate line.
%
% The default folder is the current directory.
arguments
    fld (1,1) string {mustBeFolder} = pwd
    kwargs.IncludeSubFolders = true
    kwargs.verbose = true
end

% get all files to be reformatted
if kwargs.IncludeSubFolders
    fls = dir(fullfile(fld, "**", "*.m"));
else
    fls = dir(fullfile(fld,       "*.m"));
end
fls = string(fullfile({fls.folder}, {fls.name}));
% fls = "src/ChannelData.m"; % DEBUG

%% build one-liner pattern
w1 = whitespacePattern(1,Inf); % alias
w0 = whitespacePattern(0,Inf); % alias
wc = wildcardPattern; % alias
vn = alphanumericsPattern | "_"; % variable name characters
fpat = ... function pattern
    "function" + w1 ... declaration
    + optionalPattern(("["+asManyOfPattern(vn+","+w0)+vn+"]" | vn)+w0+"="+w0) ... multiple outputs
    + asManyOfPattern(vn | ".", 1) ... valid func name
    + optionalPattern("("+asManyOfPattern(vn+","+w0)+vn+")") ... inputs
;
sep = ("," | ";"); % inline separator

decl = fpat+w0+sep;
oneL = decl + wc + "end";

%% reformat any 1-liners into multi-line with comments moved just below the
% function declaration
for fl = fls
    txt = readlines(fl); % read file line by line
    i = contains(txt,"function") & contains(txt,"end"); % mandatory
    i(i) = contains(txt(i), oneL); % see if this is truly a one-liner
    if kwargs.verbose, [~,nm] = fileparts(fl); disp(nm+": "+nnz(i)+" 1-liners."), end
    i = find(i); % we need line indices
    code = txt(i);
    if isempty(i), continue; end % short-circuit

    % extract pieces: declaration, body, "end", help-text
    dec = extract(code, decl);
    bod = extractBetween(code, decl, "end", "Boundaries", "exclusive");
    % hlp_after = extractAfter(code, oneL+wc+("%"|"...")); % if help is one-liner
    hlp_below = extractAfter(txt(i+1), lineBoundary+w0+"%");
    
    % prioritize help text
    hlp = hlp_below;
    % hlp = hlp_below;
    % j = ismissing(hlp);
    % hlp(j) = hlp_after(j);
    
    % re-write: insert newlines and move help text
    j = ~ismissing(hlp); % help available
    code    = replaceBetween(code   , decl, bod, newline);
    code(j) = replaceBetween(code(j), decl, newline, newline+"%"+hlp(j));
    code    = replaceBetween(code   , bod, "end",newline);
    
    % write to file
    txt(i) = code;
    writelines(txt, fl);
    
    % i = contains(txt, fpat + sep + wildcardPattern + (sep | w1) + "end" + wildcardPattern);
end

