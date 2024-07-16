function dbr(mode, rang)
% DBR - dB range for plotting
%
% DBR MODE set the colormap, colorbar, and range of the plot according
% to a series of preset MODEs.
%
% DBR MODE RANG sets the range of the colorbar.
%
% For the 'b-mode' and 'echo' presets, the range is in dB.
% For the 'phase' preset, the range is in degrees.
% 
% Example:
% figure;
% imagesc;
% dbr phase 35
% 
% pause(1);
% dbr b-mode;
% 
% pause(1); 
% dbr echo 30
% 
% pause(1);
% dbr corr 0.5
% 
% See also: SCAN/IMAGESC CHANNELDATA/IMAGESC
arguments
    mode (1,1) string {mustBeMember(mode, ["b-mode", "phase", "echo", "corr"])} = "b-mode"
    rang (1,1) string = defaultRange(mode)
end
rang = str2double(rang); % interpret as a number

colorbar; % always activate colorbar
ax = gca; h = ax.Children( arrayfun(@(h) ...
      isa(h,'matlab.graphics.primitive.Image'  ) ...
    | isa(h,'matlab.graphics.primitive.Surface'), ...
    ax.Children));
assert(~isempty(h), "Cannot identify an image for this axis."); % HACK - get the maximum value
% caxis auto; cmax = max(caxis); 
if exist('clim', 'file'), clm = @clim; else, clm = @caxis; end %#ok<CAXIS> - backwards compatibility
cmax = max(cellfun(@(x) max(x,[],'all','omitnan'),{h.CData}),[],'all','omitnan');
switch mode
    case "b-mode", colormap(ax, 'gray'); clm(cmax + [-rang   0 ]);
    case "echo"  , colormap(ax, 'jet' ); clm(cmax + [-rang   0 ]);
    case "phase" , colormap(ax, 'hsv' ); clm(0    + [-rang rang]);
    case "corr"  , colormap(ax, 'hot' ); clm(     + [rang    1 ]);
end

function r = defaultRange(mode)
switch mode
    case "b-mode", r = 40; % dB
    case "echo",   r = 60; % dB
    case "phase",  r = 180; % degrees
    case "corr",   r = 0; % proportion
    % otherwise, error("Not implemented");
end