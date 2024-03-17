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
% See also: SCAN/IMAGESC CHANNELDATA/IMAGESC
arguments
    mode (1,1) string {mustBeMember(mode, ["b-mode", "phase", "echo"])} = "b-mode"
    rang (1,1) string = defaultRange(mode)
end
rang = str2double(rang); % interpret as a number

colorbar; % always activate colorbar
ax = gca; h = ax.Children(isa(ax.Children,'matlab.graphics.primitive.Image'));
assert(isscalar(h), "Cannot identify the image for this axis."); % HACK - get the maximum value
% caxis auto; cmax = max(caxis); 
cmax = max(h.CData,[],'all', 'omitnan');
switch mode
    case "b-mode", colormap(gca, 'gray'); caxis(cmax + [-rang   0 ]);
    case "echo"  , colormap(gca, 'jet' ); caxis(cmax + [-rang   0 ]);
    case "phase" , colormap(gca, 'hsv' ); caxis(0    + [-rang rang]);
end

function r = defaultRange(mode)
switch mode
    case "b-mode", r = 40; % dB
    case "echo",   r = 60; % dB
    case "phase", r = 180; % degrees
    otherwise, error("Not implemented");
end