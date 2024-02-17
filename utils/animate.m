function [mvf, mvh] = animate(x, h, kwargs)
% ANIMATE - Animate imagesc data
%
% ANIMATE(x, h) animates the multi-dimensional data x by looping through
% the upper dimensions and iteratively updating the image handle h. The
% data will be plotted until the figure is closed.
%
% If h is an array of image handles and x is a corresponding cell array of
% image data, each image h(i) will be updated by the data x{i}.
%
% If x or x{i} is complex, the magnitude in dB will be displayed.
%
% ANIMATE(..., 'loop', false) plays through the animation once, rather than
% looping until the figure is closed.
%
% ANIMATE(..., 'fn', true) appends the frame number for each frame to the
% title. If 'loop' is false, the frame number is removed before the
% function returns.
% 
% ANIMATE(..., 'fs', fs) updates the image at a rate fs.
% 
% mvf = ANIMATE(...) returns a cell matrix of movie frames mvf for each
% unique figure. This can be used to construct a movie or gif of each
% figure.
% 
% [mvf, mvh] = ANIMATE(...) returns a cell matrix of movie frames for each
% axes. This can be used to construct a movie or gif of each axes.
% 
% NOTE: MATLAB execution will continue indefinitely while the animation is playing.
% Close the figure or press 'ctrl + c' in the command window to stop the
% animation.
%
% Example:
% % Simulate some data
% us = UltrasoundSystem(); % get a default system
% us.fs = single(us.fs); % use single precision for speed
% us.seq = SequenceRadial('type', 'PW', 'angles', -21:0.5:21);
% scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
% chd = greens(us, scat); % simulate the ChannelData
% 
% % Configure the image of the Channel Data
% figure;
% nexttile();
% h = imagesc(chd); % initialize the image
% caxis(max(caxis) + [-60 0]); % 60 dB dynamic range
% colorbar;
% title('Channel Data per Transmit');
% 
% % Animate the data across transmits 
% animate(chd.data, h, 'loop', false); % show once
%
% % Beamform the data
% b = DAS(us, chd, 'keep_tx', true); % B-mode images per tx
% 
% % Initialize the B-mode image
% nexttile();
% h(2) = imagesc(us.scan, b); % show the center tx
% colormap(h(2).Parent, 'gray');
% caxis(max(caxis) + [-60 0]); % 60 dB dynamic range
% colorbar;
% title('B-mode per Transmit');
%  
% % Animate both images across transmits
% mvf = animate({chd.data, b}, h, 'loop', false); % show once
% 
% % Create a movie
% vobj = VideoWriter('tmp', 'Motion JPEG AVI');
% vobj.FrameRate = 20; % set frame rate to 20 Hz
% vobj.open();
% vobj.writeVideo([mvf{:,1}]);
% vobj.close();
% 
% See also IMAGESC FRAME2GIF
arguments
    x {mustBeA(x, ["cell","gpuArray","double","single","logical","int64","int32","int16","int8","uint64","uint32","uint16","uint8"])} % data
    h (1,:) matlab.graphics.primitive.Image = inferHandles(x)
    kwargs.fs (1,1) {mustBePositive} = 20; % refresh rate (hertz)
    kwargs.loop (1,1) logical = true; % loop until cancelled
    kwargs.fn (1,1) logical = false; % whether to add frame number
end

% place data in a cell if it isn't already
if isnumeric(x) || islogical(x), x = {x}; end

% argument type checks - x must contain data that can be plotted
% cellfun(@mustBeReal, x);

% Get sizing info
I = numel(h);
Mi = cellfun(@(x) prod(size(x,3:max(3,ndims(x)))), x); % upper dimension sizing
M = unique(Mi);

% validity checks
assert(isscalar(M), "The number of images must be the same for all images (" + (Mi + ",") + ").");
assert(numel(h) == numel(x), "The number of image handles (" ...
    + numel(h) + ...
    ") must match the number of images (" ...
    + numel(x) + ")." ...
    );
hf = unique([h.Parent]); % parent of the image is an axes
% we need to keept taking the parent until we reach the parent figure
isfig = @(hf) arrayfun(@(hf) isa(hf, 'matlab.ui.Figure'), hf);
i = ~isfig(hf);
while any(i) % any handles are not figure handles ...
    hf(i) = [hf(i).Parent]; % get parent of the axes: figure or layout
    hf = unique(hf); % keep unique figures only
    i = ~isfig(hf); % find which handles are still not figures
end
F = numel(hf); % number of figures
if nargout >= 2, mvh = cell([M,I]); else, mvh = cell.empty; end
if nargout >= 1, mvf = cell([M,F]); else, mvf = cell.empty; end
if kwargs.fn, ttls = arrayfun(@(h) string(h.Parent.Title.String), h); end

while(all(isvalid(h)))
    for m = 1:M
        if ~all(isvalid(h)), break; end
        for i = 1:I, if isreal(x{i}), h(i).CData(:) = x{i}(:,:,m); else, h(i).CData(:) = mod2db(x{i}(:,:,m)); end, end% update image
        if kwargs.fn, for i = 1:I, h(i).Parent.Title.String = ttls(i) + " ("+m+")"; end, end% add index to the title
        % if isa(h, 'matlab.graphics.chart.primitive.Surface'), h(i).ZData(:) = h(i).CData(:); end % TODO: integrate surfaces
        if m == 1, drawnow; getframe(); end % toss a frame to avoid bug where the first frame has a different size
        drawnow limitrate; 
        tic;
        if ~isempty(mvh), for i = 1:I, mvh{m,i} = getframe(h (i).Parent); end, end %  get the frames 
        if ~isempty(mvf), for f = 1:F, mvf{m,f} = getframe(hf(f)       ); end, end %  get the frames 
        pause(max(0, (1/kwargs.fs) - toc)); % pause for at least 0 sec, at most the framerate interval 
    end
    if ~kwargs.loop, break; end
end
if kwargs.fn, for i = 1:I, h(i).Parent.Title.String = ttls(i); end, end % reset titles

function him = inferHandles(x)
% get the current figure
hf = gcf();
if isempty(hf.Children) % new figure; no axes
    % create images for this data
    if isnumeric(x) || islogical(x), x = {x}; end % -> cell

    % use a tiledlayout by default
    htl = tiledlayout(hf, 'flow');

    % squeeze data into first 2 dims
    x = cellfun(@squeeze, x, 'UniformOutput', false);

    % take modulos of complex data
    val = cellfun(@isreal, x);
    x(~val) = cellfun(@mod2db, x(~val), 'UniformOutput', false);

    % make images
    him = cellfun(@(x) imagesc(nexttile(htl), x(:,:,1)), x);

    % done
    return;

elseif isa(hf.Children,  'matlab.graphics.layout.TiledChartLayout')
    % parse the tree structure to extract axes handles
    hax = hf.Children.Children; 
elseif any(arrayfun(@(h) isa(h, 'matlab.graphics.axis.Axes'), hf.Children)) % (sub)plot(s)
    hax = hf.Children;
else
    error("Unable to infer plot handle; please explictly pass the handle.")
end

% parse to grab image handles in same order as the data
hax = hax(arrayfun(@(hax)isa(hax, 'matlab.graphics.axis.Axes'), hax)); % axes only
him = {hax.Children}; % image and plot handles
him = cellfun(@(h) {h(arrayfun(@(h)isa(h, 'matlab.graphics.primitive.Image'),h))}, him);
him = flip([him{:}]); % image handle array, in order of creation

% TODO: check sizing of data versus image handles?
