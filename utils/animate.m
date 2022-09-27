function animate(h, x, kwargs)
% ANIMATE - Animate imagesc data
%
% ANIMATE(h, x) animates the multi-dimensional data x by looping through
% the upper dimensions and iteratively updating the image handle h. The
% data will be plotted until the figure is closed. 
%
% ANIMATE(..., 'loop', false) plays through the animation once, rather than
% looping until the figure is closed.
%
% ANIMATE(..., 'fs', fs) updates the image at a rate fs.
% 
% Note: MATLAB execution will be paused while the animation is playing.
% Close the figure to stop the animation.
%
% Example:
% % Simulate some data
% us = UltrasoundSystem(); % get a default system
% us.fs = single(us.fs); % use single precision for speed
% us.sequence = SequenceRadial('type', 'PW', 'angles', -21:0.5:21);
% scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
% chd = greens(us, scat); % simulate the ChannelData
% 
% % Configure the image of the Channel Data
% figure;
% chd_im = mod2db(chd);
% h = imagesc(chd_im); % initialize the image
% caxis(max(caxis) + [-60 0]); % 60 dB dynamic range
% colorbar;
% title('Channel Data per Transmit');
% 
% % Animate the data across transmits 
% animate(h, chd_im.data, 'loop', false); % show once
%
% % Beamform the data
% b = DAS(us, chd, 'keep_tx', true); % B-mode image
% bim = mod2db(b); % log-compression / envelope detection
% 
% % Initialize the B-mode image
% figure;
% h = imagesc(us.scan, bim(:,:,1)); % show the first image
% colormap gray;
% caxis(max(caxis) + [-60 0]); % 60 dB dynamic range
% colorbar;
% title('B-mode per Transmit');
%  
% % Animate the data across transmits
% animate(h, bim, 'loop', false); % show once
% 
% See also IMAGESC
arguments
    h (1,1) matlab.graphics.primitive.Image
    x {mustBeNumericOrLogical} % data
    kwargs.fs = 20; % hertz
    kwargs.loop = true; % loop until cancelled
end

M = prod(size(x,3:max(3,ndims(x))));
while(isvalid(h))
for m = 1:M
    if ~isvalid(h), break; end
    h.CData(:) = x(:,:,m);
    drawnow limitrate;
    pause(1/kwargs.fs);
end
if ~kwargs.loop, break; end
end

end