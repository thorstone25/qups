function r = cohfac(b, dim)
% COHFAC - Compute the coherence factor
% 
% r = COHFAC(b, dim) computes the coherence factor of b in dimension
% dim.
%
% r = COHFAC(b) computes the coherence using the last non-singular
% dimension of b.
% 
% About: The coherence factor is a metric of the similarity between two
% signals. In medical ultrasound, it has been used to measure and improve
% image quality.
% 
% Example:
% % Define the setup - we'll use plane-waves
% us = UltrasoundSystem(); % get a default system
% us.sequence = SequenceRadial('type', 'PW', 'angles', -21 : 0.5 : 21);
% us.fs = single(us.fs); % accelerate the processing
% 
% % set the imaging resolution
% lambda = us.sequence.c0 / us.xdc.fc;
% [us.scan.dz, us.scan.dx] = deal(lambda / 4);
% 
% % Generate Scatterers
% N = 1e5;
% pos = 1e-3*([0 0 1]' ...
%   + [scan.xb(1); scan.yb(1); scan.zb(1)] ...
%   + rand([3 N]) .* range([scan.xb; scan.yb; scan.zb],2)...
% );
% scat = Scatterers('pos', pos, 'amp', rand([1,N]), 'c0', us.sequence.c0);
% 
% % Compute the image
% chd = greens(us, scat); % compute the response
% brx = DAS(us, chd, 'keep_rx', true); % beamform the data, keeping the receive dimension
% btx = DAS(us, chd, 'keep_tx', true); % beamform the data, keeping the transmit dimension
% 
% % Compute the Coherence factor across transmit and receive
% rrx = cohfac(brx);
% rtx = cohfac(btx);
% 
% % Display the images
% figure;
% ax = nexttile(); 
% imagesc(us.scan, real(rrx));
% colormap gray; colorbar;
% title('Receive coherence');
% 
% ax(2) = nexttile(); 
% imagesc(us.scan, real(rtx));
% colormap gray; colorbar;
% title('Transmit coherence');
%
% linkaxes(ax); % set both axes to scroll together
% 
% See also DMAS SLSC PCF

arguments
    b {mustBeNumeric}
    dim {mustBeInteger, mustBePositive} = find(size(b)~=1, 1, 'last');
end

% compute coherence
r = (abs(sum(b, dim)).^2) ./ sum(abs(b).^2, dim) / size(b, dim);