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
% seq = SequenceRadial('type', 'PW', 'angles', -21 : 0.5 : 21,'c0',1500);
% us = UltrasoundSystem('seq', seq); % get a default system
% us.fs = single(us.fs); % accelerate the processing
% 
% % Generate Scatterers
% spw = 10; % scatterers / wavelength
% P = 5; % wavelengths in region
% N = (spw * P) ^ 2; % total number of scatterers
% dp = P .* [1 0 1] * us.lambda; % P x P wavelengths (lateral x axial)
% p0 = [0 0 1e2*us.lambda] - dp/2; % offset to depth of 100 wavelengths
% pos = p0' + dp' .* rand([3 N]); % scatterers through the grid
% scat = Scatterers('pos', pos, 'amp', randn([1,N]), 'c0', us.seq.c0);
% 
% % set the imaging region via boundaries, then resolution
% us.scan = ScanCartesian('xb', p0(1) + [0 dp(1)], 'zb', p0(3) + [0 dp(3)]);
% [us.scan.dz, us.scan.dx] = deal(us.lambda / 8);
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