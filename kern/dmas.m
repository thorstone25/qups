function b = dmas(bn, dim, L)
% DMAS - Delay-Multiply-And-Sum (DMAS)
% 
% b = DMAS(bn) computes the delay-multiply-and-sum image b from the b-mode 
% image per receive element bn.
%
% b = DMAS(bn, dim) specifies the receive dimension of bn. The default is 
% the last non-singular dimension of bn.
%
% b = DMAS(bn, dim, lags) specifies the element lags to be included. The
% default is (1 : size(bn,dim) - 1).
% 
% b = DMAS(bn, dim, L) where L is a scalar selects lags 1:L. To select only
% lag L, specify an array which includes 0 as in lags = [0, L]. 
%
% About: Delay-multiply-and-sum is a contrast enhancement method based on
% correlations between pairs of signals across the aperture.
%
% For complex numbers, this implementation preserves the phase in the
% multiplication stage, which is the complex analog of preserving the sign.
%
% References:
% [1] G. Matrone, A. S. Savoia, G. Caliano and G. Magenes, 
% "The Delay Multiply and Sum Beamforming Algorithm in Ultrasound B-Mode Medical Imaging," 
% in IEEE Transactions on Medical Imaging, vol. 34, no. 4, pp. 940-949, April 2015, 
% doi: <a href="matlab:web('https://doi.org/10.1109/TMI.2014.2371235')"> 10.1109/TMI.2014.2371235</a>.
% 
% Example:
% % Define the system
% us = UltrasoundSystem(); % get a default system
% 
% % Generate Scatterers
% scat = Scatterers('pos', 1e-3*[0,0,30]', 'c0', us.seq.c0);
% 
% % Compute the image
% chd = greens(us, scat); % compute the response
% fmod = us.xdc.fc; % baseband demodulation frequency
% chd = downmix(chd, fmod); % baseband the data
% brx = DAS(us, chd, 'keep_rx', true, 'fmod', fmod); % beamform the data at baseband, keeping the receive dimension
% ndim = ndims(brx); % get receive dimension
% bdas = sum(brx, ndim); % DAS
% bdmas = dmas(brx, ndim); % DMAS
% 
% % Display the images
% figure;
% ax = nexttile(); 
% imagesc(us.scan, bdas);
% colormap gray; colorbar;
% title('Delay-and-Sum');
% cax1 = caxis;
% 
% ax(2) = nexttile(); 
% imagesc(us.scan, bdmas);
% colormap gray; colorbar;
% title('Delay-Multiply-and-Sum');
% cax2 = caxis;
%
% linkaxes(ax); % set both axes to scroll together
% linkprop(ax, 'CLim'); % set magnitude scaling
% xlim([-02.5 02.5]*1e-3);
% ylim([ 27.5 32.5]*1e-3);
% caxis(max(cax1(2), cax2(2)) + [-60 0]); % 60dB dynamic range
% 
% See also SLSC COHFAC PCF

arguments
    bn {mustBeNumeric}
    dim {mustBeInteger, mustBePositive} = find(size(bn) ~= 1, 1, 'last');
    L  (1,:) {mustBeInteger, mustBeNonnegative} = 1:size(bn, dim)-1 % lags
end

b = 0; % init
N = size(bn, dim); % aperture length
if isscalar(L), lags = 1:L; else; lags = intersect(1:N-1, L); end
for i = lags % multiply and sum along all non-identical pairs
    b = b + sum(sub(bn,1:N-i,dim) .* sub(bn,1+i:N,dim), dim);
end

% re-scale the amplitude, preserve the phase (sign)
b = exp(1j*angle(b)) .* sqrt(abs(b));
