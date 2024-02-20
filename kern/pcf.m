function [w, sf] = pcf(b, dim, gamma, kwargs)
% PCF - Phase Coherence Factor
% 
% w = PCF(b) computes the phase coherence factor w from the pre-sum
% image b.
%
% w = PCF(b, dim) operates across dimension dim. The default is the last
% non-singular dimension of b.
% 
% w = PCF(b, dim, gamma) adjusts the out-of-focus sensitivity. Larger
% values of gamma provide more suppresion. The default is 1.
% 
% w = PCF(..., 'unwrap', 'auxiliary') uses the auxiliary phase described in
% [1] to unwrap the phase prior to computing the phase diversity. This is
% the default.
% 
% [w, sf] = PCF(...) additionally returns the phase diversity estimate sf
% in radians.
% 
% About: The phase coherence factor is a metric of the variance in the 
% phase across the aperture. In medical ultrasound, it can been used to 
% weight the aperture to improve image quality.
% 
% References:
% [1] J. Camacho, M. Parrilla and C. Fritsch, "Phase Coherence Imaging," 
% in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
% vol. 56, no. 5, pp. 958-974, May 2009,
% doi: <a href="matlab:web('doi.org/10.1109/TUFFC.2009.1128')">10.1109/TUFFC.2009.1128</a>
% 
% Example:
% % Choose a transducer
% xdc = TransducerConvex.C5_2v();
% 
% % Create a sector scan sequence
% th  = asind(linspace(sind(-35), sind(35), 192)); % linear spacing in sin(th)
% rf  = xdc.radius + 50e-3; % 50mm focal depth (w.r.t. center of convex radius)
% atx = abs(xdc.orientations()' - th) <= 30; % active aperture of +/- 30 degrees
% seq = SequenceRadial('type', 'FC', 'ranges', rf, 'angles', th, 'apex', xdc.center, 'apd', atx);
% 
% % Create a sector scan imaging region
% scan = ScanPolar('rb', xdc.radius + [0, rf], 'a', th, 'origin', xdc.center); % along the scan lines
% us = UltrasoundSystem('seq', seq, 'xdc', xdc, 'scan', scan, 'fs', single(4*xdc.fc)); % get a default system
% us.scan.dr = us.lambda / 4; % set the imaging range resolution
%
% % Generate Scatterers
% pos  = 1e-3.*[sind(-10) 0 cosd(-10)]'*(5 : 5 : 30); % scatterer positions
% scat = Scatterers('pos', pos, 'c0', seq.c0);
% 
% % Compute the image
% chd = greens(us, scat); % compute the response
% b = DAS(us, chd, 'keep_rx', true); % beamform the data, keeping the receive aperture
% rxdim = ndims(b); % rx is the last dimension
% 
% % Compute the Phase Coherence Factor across the receive aperture
% w = pcf(b, rxdim);
%
% % Display the images
% bs = {sum(b,rxdim), w, sum(b.*w,rxdim)}; % images
% ttls = ["B-mode", "Phase Coherence Factor (PCF)", "PCF-weighted B-mode"]; % titles
% figure;
% colormap gray;
% for i = 1 : numel(bs)
%     nexttile();
%     imagesc(us.scan, bs{i});
%     title(ttls(i));
%     colorbar;
% end
% 
% See also SLSC DMAS COHFAC

arguments
    b {mustBeNumeric, mustBeComplex}
    dim (1,1) double {mustBeInteger, mustBePositive} = max([1, find(size(b)~=1, 1, 'last')]);
    gamma (1,1) double {mustBeReal} = 1
    kwargs.unwrap (1,1) string {mustBeMember(kwargs.unwrap, ["auxiliary"])} = "auxiliary";
end

switch kwargs.unwrap
    case "auxiliary"
        % compute the phase and it's standard deviation
        phi = angle(b); % phase (radians)
        s0  = std(phi,1,dim,"omitnan"); % ref std.

        % compute the auxiliary phase and it's standard deviation
        phi = phi - pi * sign(phi);
        sa  = std(phi,1,dim,"omitnan");

        % get scaling factor as the lesser standard deviation
        sf = min(s0, sa, "omitnan");
    otherwise
        % TODO: add smarter methods of accounting for wrap-around
        % 
        % Method 1) 
        % a) sort data then 
        % b) iteratively
        %     b-i) wrap min value 
        %     b-ii) get std 
        % c) then take min std over all from (b)

end

% standard deviation of the distribution U(-pi, pi)
sg0 = sqrt(pi/3); 

% phase coherence factor for scaling the image
w = max(0, 1 - (gamma / sg0) .* sf);

function mustBeComplex(b), if isreal(b), throwAsCaller(MException("QUPS:pcf:realInput","Input must be complex.")); end
