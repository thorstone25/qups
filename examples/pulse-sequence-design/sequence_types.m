%% Sequence types
%
% This example demonstrates the types of pre-defined pulse sequences
% available in QUPS. QUPS's optimized beamformer (DAS) requires one of 
% these formats. The general beamformers (bfDAS, bfEikonal) are not as 
% strict. The time axis of ChannelData returned by simulation routines will
% also conform to these formats.
% 
% Pulse sequence design affects the point spread function (PSF) of the
% image by controlling the synthetic transmit aperture as well as the shape
% of the pulse.
% 

% setup acceleration (optional) 
try parpool Threads; end %#ok<TRYNC> 

%% Create multiple imaging configurations to compare/contrast
% ambient parameters
c0 = 1500; % ambient sound speed

% make an array of different transducers
xdcs = { ...
    TransducerConvex.C5_2v(), ...
    TransducerArray.L11_2v(), ...
    TransducerArray.L12_3v(), ...
    TransducerArray.L12_5v(), ...
    };

% Transducer names (for plotting later)
nmxdc = ["C5-2v", "L11-2v", "L12-3v", "L12-5v"];

% Set sequence params
th = -10 : 0.5 : 10; % PW angles
M = 64; % active aperture element size
S = [1 1 2 4]; % aperture stride: increased for larger apertures
zf  = +50e-3; %  focused  pulse focal depth
zdv = -20e-3; % diverging pulse focal depth

% create pulse sequences for each transducer
for i = numel(xdcs):-1:1
    % select transducer
    xdc = xdcs{i};

    % get the walking aperture apodization
    aptx = activeLinearApertureApodization(xdc.numel, M, S(i)); 
    
    % get the focused and diverging wave focal positions
    pf  = activeApertureFocalPositions(xdc, aptx, zf ); % focused
    pdv = activeApertureFocalPositions(xdc, aptx, zdv); % diverging

    % construct a FSA pulse sequence
    seq0(i) = Sequence('type','FSA', 'numPulse',xdc.numel);

    % construct a plane wave pulse sequence
    seqp(i) = SequenceRadial('type','PW','angles', -10 : 0.5 : 10);

    % construct a focus pulse sequence
    seqf(i) = Sequence('type','VS', 'focus', pf);
    seqf(i).apodization_ = aptx; % set the apodization explicitly

    % construct a diverging pulse sequence
    seqv(i) = Sequence('type','VS', 'focus', pdv);
    seqv(i).apodization_ = aptx; % set the apodization explicitly
end

% set the sound speed for all sequences
[seq0.c0, seqp.c0, seqf.c0, seqv.c0] = deal(c0);
seqs = {seq0, seqp, seqf, seqv}; % all sequences

% create a single scan to be used across all sequences
scan = ScanCartesian('xb', [-5e-3 5e-3], 'zb', [25e-3 35e-3]);

% set plot units
scan.xlabel = scan.xlabel + " (m)";
scan.zlabel = scan.zlabel + " (m)";

% construct UltrasoundSystems: transducers are along rows, sequences along
% columns
for i = numel(xdcs):-1:1
    for j = numel(seqs):-1:1
        us(i,j) = UltrasoundSystem( ...
            'xdc', xdcs{i}, ...
            'seq', seqs{j}(i), ...
            'scan', scan, ...
            'fs', single(25e6) ...
            );
    end
end

% set the resolution for all scans
[scan.dx, scan.dz] = deal(min([us.lambda]) / 4);


%% Simulate data for each configuration
% define a single point scatterer
zs = 30e-3; % scatterer depth
scat = Scatterers('pos', [0;0;zs], 'c0',c0);

% simulate data for each system
chd = repmat(ChannelData(), size(us)); % pre-allocate
for i = 1:numel(us)
    tic, 
    chd(i) = greens(us(i), scat); 
    toc, 
end


%% Beamform for each configuration
% beamforming options
fnbr = 1; % set the acceptance angle with an f#

bi = cell(size(us)); % pre-allocate
for i = 1:numel(us)
    tic, 
    % create the receive apodization matrix
    a = us(i).apAcceptanceAngle(atand(0.5/fnbr));

    % beamform
    bi{i} = DAS(us(i), chd(i), 'apod', a);

    % normalize to the number of pulses (for comparable scaling)
     bi{i} = bi{i} ./ sqrt(us(i).seq.numPulse);
    toc,
end

%% Display the images, showing the effect of the Transducer and Sequence
% plot options
dbr = 60; % dbRange
scale_global = true; % whether to use the same intensity scale

% make plots
figure('Name', 'PSF');
htl = tiledlayout(size(us,2), size(us,1));
title(htl, 'PSF vs. Transducer | Sequence');
clear him;
for i = numel(us):-1:1
    % make plot
    hax(i) = nexttile(htl, i);
    him(i) = imagesc(us(i).scan, bi{i}, hax(i));

    % title
    ixdc = cellfun(@(xdc) us(i).xdc == xdc, xdcs); % get transducer index
    title(hax(i), nmxdc(ixdc) + " | " + us(i).seq.type);

    % scaling
    caxis(max(caxis) + [-dbr 0]);
end

% formatting
linkaxes(hax); % make scrolling identical on all plots
colormap gray; % use grayscale for the whole figure

% set axis limits around the point scatterers
xlim([-2.5e-3, 2.5e-3] + mean(scat.pos(1,:)));
ylim([-2.5e-3, 2.5e-3] + mean(scat.pos(3,:)));

% use global scaling
if scale_global
linkprop(hax, 'CLim');
bmax = gather(max(cellfun(@(b)max(b,[],'all'), bi(:))));
caxis(mod2db(bmax) + [-dbr 0]);
end

%% Helper functions (only available when run inside this script)
function aptx = activeLinearApertureApodization(N,M,S)
% ACTIVELINEARAPERTUREAPODIZATION
%
% aptx = ACTIVELINEARAPERTUREAPODIZATION(N, M, S) makes an apodization matrix for an active
% aperture of size M with stride S for an aperture with N elements in
% total. This can be used to construct focal sequences.
%

arguments
    N  % total elements
    M = floor(N/2); % active elements
    S = 1; % stride
end

% apodization matrix (elems x txs)
aptx = cell2mat(arrayfun(@(i) {circshift([true(M,1);false(N-M,1)],i,1)}, 0 : S : N - M));

end

function pf = activeApertureFocalPositions(xdc, aptx, zf)
% ACTIVEAPERTUREFOCALPOSITIONS
%
% pf = ACTIVEAPERTUREFOCALPOSITIONS(xdc, aptx, zf) returns focal positions
% where the focal points pf are projected to a depth of zf through the
% center of the active aperture defined by the non-zero entreis of the
% apodization matrix aptx. For a TransducerArray, the projection is axial.
% For a TransducerConvex, the projection is radial, xdc.center.
%
%
% for Transducer xdc, apodization matr
arguments
    xdc Transducer {mustBeA(xdc, ["TransducerArray", "TransducerConvex"])}
    aptx (:,:) % (elems x txs)
    zf (1,1) double
end

% size check
if size(aptx,1) ~= xdc.numel
    error( ...
        "QUPS:example_sequence_types:sizeMismatch", ...
        "Expected an apodization matrix for " + xdc.numel ...
        + " elements, but instead there are " + size(aptx,1) ...
        + " elements for " + size(aptx,2) + " transmits." ...
        );
end

if isa(xdc, 'TransducerArray')
    xn = sub(xdc.positions(),1,1); % lateral element position
    clear xf;
    for j = size(aptx,2):-1:1
        xf(j) = mean(xn(aptx(:,j))); % lateral center
    end
    pf  = [0;0; zf] + [1;0;0].*xf; % focal positions

elseif isa(xdc, 'TransducerConvex')
    az = xdc.orientations(); % element azimuth angles
    clear th;
    for j = size(aptx,2):-1:1
        th(j) = mean(az(  aptx(:,j))  ); % beam center angle
    end
    nf  = [sind(th); 0*th; cosd(th)]; % beam vector
    pf  = xdc.center + (xdc.radius + zf ) * nf; % focal positions
end
end

