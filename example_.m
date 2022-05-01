%% QUPS - Quick Ultrasound Processing & Simulation
%% What?
% QUPS is designed to be an accessible, verastile, shareable, lightweight codebase 
% designed to make developing and running ultrasound algorithms quick and easy! 
% It offers a standardization of data formatting that serves most common pulse-echo 
% ultrasound systems and supports simulation programs such as k-Wave or Fullwave.
%% Why?
% This package seeks to lower the barrier to entry for doing research in ultrasound 
% by offering a lightweight, easy to use package. Most of the underlying implementations 
% are here so that users can get to their work more efficiently! 
% 
% Currently you can:
%% 
% * Create linear and curvilinear (convex) transducers
% * Define full-synthetic-aperture, focused, diverging, or plane-wave transmits
% * Beamform to create a b-mode image
% * Simulate point targets (via MUST or FieldII) or distributed media (via k-Wave)
%% 
% Please submit issues or feature requests via <https://github.com/thorstone25/qups/issues 
% github>!
% 
% 
% 
% 
% Setup the workspace
% (Run once)

 
%#ok<*UNRCH> ignore unreachable code due to constant values
%#ok<*BDLGI> ignore casting numbers to logical values
device = -logical(gpuDeviceCount); % select 0 for cpu, -1 for gpu if you have one
% setup;  % add all the necessary paths
setup parallel cache; % start a parpool for faster processing | recompile and cache local mex/CUDA binaries

try hw = waitbar(0, 'test', 'Visible', 'off'); close(hw); end %#ok<TRYNC> % R2021b bug - errors on first call
%% 
%% Create a simple simulation
% Choose a Target

target_depth = 1e-3 * 30;
switch "single"
    case 'single', targ = Target('pos', [0;0;target_depth], 'c0', 1500); % simple point target
    case 'array' , targ = Target('pos', [0;0;1] * 1e-3*(10:5:50), 'c0', 1500); % targets every 5mm
    case 'diffuse', targ = Target('pos', 1e-3*[40;10;60].*([1;0;1].*rand(3,2501)-[0.5;0;0]), 'c0', 1500); % diffuse scattering
end
targ.rho_scat = 2; % make density scatterers at 2x the density
targ.scat_mode = 'ratio'; 
% Choose a transducer

switch "L11-5V"
    case 'L11-5V', xdc = TransducerArray.L11_5V(); % linear array
    case 'L12-3V', xdc = TransducerArray.L12_3V(); % another linear array
    case 'C5-2V' , xdc = TransducerConvex.C5_2V(); % convex array
end
% Choose the simulation region
% Make sure it covers the target and the transducer!

switch "large"
    case 'small', tscan = ScanCartesian('x',linspace(-20e-3, 20e-3, 16), 'z', linspace(  0e-3, 60e-3, 16));
    case 'large', tscan = ScanCartesian('x',linspace(-50e-3, 50e-3, 16), 'z', linspace(-30e-3, 60e-3, 16));
end
% Choose a transmit sequence
% For linear transducers only!

if isa(xdc, 'TransducerArray') 
    switch "Focused"
        case 'FSA', seq = Sequence('type', 'FSA', 'c0', targ.c0, 'numPulse', xdc.numel); % set a Full Synthetic Aperture (FSA) sequence
        case 'Plane-wave', 
            [amin, amax, Na] = deal( -25 ,  25 , 26 );
            seq = SequenceRadial('type', 'PW', ...
                'ranges', 1, 'angles',  linspace(amin, amax, Na), 'c0', targ.c0); % Plane Wave (PW) sequence
        case 'Focused', 
            [xmin, xmax, Nx] = deal( -10 ,  10 , 21 );
            seq = Sequence('type', 'VS', 'c0', targ.c0, ...
        'focus', [1;0;0] .* 1e-3*linspace(xmin, xmax, Nx) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 2mm
        ...'focus', [1;0;0] .* 1e-3*(-10 : 0.2 : 10) + [0;0;target_depth] ... % translating aperture: depth of 30mm, lateral stride of 0.2 mm
        ); 
    end

    
%% 
% For convex transducers only!

elseif isa(xdc, 'TransducerConvex') 
    switch "sector"
        case 'sector'
            [amin, amax, Na] = deal( -40 ,  40 , 41 );
            seq = SequenceRadial('type', 'VS', 'c0', targ.c0, ...
                'angles', linspace(amin, amax, Na), ...
                'ranges', norm(xdc.center) + 35e-3, 'apex', xdc.center ...
                ); % sector scan sequence
        case 'FSA'
            seq = SequenceRadial('type', 'FSA', 'numPulse', xdc.numel, 'c0', targ.c0 ...
                , 'apex', xdc.center ... % for convex transducer only!
                ); % FSA - convex sequence
    end
end
% Choose an imaging region

% For linear transducers only!
if isa(xdc, 'TransducerArray') 
    % set the scan at the edge of the transducer
    pn = xdc.positions(); % element positions
    xb = pn(1,[1,end]); % x-limits are the edge of the aperture
    zb = [-10e-3, 10e-3] + [min(targ.pos(3,:)), max(targ.pos(3,:))]; % z-limits surround the point target
    
    scan = ScanCartesian(...
        'x', linspace(xb(1), xb(end), 2^9), ...
        'z', linspace(zb(1), zb(end), 2^9) ...
        ); % X x Z scan

% For convex transducers only!
elseif isa(xdc, 'TransducerConvex') 

    % use with a SequenceRadial!
    scan = ScanPolar('origin', seq.apex, 'a', -40:0.5:40, ...
        'r', norm(seq.apex) + linspace(0, 2*max(vecnorm(targ.pos,2,1)), 2^9)...
        ); % R x A scan

end
%% 
% 

% Show the transducer's impulse response
% figure; plot(xdc.impulse, '.-'); title('Element Impulse Response'); 
%  Plot configuration of the simulation

figure; hold on; title('Geometry');
imagesc(targ, tscan); colorbar; % show the background medium for the simulation/imaging region
plot(xdc, 'r+', 'DisplayName', 'Elements'); % elements
hps = plot(scan, 'y.', 'MarkerSize', 0.5, 'DisplayName', 'Image'); % the imaging points
switch seq.type % show the transmit sequence
    case 'PW', plot(seq, 3e-2, 'k.', 'DisplayName', 'Tx Sequence'); % scale the vectors for the plot
    otherwise, plot(seq, 'k.', 'DisplayName', 'Tx Sequence'); % plot focal points, if they exist
end
plot(targ, 'k', 'LineStyle', 'none', 'Marker', 'diamond', 'MarkerSize', 5, 'DisplayName', 'Scatterers'); % point scatterers
hl = legend(gca, 'Location','bestoutside'); 
set(gca, 'YDir', 'reverse'); % set transducer at the top of the image
%% Simulate the Point Scatterer(s)

 
%% 
% 

% Construct an UltrasoundSystem object, combining all of these properties
us = UltrasoundSystem('xdc', xdc, 'sequence', seq, 'scan', scan, 'fs', 40e6);

% Simulate a point target
% run on CPU to use spline interpolation
switch "SIMUS"
    case 'FieldII', chd0 = calc_scat_all(us, targ, [1,1], 'device', device, 'interp', 'cubic'); % use FieldII, 
    case 'SIMUS'  , chd0 = simus(us, targ, 'periods', 1, 'dims', 3, 'interp', 'freq'); % use MUST: note that we have to use a tone burst or LFM chirp, not seq.pulse
    case 'Greens' , chd0 = comp_RS_FSA(us, targ, [1,1], 'method', 'interpd', 'device', device, 'interp', 'cubic'); % use a Greens function with a GPU if available!
end
chd0

% display the channel data across the transmits
chd = mod2db(chd0); % == 20*log10(abs(x)) -> the power of the signl in dB
figure; h = imagesc(chd, 1); colormap jet; colorbar; caxis(gather([-80 0] + (max(chd.data(chd.data < inf)))))
xlabel('Channel'); ylabel('Time (s)'); ylim([min(chd.time(:)), max(chd.time(:))]);
for m = 1:size(chd.data,3), if isvalid(h), h.CData(:) = chd.data(:,:,m); h.YData(:) = chd.time(:,:,min(m,size(chd.time,3))); drawnow limitrate; title(h.Parent, "Tx " + m); pause(1/10); end, end
%% Create a B-mode Image

 
%% 
% 

% Precondition the data
chd = chd0;
chd = single(chd); % use less data
chd.data = chd.data - mean(chd.data, 1, 'omitnan'); % remove DC 
if device, chd = gpuArray(chd); end % move data to GPU
if isreal(chd.data), chd = hilbert(chd, 2^nextpow2(chd.T)); end % apply hilbert on real data
chd = rectifyt0(chd, 'linear'); % linearly resample if needed to align t0 time axis

% Run a simple DAS algorithm
switch class(xdc)
    case 'TransducerArray' , scale = 1e-3;
    case 'TransducerConvex', scale = 1;
end
switch seq.type
    case "VS", 

%% 
% Choose an apodization method for Virtual Source (VS) transmits

        switch "accept"
            case 'multiline', apod = multilineApodization(us.scan, us.sequence);
            case 'scanline', apod = scanlineApodization(us.scan, us.sequence);
            case 'translating', apod = translatingApertureApodization(us.scan, us.sequence, us.rx, 19.8*scale);
            case 'aperture-growth', apod = apertureGrowthApodization(us.scan, us.sequence, us.rx, 1.8);
            case 'accept', apod = acceptanceAngleApodization(us.scan, us.sequence, us.rx, 55); 
            case 'none', apod = 1;
        end        
    otherwise, apod = 1; % apodization profiles for plane-wave/FSA not implemented :(
end
b = DAS(us, chd, struct('c0', targ.c0), [], 'device', device, 'interp', 'cubic', 'apod', apod);

% show the image
b_im = mod2db(b); % convert to power in dB
figure; imagesc(scan, b_im, [-80, 0] + max(b_im(:))); % display with 60dB dynamic range
colormap gray; colorbar;
%% Scan Convert for a Sector Scan

if ~isa(scan, 'ScanCartesian')
    scanc = scanCartesian(scan); % mind the caps! this gives us a ScanCartesian
    [scanc.nx, scanc.nz] = deal(2^9); % don't change boundaries, but increase resolution
    bc_im = (scanConvert(scan, mod2db(b), scanc)); % interpolate in dB - could also do abs, as long as it's not complex!

    imagesc(scanc, bc_im, [-80, 0] + max(bc_im(bc_im < Inf)));
    colormap gray; colorbar;
end
%% Run a k-Wave Sim (requires k-Wave)

 
% Note: For a FSA acquisition, that means 1 sim per element!

% use some low-res simulation parameters 
% Takes about 1 min per sim, 1 sim per tx
bds = cat(1,tscan.xb, tscan.yb, tscan.zb);
kwv_args = {'dims', 2, 'CFL_max', 0.5, 'resolution_ratio', 0.25, 'PML_min', 2^4, 'PML_max', 2^6, ... 
    'bounds', bds, 'PlotSim', false ...
    };
if device, chd1 = kspaceFirstOrderND(us, targ, [1,1], kwv_args{:}); % use just a single sub-element for speed
else,      chd1 = kspaceFirstOrderND(us, targ, [1,1], kwv_args{:}, 'DataCast', 'single'); end
% Show the data

chd = copy(chd1);
% chd.data(chd.time < 20e-6,:,:) = nan; % normalize display to data after 20us 
figure; h = imagesc(chd, 1); colormap jet; colorbar; caxis([-80 0] + mod2db(max(chd.data(:), [], 'omitnan')))
xlabel('Channel'); ylabel('Time (s)');
for m = 1:size(chd.data,3), if isvalid(h), h.CData(:) = mod2db(chd.data(:,:,m)); drawnow limitrate; title(h.Parent, "Tx " + m); pause(1/10); end, end
%% Create an image 

 
% Precondition the data
chd = copy(chd1);
chd = single(chd); % use less data
chd.data(chd.time < 10e-6 | chd.time > 50e-6,:,:) = 0; % clear direct feedback
chd.data = chd.data - mean(chd.data, 1, 'omitnan'); % remove DC bias
if device, chd = gpuArray(chd); end % move data to GPU
if isreal(chd.data), chd = hilbert(chd, 2^nextpow2(chd.T)); end % apply hilbert on real data
chd = rectifyt0(chd, 'linear'); % linearly resample if needed to align t0 time axis

% Run a simple DAS algorithm
b = DAS(us, chd, struct('c0', targ.c0), [], 'device', -1, 'interp', 'linear');

% show the image
b_im = mod2db(b); % convert to power in dB
figure; imagesc(scan, b_im, [-60, 0] + max(b_im(:))); % display with 60dB dynamic range
colormap gray; colorbar;
%% Scan Convert for a Sector Scan

if ~isa(scan, 'ScanCartesian')
    scanc = scanCartesian(scan); % mind the caps! this gives us a ScanCartesian
    [scanc.nx, scanc.nz] = deal(2^9); % don't change boundaries, but increase resolution
    bc_im = (scanConvert(scan, mod2db(b), scanc)); % interpolate in dB - could also do abs, as long as it's not complex!

    imagesc(scanc, bc_im, [-60, 0] + max(bc_im(bc_im < Inf)));
    colormap gray; colorbar;
end