%% Create a system definition
% ambient parameters
c0     = 1500;      % ambient / beamforming sound speed for this example
rho0   = 1000;      % ambient density
fc     = 7.5e6;     % L12-3v central frequency
lambda = c0 / fc;   % (spatial) wavelength

xdc = TransducerArray.L12_3v(); % Verasonics L12-3v linear array
seq = SequenceRadial('type', 'PW', 'angles', -10 : 10, 'c0', c0); % -10 -> 10 degree plane wave sequence
x = (-20e-3 : lambda / 4 : 20e-3); % simulation / imaging x-axis
z = (-01e-3 : lambda / 4 : 39e-3); % simulation / imaging z-axis
scan = ScanCartesian('x', x, 'z', z, 'y', 0); % simulation / imaging grid
us = UltrasoundSystem( ...
    'fs', single(25e6), ... request a (multiple of) 25MHz sampling frequency
    'xdc', xdc, 'seq', seq, 'scan', scan ... include properties
    ); % full system definition

% NOTE: the object 'scan' and 'us.scan' are the same object. Modifying
% 'scan' e.g. with `scan.dx = 2e-4` will also modify 'us.scan'.
% 
% If you need to treat 'scan' and 'us.scan' separately, make a copy e.g.: 
%
%     `scan = copy(scan)`
%


%% Define the scaterring medium
% define a grid of point scaterrers for resolution test
[~, xs, ys, zs] = ndgrid(0, 1e-3*(-10 : 5 : 10), 0, 1e-3*(10 : 5 : 30));
pos = [xs; ys; zs]; % point target positions
scat = Scatterers('pos', pos(:,:), 'amp', ones([1,numel(xs)]), 'c0', c0); % construct a Scaterrer object

% background
c   = c0   * ones(scan.size);

% add sound speed layers
c(07.5e-3 <= us.scan.z & us.scan.z <= 12.5e-3,:) = 1400;
c(12.5e-3 <= us.scan.z & us.scan.z <= 17.5e-3,:) = 1450;
c(17.5e-3 <= us.scan.z & us.scan.z <= 22.5e-3,:) = 1600;
c(22.5e-3 <= us.scan.z & us.scan.z <= 27.5e-3,:) = 1550;
c(27.5e-3 <= us.scan.z & us.scan.z <= 39.5e-3,:) = 1500;

% make the impedance constant
rho = rho0 * (c0 ./ c);

% add density scatterers
psn = scan.positions();
for j = 1:scat.numScat
    d = vecnorm(scat.pos(:,j) - psn,2,1); % distance to the scatterer
    i = argmin(d, [], 'all', 'linear'); % closest pixel
    rho(i) = rho(i) + rho0*scat.amp(j); % interpret amplitude as additive w.r.t. ambient
end

% create a sampled medium
med = Medium.Sampled(scan, c, rho, 'c0', c0, 'rho0', rho0);


%% Display the simulation
hf = figure('Name', 'Simulation Geometry');
imagesc(med, scan, 'props', 'c');
hold on; 
plot(us);
title(hf.Name);

hf = figure('Name', 'Excitation Pulse');
plot(us.seq.pulse, '.-');
title(hf.Name);

hf = figure('Name', 'Transducer Impulse Response');
plot(us.xdc.impulse, '.-');
title(hf.Name);


%% Run a FDTD simulation (requires k-Wave)
% optional: start a parallel pool to accelerate processing
parallelize = true;
if parallelize
    gpus = gpuDeviceCount(); % number of gpus available locally
    if gpus
        clu = parcluster('local'); % get default ProcessPool
        clu.NumWorkers = 2 * gpus; % 2 sims per gpu
        clu.NumThreads = 2; % 2 threads per gpu
        try parpool(clu, clu.NumWorkers); end
    else
        try parpool Threads; end % run a simulation per compute thread
    end
end

% launch the simulation, increasing the sampling frequency to satisfy the CFL
chd0 = kspaceFirstOrder(us, med, scan, "CFL_max", 0.25, "PML", [20 64]);


%% Show the channel data
figure('Name', 'Simulation Channel Data'); 
him = imagesc(hilbert(chd0)); % show magnitude of the data in dB
animate(him, hilbert(chd0).data, 'loop', false, 'fs', 5); % make a short animation across transmits


%% Beamform into a b-mode image
% construct an apodization matrix for receive sensitivity windowing
a = us.apAcceptanceAngle(25); % 25 degrees is close to an f# of 1

% preprocess
chd = chd0;
chd = downmix(chd, floor(chd.fs / 50e6)); % at least 25MHz sampling frequency
chd = hilbert(chd); % get analytic signal
if gpuDeviceCount(), chd = gpuArray(chd); end % move to gpu if available

% use a naive beamformer
b_naive = DAS(us, chd, c0, "apod", a);

% display the image
figure('Name', 'B-mode');
nexttile();
him    = imagesc(us.scan, b_naive); colorbar; 
caxis(max(caxis) + [-40 0]);
colormap gray;
title(him.Parent, 'B-mode')

%% Eikonal beamforming (requires mex files)
% extract eikonal delays on receive
[~, tau_rx, ~] = bfEikonal(us, chd, med, 'delay_only', true);

% extract standard delays on tx
[~, ~, tau_tx] = bfDAS(us, chd, 'delay_only', true); 

% beamform using an eikonal beamformer (requires mex file compilation)
b_eik = bfDASLUT(us, chd, tau_rx, tau_tx, "apod", a);

nexttile();
him(2) = imagesc(us.scan, b_eik); colorbar;
linkaxes([him.Parent])
linkprop([him.Parent], 'CLim');
caxis(max(caxis) + [-40 0]);

title(him(1).Parent, 'Naive B-mode')
title(him(2).Parent, 'Eikonal Rx. B-mode')
