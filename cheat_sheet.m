% CHEAT_SHEET - "Copy & paste" QUPS syntax and usage examples
%
% This script provides basic syntax and usage examples for a variety of
% QUPS classes, methods, and compute kernels (functions). It is meant to
% provide a list of "copy & paste" ready lines of code.

%% Transducers
c0 = 1500; 
fc = 6e6;

% define a pre-defined linear or phased array transducer
xdc = TransducerArray.L12_3v();

% define a custom linear array or phased array transducer
xdc = TransducerArray('fc', fc, 'pitch', c0/fc, 'numel', 128);

% define a pre-defined curvilinear array 
xdc = TransducerConvex.C5_2v();

% define a custom curvilinear array
xdc = TransducerConvex('fc', fc, 'pitch', c0/fc, 'numel', 128);

% define a pre-defined matrix transducer
xdc = TransducerMatrix.PO192O();

% define a custom matrix transducer
xdc = TransducerMatrix('fc', fc, 'pitch', [c0/fc c0/fc], 'numd', [16 16]);

% change number of elements
xdc.numel = 256;

% change element spacing
xdc.pitch = c0 / fc / 4; % lambda / 4

% (re)set the impulse response waveform (affects simulation only)
xdc.impulse = xdc.xdcImpulse(); % ideal impulse response
% or
xdc.impulse = Waveform.Delta(); % no impulse response

% Create a generic transducer with slightly random locations
pos = 1e-3 * [
    sort(rand([1 128]) - 0.5 + linspace(-20, 20, 128));  % x
    zeros([2 128]) % y/z
    ]; % random positions
xdc = TransducerGeneric('fc', fc, 'pos', pos, 'az', 0, 'el', 0);

%% Pulse Sequences
% Pulse sequences define groups of signals transmitted in quick 
% succession that are used to form an image. 
% 
% Typically, an excitation signal is sent to all transducer elements 
% with relative delays and amplitudes across the aperture to form a single 
% larger wavefront in the medium. This wavefront then propagates, creating 
% reflections at a known time and location. Beamforming uses wavefront 
% models to combine images from each transmit into a single sharper image.
% Pulse sequence design is a critical part of this beamforming process, and
% has a large effect on image quality.

% ------ FSA sequences ---------- %
% setup a full-synthetic-aperture (FSA) sequence
seq = Sequence('type', 'FSA', 'numPulse', xdc.numel);

% setup a Hadamard encoded transmit
% NOTE: array size N is supported only where N, N/12 or N/20 is a power of 2.
seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'apd', hadamard(xdc.numel));

% ---------- Plane Wave Sequences -------- %
% setup a plane wave (PW) sequence
seq = SequenceRadial('type','PW', 'angles', -20 : 0.5 : 20);

% ---------- Focused Sequences ---------- %
% setup a focused pulse (FC) sequence
zf = 60e-3; % focal depth 
xf = (-40 : 1 : 40) * 1e-3; % focal point lateral positions
pf = [0,0,1]'.*zf + [1,0,0]'.*xf; % focal points 
seq = Sequence('type','FC', 'focus', pf);

% setup a walking transmit aperture focused pulse (FC) for a linear transducer
xdc = TransducerArray();
pn = xdc.positions(); % element positions
Na = floor(xdc.numel/2); % active aperture size
Nv = xdc.numel - Na + 1; % number of transmits
apd = [ones([Na, Nv]); zeros([xdc.numel - Na, Nv])]; % transmit apodization (elems x txs)
for i = 1 : Nv
    % shift the active aperture over by i elements for the ith transmit
    apd(:,i) = circshift(apd(:,i), i-1, 1);

    % get focal positions centered on the active aperture
    pf(:,i) = mean(pn(:, logical(apd(:,i))),2); 
end
seq = Sequence('type','FC','focus',pf, 'apd', apd);

% setup a walking transmit aperture focused pulse (FC) for a curvilinear array
xdc = TransducerConvex();
th = xdc.orientations(); % element azimuth angles (deg)
Na = floor(xdc.numel/2); % active aperture size
Nv = xdc.numel - Na + 1; % number of transmits
apd = [ones([Na, Nv]); zeros([xdc.numel - Na, Nv])]; % transmit apodization (elems x txs)
for i = 1 : xdc.numel - Na + 1
    % shift the active aperture over by i elements for the ith transmit
    apd(:,i) = circshift(apd(:,i), i-1, 1);

    % get focal positions centered on the active aperture
    tha(:,i) = mean(th(logical(apd(:,i))),2); % get focal positions (centered)
end
rfocal = 60e-3; %% focal range
seq = SequenceRadial( ...
    'type','FC', ...
    'angles',tha, ...
    'ranges',norm(xdc.center) + rfocal, ...
    'apex',xdc.center ...
    ,'apd', apd ...
    );

% --------- Arbitrary Delay Sequences ----------- %
% Create an arbitrary-delay sequence
tau = randn(xdc.numel); % (elems x txs)
apd = hadamard(xdc.numel); % (elems x txs)
seq = SequenceGeneric('del', tau, 'apd', apd);

% ---------- Modify Sequence Parameters ---------- %
% set the beamforming sound speed (m/s)
seq.c0 = 1500; 

% set the transmit pulse waveform (affects simulation only)
seq.pulse = Waveform.Delta();

%% Imaging Regions (/ Simulation Regions)
% setup a 2D Cartesian imaging region
x = -20e-3 : 0.2e-3 : 20e-3;
z = 0      : 0.2e-3 : 40e-3;
scan = ScanCartesian('x',x,'z',z);

% setup a 3D Cartesian imaging region
x = -20e-3 : 0.2e-3 : 20e-3;
y = -20e-3 : 0.2e-3 : 20e-3;
z = 0      : 0.2e-3 : 40e-3;
scan = ScanCartesian('x',x,'y',y,'z',z);

% set the resolution for a Cartesian region
scan.dx = 0.1e-3;
scan.dz = 0.1e-3;

% setup a 2D polar imaging region
r =   0 : 0.2e-3 : 40e-3;
a = -40 : 0.5    : 40   ;
scan = ScanPolar('r',r,'a',a);

% setup a 3D cylindrical imaging region
r =   0 : 0.2e-3 : 40e-3;
a = -40 : 0.5    : 40   ;
y = -20e-3 : 0.2e-3 : 20e-3;
scan = ScanPolar('r',r,'a',a,'y',y);

% set the resolution for a Polar region
scan.dr = 0.1e-3;
scan.da = 0.25;

% create a Generic Scan with an arbitrary axes-to-pixel transform function
scan = ScanGeneric('u', 1:100, 'v', -90:90, 'w', 1 : 100);
scan.trans = @(u,v,w) cat(1,u,cosd(v), exp(w));

%% Beamforming
% create some example objects and data for this section
us = UltrasoundSystem(); % default system
scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % point scatterer(s)
med = Medium('c0', scat.c0); % FDTD medium
chd = greens(us, scat); % get some data for this example 

% ------------- B-mode Images ----------------- %
% standard delay-and-sum (DAS)
b = bfDAS(us, chd);

% compute optimized standard delay-and-sum (DAS)
% (incompatible with non-standard Sequences)
b = DAS(us, chd);

% look-up table (LUT) delay-and-sum (DAS)
tau_tx = zeros([us.scan.size, chd.M]); % delay tensor: pixels x transmits
tau_rx = zeros([us.scan.size, chd.N]); % delay tensor: pixels x receives
b = bfDASLUT(us, chd, tau_tx, tau_rx);

% eikonal equation beamformer (FSA only)
b = bfEikonal(us, chd, med);

% eikonal equation beamformer with a high sound speed resolution grid
% currently lambda/10 or finer recommended (FSA only)
cscan = copy(us.scan); 
cscan.dx = us.lambda / 10; 
cscan.dz = us.lambda / 10;
b = bfEikonal(us, chd, med, cscan);

% frequency-domain adjoint Green's function beamformer
% Note: you may get poor performance on virtual source ('FC'/'DV'/'VS')
% sequences, convex arrays, or rotated/offset transducers
b = bfAdjoint(us, chd);

% Stolt's f-k-migration with FFT padding and output scan (PW only)
uspw = copy(us); % same system
uspw.seq = SequenceRadial('type', 'PW','angles',-10:1/4:10); % use plane waves instead
chdpw = focusTx(uspw, chd); % focus FSA into PW pulses
[b, bscan] = bfMigration(uspw, chdpw, [2*chd.T, 4*chd.N]);
% image using this scan i.e. with `imagesc(bscan, b);`

% ----------------- Aperture Reduction Functions --------------- %
bn = DAS(us, chd, 'keep_rx', true); % first, preserve the receive dimension
ndim = ndims(bn); % the last dimension is the receive dimension, after combining transmits

% coherence factor
cf = cohfac(bn, ndim);

% short-lag spatial-coherence (SLSC) or short-lag angular-coherence (SLAC)
s = slsc(bn, ndim);

% delay-multiply-and-sum (DMAS)
bd = dmas(bn, ndim);

% ---------------- Apodization windows -------------------- %
% 45 deg acceptance angle weights on receive
a = us.apAcceptanceAngle(45); 

% set aperture growth to limit f# >= 1
a = us.apApertureGrowth(1);

% apply apodization when beamforming 
% (works in most cases with most beamformers)
b = DAS(us, chd, 'apod', a);

%% Channel Data
% This section demonstrates operations that can be efficiently done on
% ChannelData objects. All operations that affect the time axes will modify
% the 't0' and 'fs' properties appropriately. 
% 
% The size of the 't0' property should be scalar in the time and receive 
% dimensions i.e. all(size(chd.t0, [chd.tdim chd.ndim]) == 1) should be
% true. QUPS will not stop you from breaking this rule, but unexpected
% errors later on may occur.


% create some example objects and data for this section
us = UltrasoundSystem(); % default system
seq = SequenceRadial('type','PW','angles',-10:10); % alternate pulse sequence
scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % point scatterer(s)
chd_fsa = greens(us, scat); % get some FSA data for this example 
chd = chd_fsa; % alias

% -------------- focusing (/ beamforming) ---------------- %

% focus full synthetic aperture (FSA) data into another pulse sequence
us_foc = copy(us); % make a new system
us_foc.seq = seq; % set the desired pulse sequence
chd_foc = focusTx(us_foc, chd_fsa);

% return an arbitrary pulse sequence to full synthetic aperture (FSA) data
% 'us_foc.seq' must contain the matching pulse sequence for chd_foc
chd_rfc = refocus(us_foc, chd_foc);

% -------------- sampling ------------------ %

% get sampling delay tensors
% NOTE: we assume chd.ord == 'TNM', which means time(T) x receives(N) x transmits(M)
tau_rx = zeros([64, chd.N,   1  ]); % delay matrix: pixels x receives x     1
tau_tx = zeros([64,   1  , chd.M]); % delay matrix: pixels x     1    x transmits
tau = tau_tx + tau_rx;              % delay tensor: pixels x receives x transmits

% sample the data applying transmit delays only, for all transmits
y = sample(chd, tau_tx);

% sample the data applying receive delays only, for all receives
y = sample(chd, tau_rx);

% sample the data using a full sized delay tensor
y = sample(chd, tau); 

% sample the data using separable transmit and receive delay tensors
% (more memory efficient)
y = sample2sep(chd, tau_tx, tau_rx);

% -------- signal processing ------------ %

% use the hilbert transform to get the analytical signal
chd = hilbert(chd);

% 0-pad the data with 128 samples to reduce temporal FFT artefacts
chd = zeropad(chd, 128);

% compute the fft of the data, and the corresponding fft axes
% NOTE: this function does NOT modify the time axes
chdf = fft(chd);
fc = chdf.fftaxis;

% downmix (demodulate) the data to a baseband frequency
% NOTE: the modulation frequency must be passed to the beamformer
fmod = us.xdc.fc; % modulation frequency
chd = downmix(chd, fmod);

% downsample the data by using every 3rd sample
chd = downsample(chd, 3);

% resample the data at a new frequency (25 MHz)
chd = resample(chd, 25e6);

% create a simple passband filter
D = chd.getPassbandFilter(us.xdc.bw);

% filter the ChannelData
chd = filter(chd, D);

% Create a custom digitalFilter directly
% Note: running 'designfilt' without arguments launches a GUI to construct
% a filter interactively
D = designfilt('bandpassfir', ...
    'SampleRate',chd.fs, ...
    'StopbandFrequency1',2e6, ...
    'PassbandFrequency1',3e6, ...
    'PassbandFrequency2',7e6, ...
    'StopbandFrequency2',8e6, ...
    'StopbandAttenuation1',60, ...
    'PassbandRipple',1, ...
    'StopbandAttenuation2',60 ...
    );

% ------------ data typing --------------- %

% change precision to single
% NOTE: double and integer types are supported as well
chd = singleT(chd); 

% extract the datacube as a native MATLAB type
x   = single(chd);

% send data to the GPU (incompatible with tall arrays)
chd = ChannelData();
if gpuDeviceCount, chd = gpuArray(chd); end

% make data a tall type (incompatible with gpuArrays)
chd = ChannelData();
chd = tall(chd);

%% Display and Inspection
% create some example objects and data for this section
us = UltrasoundSystem();
[scan, seq, xdc] = deal(us.scan, us.seq, us.xdc);
us.seq.pulse = Waveform('t0',-1/xdc.fc, 'tend',1/xdc.fc, 'fun',@(t)sinpi(2*xdc.fc*t));
med = Medium();
scat = Scatterers('pos', [mean(us.scan.x),0,mean(us.scan.z)]');
chd = greens(us, scat);
b = DAS(us, chd, 'keep_tx', true);

% --------- Ultrasound System Definitions ----------- %
% plot the entire system 
plot(us);

% plot the pixel locations as dots
plot(scan, '.');

% plot the pulse sequence foci
plot(seq);

% plot the transducer locations as red crosses
plot(xdc, 'r+');

% plot the surface of the transducer elements
patch(xdc); shading faceted;

% plot the impulse response (affectssimulation only)
plot(xdc.impulse);

% plot the excitation pulse (affects simulation only)
plot(us.seq.pulse);

% --------------- Simulation Media ---------------------- %
% plot scatterers as dots
plot(scat, '.');

% show the simulation medium (sound speed)
imagesc(med, scan);

% show the simulation medium density
imagesc(med, scan, 'props', 'rho');

% ------------------ Images & Channel Data ------------------ %
% show the channel data with linear scaling (when data is real)
imagesc(real(chd));

% show the channel data with log-magnitude scaling (when data is complex)
imagesc(hilbert(chd));

% loop through transmits of channel data
h = imagesc(hilbert(chd)); colormap(h.Parent, 'jet');
animate(chd.data, h, 'loop', false);

% display a b mode image
h = imagesc(us.scan, b); colormap(h.Parent,'gray');

% loop through frames/transmits/receives of images data
animate(b, h, 'loop', false);

% animate multiple plots together
nexttile(); h(1) = imagesc(hilbert(chd)); colormap(h(1).Parent,'jet');
nexttile(); h(2) = imagesc(us.scan, b  ); colormap(h(2).Parent,'gray');
hmv = animate({chd.data, b}, h, 'loop', false);

% save animation to a gif
% NOTE: MATLAB may have a bug causing frame sizing to be inconsistent
frame2gif(hmv, 'tmp.gif');

% save animation to an .avi file (higher resolution)
vobj = VideoWriter('tmp', 'Motion JPEG AVI');
vobj.open();
vobj.writeVideo(hmv);
vobj.close();


%% Waveforms (affects Simulation only)
% Waveforms are used when simulating with an excitation function and/or a 
% transducer impulse response function.
 
% --------------- Creation ---------------- %
% Create a delta dirac function
wv = Waveform.Delta();

% create a 2-cycle 5 MHz sine wave tone burst
fc = 5e6;
wv = Waveform('fun', @(t) sin(2*pi*fc*t), 't0', -1/fc, 'tend', 1/fc);

% create an ideal transducer gaussian impulse response 
wv = us.xdc.xdcImpulse();

% import a sampled signal
t = -2e-6 : 1/50e6 : 2e-6;
x = sin(2*pi*fc*t) + sin(2*pi*1/2*fc*t) + sin(2*pi*1/3*fc*t);
wv = Waveform('t', t, 'samples', x);


% ----------- Signal processing ------------ %
% sample the waveform 
tau = (0 : 1023) .* 0.1e-6;
y = wv.sample(tau);

% convolve waveforms
wv2 = conv(wv, wv);

%% Simulation
% create some example objects and data for this section
us = UltrasoundSystem('fs', single(50e6)); % set output sampling frequency (affects simulation only)
scat = Scatterers('pos', [0;0;30e-3], 'c0', 1500);
med = Medium('c0', scat.c0);
dx = 2^(nextpow2(us.xdc.pitch*1e3) - 2); % (optional) compute < 1/2 pitch spacing
cscan = ScanCartesian('x', 1e-3*(-20 : dx : 20), 'z', 1e-3*(-2 : dx : 43));
us.seq = SequenceRadial('type','PW','angles',[-10 0 10]);

% simple lossless delayed signal
chd = greens(us, scat); % NOTE: setting us.fs to single precision is much faster on the GPU.

% FieldII
chd = calc_scat_all(us, scat);
% or 
chd = calc_scat_multi(us, scat);

% MUST
chd = simus(us, scat);

% k-Wave
chd = kspaceFirstOrder(us, med, cscan);

% run simulations on a local or remote cluster
clu = parcluster(); % your default cluster
[job, rfun] = kspaceFirstOrder(us, med, cscan, 'parenv', clu);
% submit(job);
% wait(job);
% chd = rfun(job);

[job, rfun] = calc_scat_multi(us, scat, 'parenv', clu, 'job', true);
% submit(job);
% wait(job);
% chd = rfun(job);



