% ULTRASOUNDSYSTEM - Complete ultrasound system description
%
% The ULTRASOUNDSYSTEM class is a synthesis class containing the properties
% describing a medical ultrasound system and providing methods to simulate
% channel data or beamform channel data. The complete system is described
% by the transmit and receive Transducer, the transmit Sequence, and the
% Scan defining the region for simulation or beamforming.
%
% Multiple simulators are supported, but must be installed by the user. 
% They include:
% 
% * simus via MUST
% * calc_scat_all via FieldII
% * kspaceFirstOrderND via K-wave
% * fullwaveSim via Fullwave
% 
% 
% See also TRANSDUCER SEQUENCE SCAN TARGET CHANNELDATA

classdef UltrasoundSystem < handle
    
    % objects
    properties
        tx          % Transducer object (transmit)
        rx          % Transducer object (receive)
        sequence    % Sequence object
        scan        % Scan object
    end
    
    % parameters
    properties
        fs = 40e6   % simulation sampling frequency (Hz)
    end
    
    properties(Dependent)
        fc                  % central operating frequency (from the Transducer) (Hz)
        xdc                 % Transducer object (if receive and transmit are identical)
        pulse               % Waveform object (from the Sequence)
    end
    
    properties(Hidden,SetAccess=protected)
        tmp_folder          % temporary folder for compiled binaries
    end
        
    % get/set & constructor
    methods
        % constructor
        function self = UltrasoundSystem(varargin)
            % UltrasoundSystem - Construct an UltrasoundSystem
            %
            % us = UltrasoundSystem(Name,Value,...) constructs an
            % UltrasoundSystem object using name value pairs. It's the
            % recommended method of construction.
            %

            
            % initialize Target / Transducer array
            xdc_args = {};
            for i = 1:2:nargin
                switch lower(varargin{i})
                    case {'width', 'radius', 'height', 'impulse','beampattern',...
                            'pitch','numel','kerf', 'fc',...
                            'bandwidth','fractional_bandwidth'}
                        xdc_args{end+1} = varargin{i}; %#ok<AGROW>
                        xdc_args{end+1} = varargin{i+1}; %#ok<AGROW>
                end
            end
            
            % default receiving and transmitting tranducer
            self.rx = TransducerArray(xdc_args{:});
            self.tx = self.rx;

            % set US arguments
            for i = 1:2:nargin
                switch lower(varargin{i})
                    case 'xdc'
                        if(ischar(varargin{i+1}))
                            switch varargin{i+1}
                                case 'IVUS'
                                    self.rx = TransducerIVUS(xdc_args{:});
                                case 'linear'
                                    self.rx = TransducerArray(xdc_args{:});
                                case 'convex'
                                    self.rx = TransducerConvex(xdc_args{:});
                                otherwise
                                    self.rx = TransducerArray(xdc_args{:});
                            end
                        elseif(isa(varargin{i+1},'Transducer'))
                            self.rx = varargin{i+1};
                        end
                        
                        % set the transmitter and receiver to be identical objects
                        self.tx = self.rx;
                        
                    case 'fc'
                        self.fc = varargin{i+1};
                    case 'fs'
                        self.fs = varargin{i+1};
                end
            end
            
            % set the default pulse
            excitation = @(t)exp(-2j*pi*self.fc*t); % functional form
            P = 2; % 2 periods
            
            % set the default pulse sequence
            self.sequence = Sequence(...
                'type','FSA',...
                'focus', [0;0;0], ...
                'pulse', Waveform('t0', 0, 'tend', P / self.fc, 'fun', excitation), ...
                'numPulse', self.tx.numel ...
                );
            
            % set default resolution
            c = 1540; % assuming sound speed == 1540 m/s
            lambda = c ./ self.fc;
            self.scan = ScanCartesian('res', lambda / 4 * ones([1,3]));
            
            % check for sequence, pulse or scan
            for i = 1:2:nargin
                switch lower(varargin{i})
                    case 'sequence'
                        self.sequence = varargin{i+1};
                    case 'pulse'
                        self.sequence.pulse = varargin{i+1};
                    case 'scan'
                        self.scan = varargin{i+1};
                end
            end

            % check for scan options?
            for i = 1:2:nargin
                switch lower(varargin{i})
                    case 'res'
                        self.scan.res = varargin{i+1};
                    case 'xlim'
                        self.scan.xb  = varargin{i+1};
                    case 'ylim'
                        self.scan.yb  = varargin{i+1};
                    case 'zlim'
                        self.scan.zb  = varargin{i+1};
                end
            end
            
            % TODO: use a bin folder locally if we can create one to avoid
            % recompiling to temp when we can instead use a cached version
            % get a temp folder for binaries or code that needs to be
            % recompiled
            self.tmp_folder = tempname; % gives a folder usually in /tmp
            mkdir(self.tmp_folder); % make the folder
            addpath(self.tmp_folder); % this should let us shadow any other binaries

            % copy code or recompile it
            if gpuDeviceCount % only do CUDA stuff if there's a MATLAB-compatible GPU
            defs = self.getCUDAFileDefs();
            fls = arrayfun(@(d) string(strrep(d.Source, 'cu', 'ptx')), defs);
            s = arrayfun(@(fl) copyfile(which(fl), fullfile(self.tmp_folder, fl)), fls);
            if any(~s), self.recompileCUDA(); end % attempt to recompile code
            end
            
            % copy code or recompile it
            % TODO: generalize to mex extension for other machines (with 
            % isunix or iswindows or ismac)
            defs = self.getMexFileDefs();
            fls = arrayfun(@(d) string(strrep(d.Source, 'c', mexext())), defs);
            s = arrayfun(@(fl) copyfile(which(fl), fullfile(self.tmp_folder, fl)), fls);
            if any(~s), self.recompileMex(); end % attempt to recompile code
        end

        function delete(self)
            % DELETE - Destroy an UltrasoundSystem ... programatically.
            %
            % On object destruction, any temporary directories are removed.
            %
            % See also HANDLE

            % if we made a temp folder, clean it up
            if ~isempty(self.tmp_folder) && exist(self.tmp_folder, 'dir')
                rmpath(self.tmp_folder) % remove from the path
                rmdir(self.tmp_folder, 's'); % recursive deletion
            end
        end

        function sub_div = getLambdaSubDiv(self, p, c, ap)
            % GETLAMBDASUBDIV - Get subelement divisions
            % 
            % sub_div = GETLAMBDASUBDIV(self, p, c) returns the element 
            % subdivision vector corresponding to a proportion p of the 
            % wavelength given sound speed c (m/s).
            % 
            % sub_div = GETLAMBDASUBDIV(self, p, c, ap) specifies the
            % aperture ap. Must be one of {'rx'*, 'tx', 'xdc'}
            %
            % 
            
            if(nargin < 4), ap = 'rx'; end
            if isa(c, 'Medium'), c = c.c0; end
            
            % get wavelength
            lam = c / self.(ap).fc;
            
            % get length of proportion of lam
            dp = lam * p;
            
            % get divisions
            sub_div = ceil([self.(ap).width, self.(ap).height] ./ dp);
            
            % ensure number of divisions is odd
            sub_div = sub_div + 1 - rem(sub_div, 2);
        end
    end

    % Modified Green's function based direct computations
    methods
        function [chd, wv] = greens(self, target, element_subdivisions, varargin)
            % GREENS - Simulate ChannelData via a shifted Green's function.
            % 
            % chd = GREENS(self, target)
            % computes the full synthetic aperture data using a simple
            % Green's function kernel applied to all sub elements and all
            % point scatterers. It then applies the transmit sequence to
            % form the channel data.
            %
            % chd = GREENS(self, target, element_subdivisions) uses the 
            % length 2 array element_subdivisions to specifiy the 
            % subdivision of each element into a grid of sub-apertures 
            % in the integration. This argument is passed to FieldII to 
            % construct the subelements. FieldII must be on the path.
            % Setting the subdvisions to [1,1] avoids this behaviour. The
            % Default is [1,1].
            %
            % [chd, wv] = GREENS(...) additionally returns the final
            % waveform convolved across the point scatterers and aperture.
            %
            % [...] = GREENS(..., Name, Value, ...) provides additional
            % name/value pair arguments.
            %
            % Name-Value Arguments:
            %   device -    integer representing gpu selection. -1 selects
            %               the current GPU. 0 selects a cpu. n where n > 0
            %               selects and resets a GPU
            %
            %   interp -    interpolation method; it must be a method 
            %               supported by wsinterpd. The default is 'cubic'.
            %
            %   tall -      specifies whether to use a tall type. Set this
            %               to true if memory is an issue. The default is
            %               false.
            %
            %   bsize -     number of scatterers to compute at a time.
            %               If operating on the GPU when the interp method
            %               is one of {'nearest', 'linear', 'cubic',
            %               'lanczos3'}, this option has no effect. The 
            %               default is 1.
            %
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_MULTI WSINTERPD
            % CHD/SAMPLE
            
            % get Tx/Rx apertures subelement positions
            % (3 x N x E)
            if nargin < 3, element_subdivisions = [1,1]; end
            kwargs.device = -logical(gpuDeviceCount());
            kwargs.interp = 'linear';
            kwargs.tall = false;
            kwargs.bsize = 1;
            kwargs.verbose = false;
            
            % parse inputs
            % TODO: switch to kwargs properties
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end
            
            % get the centers of all the sub-elements
            if all(element_subdivisions == 1) % no sub-elements
                ptc_tx = self.tx.positions();
                ptc_rx = self.rx.positions();
            else % use FieldII's definitions
                ptc_tx = self.tx.getFieldIIBaryCenters(element_subdivisions);
                ptc_rx = self.rx.getFieldIIBaryCenters(element_subdivisions);
            end

            % cast dimensions to compute in parallel
            ptc_rx = permute(ptc_rx, [6,1,2,4,3,5]); % 1 x 3 x N x 1 x En x 1
            ptc_tx = permute(ptc_tx, [6,1,4,2,5,3]); % 1 x 3 x 1 x M x  1 x Em
            
            % get maximum necessary time sample (use manhattan distance and
            % sound speed to give an upper bound)
            maxdist = @(p) max(vecnorm(p,2,1), [], 'all');
            mindist = @(p) min(vecnorm(p,2,1), [], 'all');
            taumax = arrayfun(@(target)(2 * maxdist(target.pos) + maxdist(ptc_tx) + maxdist(ptc_rx)) ./ min(target.c0), shiftdim(target(:),-3));
            taumin = arrayfun(@(target)(2 * mindist(target.pos) + mindist(ptc_tx) + mindist(ptc_rx)) ./ max(target.c0), shiftdim(target(:),-3));
            
            % Directly convolve the Waveform objects to get the final
            % convolved kernel
            wv = conv(self.rx.impulse, ...
                conv(self.tx.impulse, self.sequence.pulse, 4*self.fs), ...
                4*self.fs); % transmit waveform, 4x intermediate convolution sampling
            tk = wv.getSampleTimes(self.fs);
            kern = wv.sample(tk);

            F = numel(target);
            if kwargs.verbose, hw = waitbar(0); end

            for f = F:-1:1 % for each target
            % get minimum/maximum sample times
            tmin = taumin(f) + wv.t0;
            tmax = taumax(f) + wv.tend;

            % create time vector (T x 1)
            % this formulation is guaranteed to pass through t == 0
            t = (floor(tmin * self.fs) : ceil(tmax * self.fs))';

            % pre-allocate output
            [T, N, M, E] = deal(size(t,1), self.rx.numel, self.tx.numel, prod(element_subdivisions));
            x   = complex(zeros([1 T N M]));

            % splice
            c0  = target(f).c0;
            pos = target(f).pos.'; % S x 3
            amp = target(f).amp.'; % S x 1
            fs_ = self.fs;
            if kwargs.device && exist('greens.ptx', 'file') ... % use the GPU kernel
                    && (ismember(kwargs.interp, ["nearest", "linear", "cubic", "lanczos3"]))
                % function to determine type
                isftype = @(x,T) strcmp(class(x), T) || any(arrayfun(@(c)isa(x,c),["tall", "gpuArray"])) && strcmp(classUnderlying(x), T);

                % determine the data type
                if     isftype(kern, 'double'), suffix = "" ;  cfun = @double;
                elseif isftype(kern, 'single'), suffix = "f";  cfun = @single;
                elseif isftype(kern, 'halfT'  ), suffix = "h"; cfun = @(x) alias(halfT(x));
                else,   error("Datatype " + class(kern) + " not recognized as a GPU compatible type.");
                end

                % translate the interp flag
                switch kwargs.interp
                    case "nearest", flagnum = 0;
                    case "linear",  flagnum = 1;
                    case "cubic",   flagnum = 2;
                    case "lanczos3",flagnum = 3;
                    otherwise, error('Interp option not recognized: ' + string(interp));
                end

                % cast data / map inputs
                [x, ps, as, pn, pv, kn, t0k, t0x, fs_, cinv_] = dealfun(cfun, ...
                    x, pos.', amp.', ptc_rx, ptc_tx, kern, t(1)/fs_, tk(1), fs_, 1/c0 ...
                    );

                % re-map sizing
                [QI, QS, QT, QN, QM] = deal(target(f).numScat, length(t), length(kern), N, M);

                % grab the kernel reference
                k = parallel.gpu.CUDAKernel('greens.ptx', 'greens.cu', 'greens' + suffix);
                k.setConstantMemory( ...
                    'QUPS_I', uint64(QI), 'QUPS_T', uint64(QT), 'QUPS_S', uint64(QS), ...
                    'QUPS_N', uint64(QN), 'QUPS_M', uint64(QM) ... , 'QUPS_F', uint64(QF) ...
                    );
                k.ThreadBlockSize = min(k.MaxThreadsPerBlock,QS); 
                k.GridSize = [ceil(QS ./ k.ThreadBlockSize(1)), N, M];

                % call the kernel
                x = k.feval(x, ps, as, pn, pv, kn, t0k, t0x, fs_, cinv_, [E,E], flagnum);

            else % operate in native MATLAB
                % make time in dim 2
                tvec = reshape(t,1,[]); % 1 x T full time vector
                kern_ = reshape(kern,1,[]); % 1 x T' signal kernel

                % TODO: set data types | cast to GPU/tall? | reset GPU?
                if kwargs.device > 0, gpuDevice(kwargs.device); end
                if kwargs.device && ~kwargs.tall,
                                  [pos, ptc_rx, ptc_tx, tvec, kern_] = dealfun(...
                        @gpuArray, pos, ptc_rx, ptc_tx, tvec, kern_ ...
                        );
                elseif kwargs.tall
                              [pos, ptc_rx, ptc_tx, tvec] = dealfun(...
                        @tall, pos, ptc_rx, ptc_tx, tvec ...
                        );
                end

                % for each tx/rx pair, extra subelements ...
                % TODO: parallelize where possible
                svec = num2cell((1:kwargs.bsize)' + (0:kwargs.bsize:target(f).numScat-1), 1);
                svec{end} = svec{end}(svec{end} <= target(f).numScat);
                S = numel(svec); % number of scatterer blocks

                for sv = 1:S, s = svec{sv}; % vector of indices
                    for em = 1:E
                        for en = 1:E
                            % compute time delays
                            % TODO: do this via ray-path propagation through a
                            % medium
                            % S x 1 x N x M x 1 x 1
                            r_rx = vecnorm(sub(pos,s,1) - sub(ptc_rx, en, 5),2,2);
                            r_tx = vecnorm(sub(pos,s,1) - sub(ptc_tx, em, 6),2,2);
                            tau_rx = (r_rx ./ c0); % S x 1 x N x 1 x 1 x 1
                            tau_tx = (r_tx ./ c0); % S x 1 x 1 x M x 1 x 1

                            % compute the attenuation (S x 1 x [1|N] x [1|M] x 1 x 1)
                            att = sub(amp,s,1);% .* (1 ./ r_rx) .* (1 ./ r_tx); % propagation attenuation

                            % get 0-based sample time delay
                            % switch time and scatterer dimension
                            % S x T x N x M x 1 x 1
                            t0 = gather(tk(1)); % 1 x 1

                            % S x T x N x M x 1 x 1
                            if any(cellfun(@istall, {tau_tx, tau_rx, tvec, kern_}))
                                % compute as a tall array  % S | T x N x M x 1 x 1
                                tau = tvec - (tau_tx + tau_rx + t0)*fs_; % create tall ND-array
                                s_ = matlab.tall.transform( ...
                                    @(tau, att) wsinterpd(kern_, tau, 2, att, 1, kwargs.interp, 0), ...
                                    ... @(x) sum(x, 1, 'omitnan'), ... reduce
                                    tau, att ...
                                    );
                                
                                % add contribution (1 x T x N X M)
                                x = x + s_;

                            else
                                % compute natively
                                % TODO: compute where we actually receive a response
                                tau = (tau_tx + tau_rx + t0)*fs_; % S x 1 x N x M x 1 x 1
                                it = tvec == tvec; %(1 x T')
                                it = it & T-1 > tvec - max(tau(:));
                                it = it &   0 < tvec - min(tau(:));

                                % compute only for this range and sum as we go
                                % (1 x T' x N X M)
                                s_ = nan2zero(wsinterpd(kern_, tvec(it) - tau, 2, att, 1, kwargs.interp, 0));

                                % add contribution (1 x T x N X M)
                                x(1,it,:,:,:,:) = x(1,it,:,:,:,:) + s_;
                            end
                            
                            % update waitbar
                            if kwargs.verbose && isvalid(hw), waitbar(sub2ind([E,E,S,F],en,em,sv,F-(f-1)) / (E*E*S*F), hw); end
                        end
                    end
                end
                % compute for tall arrays
                if istall(x), x = gather(x); end

                % move back to GPU if requested
                if kwargs.device, x = gpuArray(gather(x)); end
            end

            % make a channel data object (T x N x M)
            chd(f) = ChannelData('t0', sub(t,1,1) ./ fs_, 'fs', fs_, 'data', shiftdim(x,1));

            % truncate the data if possible
            iszero = all(chd(f).data == 0, 2:ndims(chd(f).data)); % true if 0 for all tx/rx/targs
            n0 = find(cumsum(~iszero, 'forward'), 1, 'first');
            T_ = find(cumsum(~iszero, 'reverse'), 1, 'last' );
            chd(f) = sub(chd(f), n0:T_, chd(f).tdim);


            % synthesize linearly
            [chd(f)] = self.focusTx(chd(f), self.sequence, 'interp', kwargs.interp);
            end

            % close the waitbar if it (still) exists
            if kwargs.verbose && isvalid(hw), delete(hw); end

            % send data (back) to CPU
            % chd = gather(chd);

            % combine all frames
            chd = join(chd, 4);
        end
    end

    % Fullwave calls
    methods
        function conf = fullwaveConf(self, target, sscan, varargin)
            % FULLWAVECONF - Generate a Fullwave simulation configuration
            %
            % conf = FULLWAVECONF(self, targ, scan) creates a simulation
            % configuration struct to be used with fullwaveJob to simulate 
            % the response from the Target targ using the simulation region
            % in the ScanCartesian sscan.
            %
            % See also ULTRASOUNDSYSTEM/FULLWAVEJOB

            % kwarg defaults
            kwargs = struct(...
                'f0', self.xdc.fc, ...    center frequency of the transmit / simulation
                'CFL_max', 0.5, ...       maximum CFL
                'txdel', 'terp' ...     delay models {'disc', 'cont', 'terp'}
                );

            % name-value pairs
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            %% Configuration variables

            % basic vars
            c0       = target.c0;       % speed of sound (m/s)
            omega0   = 2*pi*kwargs.f0;  % center radian frequency of transmitted wave
            dur      = diff(sscan.zb)*2.3/c0; % duration of simulation (s) TODO: make this more general, or an input?

            % define the spatial grid
            % TODO: we can only accept grids where dx/dy/dz exist and are 
            % identical (or infinite/nan) in all dimensions - error if not
            [~, sdims] = ismember('XZ', sscan.order);
            grid.size = sscan.size(sdims); % simulation size
            grid.step = [mode(diff(sscan.x)), mode(diff(sscan.z))]; % simulation step size
            grid.origin = [sscan.x(1), sscan.z(1)]; % first value

            % determine other grid vars
            dX   = min(grid.step);  % limit of spatial step size - will be the same in both dimensions
            fs_  = self.fs;         % data sampling frequency
            cfl0 = (c0*(1/fs_)/dX); % cfl at requested frequency
            modT = ceil(cfl0 / kwargs.CFL_max); % scaling for desired cfl
            dT   = (1/fs_)/modT;    % simulation sampling interval
            cfl  = c0 * dT / dX;    % Courant-Friedrichs-Levi condition
            ppw  = c0 / kwargs.f0 / dX; % points per wavelength - determines grid spacing
            

            % DEBUG 1: this must be consistent with itself!!!
            % cfl,c0,ppw,omega0,dX,dY all go together!
            % at the end: conf.sim  = {c0,omega0,dur,ppw,cfl,maps2,xdcfw,nTic,modT};

            %% Define the Transducer
            xdcfw = self.xdc.getFullwaveTransducer(grid);

            %% Define the Transmit Delays

            % get transmit delays
            tau_tx = -self.sequence.delays(self.tx); % M x V

            % forward to all subelement delays, synchronized
            tau_tx_pix = zeros([xdcfw.nInPx, self.sequence.numPulse]);
            i = logical(xdcfw.incoords(:,4)); % find non-zero indices - they map to an element
            tau_tx_pix(i,:) = tau_tx(xdcfw.incoords(i,4),:); % apply synchronized delays

            
            %% Define the Transmit Apodization

            % get apodization per element
            apod = self.sequence.apodization(self.tx);

            % map to sub-elements
            tx_apod = zeros(size(tau_tx_pix)); % pre-allocate
            i = logical(xdcfw.incoords(:,4)); % find non-zero indices - they map to an element
            tx_apod(i,:) = apod(xdcfw.incoords(i,4),:); % apply synchronized delays

            %% Define the Transmit Pulse
            t0_xdc = min(min(tau_tx_pix,[],1)); % reference true start time (1 x M)
            tau_tx_pix = (tau_tx_pix - t0_xdc); % 0-base the delays for the sim (I x M)
            tau_tx_max = ceil(max(tau_tx_pix, [],'all') / dT) .* dT; % maxmium additional delay
            wv_tx = conv(self.xdc.impulse, self.sequence.pulse, 10/dT); % get the waveform transmitted into the medium
            wv_tx.tend = wv_tx.tend + tau_tx_max; % extend waveform to cover maximum delay
            t = wv_tx.getSampleTimes(1/dT); % get discrete sampling times
            nTic = numel(t); % number of transmit samples
            
            % get transmit signal per input-pixel, time index, transmit for nTic time
            % indices 
            tau_tx_pix = shiftdim(tau_tx_pix,  -1); % 1 x nPxIn x nTx, nTx == M
            switch kwargs.txdel
                % apply discrete shifts in time
                case 'disc',  icmat = cell2mat(arrayfun(@(tau) circshift(wv_tx.sample(t(:)), round(tau ./ dT)), tau_tx_pix, 'UniformOutput', false));
                    % continuous resampling of the waveform
                case 'cont', icmat = cell2mat(arrayfun(@(tau) {wv_tx.sample(t(:) - tau)}, tau_tx_pix));
                case 'terp' % interpolate, upsampling by 10x in time first
                    t_up = wv_tx.getSampleTimes(10*fs_);
                    icmat = interp1(t_up, wv_tx.sample(t_up), t(:) - tau_tx_pix, 'spline', 0); 
            end

            % apply transmit apodization
            icmat = icmat .* shiftdim(tx_apod,-1); % (T' x nInPx x nTx)

            %% Define the Medium
            % get maps in some order
            maps = target.getFullwaveMap(sscan);
            assert(isscalar(sscan.y), 'Unable to simulate in 3D (y-axis is not scalar).');

            % switch to X x Z x Y(=1)
            ord = arrayfun(@(c) find(sscan.order == c), 'XZY');
            for f = string(fieldnames(maps))', maps.(f) = permute(maps.(f), ord); end

            %% Store the configuration variables
            conf.sim = {c0,omega0,dur,ppw,cfl,maps,xdcfw,nTic,modT}; % args for full write function
            conf.tx  = real(icmat);  % transmit data
            conf.t0  = shiftdim(t0_xdc, -1); % time delay offset (1 x 1 x M)
            conf.f0  = kwargs.f0;     % simulation frequency
            conf.fs  = fs_;     % sampling frequency
            conf.tstart = t(1); % signal start time (0 is the peak)
            conf.outmap = xdcfw.outcoords(:,4); % 4th column maps pixels to elements
        end
        
        function [chd, conf] = fullwaveSim(self, target, sscan, varargin)
            % FULLWAVESIM - Simulate channel data via Fullwave
            %
            % chd = FULLWAVESIM(self, target, sscan) simulates the Target 
            % target on the simulation grid sscan and returns a ChannelData
            % object chd. The simulation scan should be large and fine 
            % enough that all elements of the Transducer can be placed.
            %
            % See also ULTRASOUNDSYSTEM/KSPACEFIRSTORDERND

            % defaults
            kwargs = struct(...
                'simdir', fullfile(pwd, 'fwsim'), ...
                'parcluster', parcluster('local') ... parallel cluster
                );

            % parse inputs
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % create the configuration
            conf_args = rmfield(kwargs, setdiff(fieldnames(kwargs), {'f0', 'CFL_max', 'txdel'}));
            conf_args_ = struct2nvpair(conf_args);
            conf = fullwaveConf(self, target, sscan, conf_args_{:});

            % create a job to process it
            job = UltrasoundSystem.fullwaveJob(conf, kwargs.simdir, kwargs.parcluster);

            % submit the job
            submit(job);

            % wait for it to finish
            wait(job);

            % read in the data
            chd = UltrasoundSystem.readFullwaveSim(kwargs.simdir);
        end
    end
    methods(Static)
        function job = fullwaveJob(conf, simdir, clu)
            % FULLWAVEJOB - Create a Fullwave simulation job.
            %
            % job = ULTRASOUNDSYSTEM.FULLWAVEJOB(conf) creates a job to run the
            % fullwave simulation from conf.
            %
            % job = ULTRASOUNDSYSTEM.FULLWAVEJOB(conf, simdir) uses the directory
            % simdir to store the simulation inputs and outputs.
            %
            % job = ULTRASOUNDSYSTEM.FULLWAVEJOB(conf, simdir, clu) creates a job on
            % the parcluster clu. The parcluster must have access to the
            % simulation directory simdir. The resource configurations for
            % the parcluster clu should be setup prior to this call to
            % ensure enough RAM is present.
            %
            % Use 'submit' to submit the job.
            %
            % See also FULLWAVECONF PARCLUSTER PARALLEL.CLUSTER

            % parse inputs
            if nargin < 3, clu = parcluster('local'); end
            if nargin < 2, simdir = fullfile(pwd, 'fwsim'); end

            % make simulation directory
            mkdir(simdir);

            % Write Simulation Files
            write_fullwave_sim(simdir, conf.sim{:}); % all the .dat files
    	    src_folder = fileparts(mfilename('fullpath'));
            copyfile(fullfile(src_folder, '..', 'bin', 'fullwave2_executable'), simdir); % the executable

            % create a job
            % we'll call like this: runFullwaveTx(icmat, simdir, outdir)
            N = size(conf.tx, 3); % number of transmits
            for n = N:-1:1, args{n} = {conf.tx(:,:,n), simdir, fullfile(simdir, string(n))}; end % number the outputs
            job = clu.createJob("AutoAddClientPath",true,"AutoAttachFiles",true);
            job.createTask(@runFullwaveTx, 0, args, "CaptureDiary",true,"Name",'Fullwave-Simulation');

            % write the configuration variables we need for post-processing
            conf = rmfield(conf, ["sim", "tx"]); % remove large variables, save the rest
            save(fullfile(simdir, 'conf.mat'), "-struct", "conf");

    	end

        function chd = readFullwaveSim(simdir, conf)
            % READFULLWAVESIM - Create a ChannelData object from the fullwave simulation.
            %
            % chd = ULTRASOUNDSYSTEM.READFULLWAVESIM(simdir) creates a ChannelData object from the
            % simulation files located in the directory simdir.
            %
            % chd = ULTRASOUNDSYSTEM.READFULLWAVESIM(simdir, conf) uses the configuration file conf
            % rather than loading one from the simulation directory.
            %
            % See also RUNFULLWAVETX ULTRASOUNDSYSTEM/FULLWAVEJOB

            % see if we can find a conf file in the directory if not given to us
            if nargin < 2
                cfile = fullfile(simdir, 'conf.mat');
                if exist(cfile, 'file'), conf = load(cfile, 'fs', 't0', 'tstart', 'outmap');
                else, error('Unable to find configuration file.'); end
            end

            % get the folders in the simulation directory
            % they are numbered for each transmit, so we can just count them
            listing = dir(simdir);
            dnm = double(string({listing.name})); % directory names
            mind = dnm(~isnan(dnm));
            M = max(mind); % number of transmits - inferred from number of files

            % get output mapping
            outmap = conf.outmap;
            mPx = length(outmap);

            % number of elements
            N = max(outmap);

            % preallocate
            xm = cell(N,M);
            f = waitbar(0, 'Initializing ...');
            for m = 1:M
                % update status
                if isvalid(f), waitbar(m/M, f, sprintf('Loading Files: %d %%', floor(m/M*100))); end

                % load data (with a reference, in case it is large)
                outfile = fullfile(simdir, num2str(m), 'genout.dat');
                mdat = memmapfile(outfile, 'Writable', false, 'Format', 'single');
                xm_ = reshape(mdat.Data, mPx, []); % N'' x T'' % pixels x time steps

                % average across sub-elements to get signal for a single receive element
                for n = 1:N, xm{n,m} = mean(xm_((outmap == n),:),1)'; end % T x {N} x {M}
            end
            % set all values to have the same size
            xm = reshape(cat(1, xm{:}), [], N, M); % make datacube

            % Construct a Channel Data object
            chd = ChannelData('data', xm, 'fs', conf.fs, 't0', conf.t0 + conf.tstart);

            % cleanup
            if isvalid(f), delete(f); end
        end
    end

    % SIMUS calls
    methods
        function chd = simus(self, target, varargin)
            % SIMUS - Simulate channel data via MUST
            %
            % chd = SIMUS(self, target) simulates the Target target and
            % returns a ChannelData object chd. By default it runs on the
            % current parpool object returned by gcp
            %
            % chd = SIMUS(...,'dims', D, ...) selects the number of 
            % dimensions for the simulation. D must be one of {2, 3, []*}. 
            % If D is empty, the dimensions are chosen based on the point 
            % targets.
            %
            % chd = SIMUS(...,'periods', T, ...) selects the number of
            % periods of the tone burst. 
            %
            % chd = SIMUS(...,'interp', method, ...) selects the method of
            % interpolation to use when synthesizing transmits from the
            % full-synthetic-aperture data.
            %
            % chd = SIMUS(...,'parcluster', clu, ...) selects the
            % parcluster for parfor to run on for each transmit. It must 
            % either be a type of parcluster, parpool, or 0 for no 
            % parcluster. The default is the parpool returned by gcp. 
            %
            % chd = SIMUS(...,'device', d, ...) selects the device to run
            % on. 
            %
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_ALL
            % ULTRASOUNDSYSTEM/FOCUSTX
            
            % defaults
            % TODO: forward arguments to params or opt as appropriate
            kwargs = struct(...
                'interp', 'cubic', ...
                'parcluster', gcp('nocreate'), ...
                'periods', 1, ...
                'dims', [] ...
                );

            % load options
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end
            if isempty(kwargs.parcluster), kwargs.parcluster = 0; end % select 0 workers if empty

            % TODO: check the transmit/receive/sequence impulse: they 
            % cannot be satisfied if not a Delta or empty
            if ~ismember("periods", varargin(cellfun(@ischar,varargin) | cellfun(@isstring, varargin))) % was this an input?
                warning("QUPS:UltrasoundSystem:simus:unsatisfiable", "Transmit sequence determined by 'periods', property.");
            end
            
            % get the points and the dimensions of the simulation(s)
            [X, Y, Z, A] = arrayfun(@(target) ...
                deal(sub(target.pos,1,1), sub(target.pos,2,1), sub(target.pos,3,1), target.amp), ...
                target, 'UniformOutput',false);
            if isempty(kwargs.dims) 
                if all(cellfun(@(Y)Y == 0, 'all'),Y), kwargs.dims = 2; [Y{:}] = deal([]); % don't simulate in Y if it is all zeros 
                else, kwargs.dims = 3; end
            end
            if kwargs.dims == 2 && cellfun(@(Y)any(Y ~= 0, 'all'),Y)
                warning("QUPS:UltrasoundSystem:simus:casting", "Projecting all points onto Y == 0 for a 2D simulation.");
            end

            % get all other param struct values (implicitly force same
            % transducer)
            p = {getSIMUSParam(self.xdc)};
            p = cellfun(@struct2nvpair, p, 'UniformOutput', false);
            p = cat(2, p{:});
            p = struct(p{:});

            % set transmit sequence ... the only way we can
            % TODO: forward arguments to transmit parameters
            p.fs    = self.fs;
            p.TXnow = kwargs.periods; % number of wavelengths
            p.TXapodization = zeros([self.xdc.numel,1]); % set tx apodization
            p.RXdelay = zeros([self.xdc.numel,1]); % receive delays (none)
            
            % set options per target
            p = repmat(p, [1,1,1,numel(target)]); % (1 x 1 x 1 x F)
            pxdc = arrayfun(@(target) {getSIMUSParam(target)}, target); % properties per target
            for f = 1:numel(target), 
                for fn = string(fieldnames(pxdc{f}))', 
                    p(f).(fn) = pxdc{f}.(fn); % join the structs manually
                end 
            end

            % set options 
            % TODO: forward Name-Value pair arguments
            opt = struct( ...
                'ParPool', false, ... % parpool on pfield.m
                'FullFrequencyDirectivity', false, ... % use central freq as reference
                'ElementSplitting', 1, ... % element subdivisions
                'WaitBar', false, ... % add wait bar
                'dBThresh', -100 ... % threshold for computing each frequency
                ... 'FrequencyStep', df, ... % freuqency domain resolution
                ... 'CallFun', 'simus' ... % hack: use the simulation portion of the code
                );

            % select the computing cluster
            clu = kwargs.parcluster;
            if isempty(clu), clu = 0; end % empty pool -> 0
            isloc = ~isa(clu, 'parallel.Pool') || ~isa(clu, 'parallel.Cluster'); % local or parpool
            if isloc, [pclu, clu] = deal(clu, 0); else, pclu = 0; end % cluster or local

            % call the sim: FSA approach
            [M, F] = deal(self.xdc.numel, numel(target)); % splice
            for f = F:-1:1 % per target
                argf = {X{f},Y{f},Z{f},A{f},zeros([M,1]),p(f),opt}; % args per target
                parfor (m = 1:M, pclu) % use parallel rules, but execute on main thread
                    args = argf; % copy settings for this frame
                    args{6}.TXapodization(m) = 1; % transmit only on element m
                    if isloc, rf{m,f} = simus(args{:}); % local compute
                    else, out(m,f) = parfeval(clu, @simus, 1, args{:}); % add cluster job
                    end
                end
            end
            
            % gather outputs when finished
            if ~isloc, rf = fetchOutputs(out, "UniformOutput",false); end

            % create the output QUPS ChannelData object
            chd = cellfun(@(x) ChannelData('data', x), rf); % per transmit/frame (M x F) object array
            chd = arrayfun(@(f) join(chd(:,f), 3), 1:F); % join over transmits (1 x F) object array
            chd = join(chd, 4); % join over frames (1 x 1) object array
            chd.fs = self.fs; % set the smapling frequency
            chd.t0 = 0; % already offset within call to simus

            % synthesize transmit pulses
            chd = self.focusTx(chd, self.sequence, 'interp', kwargs.interp);
        end
    end

    % Field II calls
    methods(Access=public)
        function chd = calc_scat_all(self, target, element_subdivisions, varargin)
            % CALC_SCAT_ALL - Simulate channel data via FieldII
            %
            % chd = CALC_SCAT_ALL(self, target) simulates the Target target
            % and returns a ChannelData object chd.
            % 
            % chd = CALC_SCAT_ALL(self, target, element_subdivisions)
            % specifies the number of subdivisions in width and height for 
            % each element.
            %
            % chd = CALC_SCAT_ALL(..., 'interp', method) specifies the
            % interpolation methods for the transmit synthesis. The method
            % must be supported by focusTx.
            %
            % See also ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/FOCUSTX
            
            % helper function
            vec = @(x) x(:); % column-vector helper function

            % defaults
            kwargs = struct('parcluster', gcp('nocreate'), 'interp', 'linear');
            if nargin < 3, element_subdivisions = self.getLambdaSubDiv(0.1, target); end

            % load options
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = self.tx.impulse; % transmitter impulse
            wv_rx = self.rx.impulse; % receiver impulse
            wv_pl = self.sequence.pulse;
            
            % get the time axis (which passes through t == 0)
            t_tx = wv_tx.getSampleTimes(self.fs);
            t_rx = wv_rx.getSampleTimes(self.fs);
            t_pl = wv_pl.getSampleTimes(self.fs);

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.fun(t_tx))')));
            rx_imp = gather(double(real(vec(wv_rx.fun(t_rx))')));
            tx_pls = gather(double(real(vec(wv_pl.fun(t_pl))')));

            % choose the cluster to operate on: avoid running on ThreadPools
            clu = kwargs.parcluster;
            if isempty(clu) || isa(clu, 'parallel.ThreadPool') || isa(clu, 'parallel.BackgroundPool'), clu = 0; end

            % splice
            [M, F] = deal(self.sequence.numPulse, numel(target)); % number of transmits/frames
            [fs_, tx_, rx_] = deal(self.fs, self.tx, self.rx); % splice
            [c0, pos, amp] = arrayfun(@(t)deal(t.c0, {t.pos}, {t.amp}), target); % splice

            % Make position/amplitude and transducers constants across the workers
            if isa(clu, 'parallel.Pool'), 
                cfun = @parallel.pool.Constant;
            else, 
                cfun = @(x)struct('Value', x);
                [pos, amp] = deal({pos},{amp}); % for struct to work on cell arrays
            end
            [pos_, amp_, tx_, rx_] = dealfun(cfun, pos, amp, tx_, rx_);

            % for each target (frame)
            parfor (f = 1:F, clu) % each transmit pulse
            % reinitialize field II
            field_init(-1);
            
            % set sound speed/sampling rate
            set_field('fs', fs_);
            set_field('c',c0(f));

            % get Tx/Rx apertures
            p_focal = [0;0;0];
            Tx = tx_.Value.getFieldIIAperture(p_focal.', element_subdivisions); %#ok<PFBNS> % constant over workers
            Rx = rx_.Value.getFieldIIAperture(p_focal.', element_subdivisions); %#ok<PFBNS> % constant over workers

            % set the impulse response function / excitation function
            xdc_impulse   (Tx, tx_imp);
            xdc_impulse   (Rx, rx_imp);
            xdc_excitation(Tx, tx_pls);

            
            % call the sim
            down_sampling_factor = 1;
            [v, ts(1,f)] = calc_scat_all(Tx, Rx, ...
                 pos_.Value{f}.', amp_.Value{f}.', down_sampling_factor); %#ok<PFBNS> % constant over workers
            
            % reshape to T x N x M
            voltages{1,f} = reshape(v, [], rx_.Value.numel, tx_.Value.numel);
            end

            % create the output QUPS ChannelData object
            chd = cellfun(@(x) ChannelData('data', x), voltages);
            chd = join(chd, 4);
            chd.t0 = shiftdim(ts,-2) + t_pl(1) + t_tx(1) + t_rx(1);
            chd.fs = fs_;

            % synthesize linearly
            chd = self.focusTx(chd, self.sequence, 'interp', kwargs.interp);
        end        

        function chd = calc_scat_multi(self, target, element_subdivisions, varargin)
            % CALC_SCAT_MULTI - Simulate channel data via FieldII
            %
            % chd = CALC_SCAT_MULTI(self, target) simulates the Target target
            % and returns a ChannelData object chd.
            %
            % chd = CALC_SCAT_MULTI(self, target, element_subdivisions)
            % specifies the number of subdivisions in width and height for
            % each element.
            %
            % See also ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/FOCUSTX

            % helper function
            vec = @(x) x(:); % column-vector helper function

            % defaults
            kwargs = struct('parcluster', gcp('nocreate'));

            % load options
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % initialize field II
            try
                evalc('field_info');
                field_started = false;
            catch
                field_init(-1);
                field_started = true;
            end

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = self.tx.impulse; % transmitter impulse
            wv_rx = self.rx.impulse; % receiver impulse
            wv_pl = self.sequence.pulse;

            % get the time axis (which passes through t == 0)
            t_tx = wv_tx.getSampleTimes(self.fs);
            t_rx = wv_rx.getSampleTimes(self.fs);
            t_pl = wv_pl.getSampleTimes(self.fs);

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.fun(t_tx))')));
            rx_imp = gather(double(real(vec(wv_rx.fun(t_rx))')));
            tx_pls = gather(double(real(vec(wv_pl.fun(t_pl))')));

            % get the apodization and time delays across the aperture
            apod_tx = self.sequence.apodization(self.tx); % N x M
            tau_tx = -self.sequence.delays(self.tx); % N x M
            tau_offset = min(tau_tx, [], 1); % (1 x M)
            tau_tx = tau_tx - tau_offset; % 0-base the delays for FieldII

            % choose the cluster to operate on: avoid running on ThreadPools
            clu = kwargs.parcluster;
            if isempty(clu) || isa(clu, 'parallel.ThreadPool') || isa(clu, 'parallel.BackgroundPool'), clu = 0; end

            [M, F] = deal(self.sequence.numPulse, numel(target)); % number of transmits/frames
            [fs_, tx_, rx_] = deal(self.fs, self.tx, self.rx); % splice
            [c0, pos, amp] = arrayfun(@(t)deal(t.c0, {t.pos}, {t.amp}), target); % splice
            
            % Make position/amplitude and transducers constants across the workers
            if isa(clu, 'parallel.Pool'), 
                cfun = @parallel.pool.Constant;
            else, 
                cfun = @(x)struct('Value', x);
                [pos, amp] = deal({pos},{amp}); % for struct to work on cell arrays
            end
            [pos_, amp_, tx_, rx_] = dealfun(cfun, pos, amp, tx_, rx_);

            parfor (m = 1:M, clu) % each transmit pulse
            for (f = F:-1:1) % each target frame            
                % (re)initialize field II
                field_init(-1);

                % set sound speed/sampling rate
                set_field('fs', fs_);
                set_field('c', c0(f)); %#ok<PFBNS> % this array is small

                % get Tx/Rx apertures
                p_focal = [0;0;0];
                Tx = tx_.Value.getFieldIIAperture(p_focal.', element_subdivisions); %#ok<PFBNS> % constant over workers
                Rx = rx_.Value.getFieldIIAperture(p_focal.', element_subdivisions); %#ok<PFBNS> % constant over workers

                % set the impulse response function / excitation function
                xdc_impulse   (Tx, tx_imp);
                xdc_impulse   (Rx, rx_imp);
                xdc_excitation(Tx, tx_pls);

                % nullify the response on the receive aperture
                xdc_times_focus(Rx, 0,  zeros([1, rx_.Value.numel])); % set the delays manually

                % for each transmit, set the transmit time delays and
                % apodization
                xdc_times_focus(Tx, 0,  tau_tx(:,m)'); % set the delays
                xdc_apodization(Tx, 0, apod_tx(:,m)'); % set the apodization
                
                % call the sim
                [voltages{m,f}, ts{1,1,m,f}] = calc_scat_multi(Tx, Rx, pos_.Value{f}.', amp_.Value{f}.'); %#ok<PFBNS> % constant over workers
            end
            end
            
            % adjust start time based on signal time definitions
            t0 = cell2mat(ts) + ... % fieldII start time (1 x 1 x M x F)
                (t_pl(1) + t_tx(1) + t_rx(1)) ... signal delays for impulse/excitation
                + shiftdim(tau_offset,-1) ... 0-basing the delays across the aperture
                ; % 1 x 1 x M
            
            % create the output QUPS ChannelData object 
            chd = cellfun(@(x) ChannelData('data', x), voltages); % per transmit/frame (M x F) object array
            chd = arrayfun(@(f) join(chd(:,f), 3), 1:F); % join over transmits (1 x F) object array
            chd = join(chd, 4); % join over frames (1 x 1) object array

            % set sampling frequency and transmit times for all
            chd.fs = self.fs;
            chd.t0 = t0;

            % cleanup
            if field_started, evalc('field_end'); end
        end
    end
    
    % k-Wave calls (old)
    methods(Access=public)
        function [kgrid, PML_size, kgrid_origin, kgrid_size, kgrid_step, kmedium] ...
                = getkWaveGrid(self, target, varargin)
            % GETKWAVEGRID - Create a kWaveGrid for the Target
            % 
            % [kgrid, PML_size, kgrid_origin, kgrid_size, kgrid_step, 
            % kmedium] = GETKWAVEGRID(self, target) creates a kWaveGrid
            % kgrid for the given Target target and also returns the 
            % selected sizing for the perfectly matched layers (PML) and 
            % the corresponding sizing and offset for the absolute
            % coordinate system.
            %
            % This signature (inputs/outputs) is likely to change.
            %
            % Inputs:
            %   target -       a Target object
            %
            % Name-Value Pair Inputs:
            %   dims -         number of dimensions (default = 2)
            %
            %   PML_min -      minimum PML size (default = 4)
            %
            %   PML_max -      minimum PML size (default = 48)
            %
            %   CLF_max -      maximum CFL: the kgrid will use a physical
            %                   and temporal spacing beneath this value
            %                   that is a ratio of the sampling frequency
            %                   (default = 0.25)
            %
            %   buffer(3 x 2) - 3 x 2 array specifying computational 
            %                   buffer to provide spacing from the other 
            %                   objects including the transmitter, 
            %                   receiver, and target objects. This is the 
            %                   amount of room that the grid is expanded 
            %                   by. Given in meters as 
            %                   [xlo, xhi; ylo, yhi, zlo, zhi] 
            %                   (default = repmat(5e-3, [3,2]) )
            %
            %   resolution-ratio - sets the ratio of grid spacing to
            %                   spatial wavelength assuming a given a 
            %                   reference speed of sound from the target
            %                   (default = 0.25)
            %
            %   reference-sound-speed - overrides the reference sound
            %                   speed for the resolution ratio. 
            %                   (default = target.c0)
            %
            % Outputs:
            %   kgrid -            a kWaveGrid object
            %
            %   PML_size (3 x 1) - the chosen PML size that minimizes the
            %                       maximum prime number within PML-min and
            %                       PML-max
            %
            %   kgrid_origin (3 x 1) - cartesian grid origin to recover the
            %                       original coordinates when using
            %                       kgrid.{x|y|z} or similar functions
            %
            %   kgrid_size (3 x 1) - size of the computational grid, 
            %                       without the PML. The size is 1 for
            %                       sliced dimensions
            %
            %   kgrid_step (3 x 1) - descritization size in each dimension
            %                       in meters. The size is inf for sliced
            %                       dimensions
            %   kmedium            - kWave compatible medium structure
            %
            % 
            
            % defaults 
            lam_res = 0.25; % default resolution as proportion of wavelength
            PML_min_size = int64(4); % minimum (one-sided) PML size
            PML_max_size = int64(48); % maximum (one-sided) PML size
            max_cfl = 0.25; % maximum cfl number (for stability)
            cbuf = repmat(5e-3, [3,2]); % x/y/z computational buffer
            c0 = target.c0; % get reference sound speed
            dims = 2; % default simulation dimensions
            pb_u = zeros([3,0]); % default user provided bounds
            
            % override defaults with name-value pairs
            % TODO: switch to kwargs properties
            for i = 1:2:nargin-2
                switch varargin{i}
                    case 'dims'
                        dims = varargin{i+1};
                    case 'PML_min'
                        PML_min_size = int64(varargin{i+1});
                    case 'PML_max'
                        PML_max_size = int64(varargin{i+1});
                    case 'CFL_max'
                        max_cfl = varargin{i+1};
                    case 'buffer'
                        cbuf = varargin{i+1};
                    case 'bounds'
                        pb_u = varargin{i+1};
                    case 'resolution_ratio'
                        lam_res = varargin{i+1};
                    case 'reference_sound_speed'
                        c0 = varargin{i+1};
                end
            end
                        
            % operate at ratio of lambda min in x/y/z
            lam = c0 ./ max(max(self.rx.bw, self.tx.bw));
            [dx, dy, dz] = deal(lam_res * min(lam));
            dp = [dx;dy;dz]; % vector difference
            
            % finds the smallest PML buffer from PML_min_size to search_max that
            % minimizes the factoring size
            PML_buf = @(n) (argmin(arrayfun(@(a)max(factor(n+2*a)), PML_min_size:PML_max_size)) + PML_min_size - 1);
            
            % get total min/max bounds for the transducers and the target
            pb_t = target.getBounds(); % get min/max bounds of the target
            pb_tx = self.tx.bounds(); % get min/max bounds of the tx
            pb_rx = self.rx.bounds(); % get min/max bounds of the rx
            pb_all = cat(2, pb_t, pb_tx, pb_rx, pb_u); % min/max bounds for all objects
            pb = [min(pb_all,[],2), max(pb_all,[],2)]; % reduce
            
            % get kgrid and outer PML sizing with computational buffer
            pb_buf = pb + [-1, 1] .* cbuf;
            Np_g = ceil(diff(pb_buf,1,2) ./ dp); % minimum grid size
            Np_g = Np_g + (1-mod(Np_g,2)); % enforce odd grid size
            Np_p = double(arrayfun(PML_buf, Np_g)); % get complimentary PML buffer with good FFT sizing
                        
            % get cartesian grid origin to interface with k-wave's centered grid
            p_kw_og = min(pb_buf,[],2) + dp .* ((Np_g - 1) ./ 2); % k-wave origin
            
            % translate into additional outputs in kwave format
            ind_map = [3,1,2]; % UltrasoundSystem to k-wave coordinate mapping
            PML_size = Np_p(ind_map); % PML size
            kgrid_origin = p_kw_og(ind_map); % grid origin
            kgrid_size = Np_g(ind_map); % grid size - size 1 in sliced dimensions
            kgrid_step = dp(ind_map); % step size - inf in sliced dimensions

            % apply k-wave <-> UltrasoundSystem dimension mapping: x/y/z -> z/x/y
            switch dims
                case 1
                    kgrid = kWaveGrid(kgrid_size(1), kgrid_step(1));
                    [kgrid_size(2:3), PML_size(2:3), kgrid_step(2:3), kgrid_origin(2:3)] = deal(1, 0, inf, 0);
                case 2
                    kgrid = kWaveGrid(kgrid_size(1), kgrid_step(1), kgrid_size(2), kgrid_step(2));
                    [kgrid_size(3), PML_size(3), kgrid_step(3), kgrid_origin(3)] = deal(1, 0, inf, 0);
                case 3
                    kgrid = kWaveGrid(kgrid_size(1), kgrid_step(1), kgrid_size(2), kgrid_step(2), kgrid_size(3), kgrid_step(3));
            end
            
            
            % get sound speed bounds
            kmedium = target.getKWaveMedium(kgrid, kgrid_origin);
            c = kmedium.sound_speed;
            [c_min, c_max] = deal(min(c,[],'all'), max(c, [], 'all'));
            
            % get the sampling interval as the largest value within the cfl
            dt_cfl_max = max_cfl * min(kgrid_step) / c_max;
            
            % use a time step that aligns with the sampling frequency
            dt_us = inv(self.fs);
            cfl_ratio = dt_cfl_max / dt_us;
            if cfl_ratio  >= 1
                time_step_ratio = inv(floor(cfl_ratio));
            else
                warning('Upsampling the kWave time step for an acceptable CFL.');
                time_step_ratio = ceil(inv(cfl_ratio));
            end
            dt = dt_us / time_step_ratio;
            
            % get the source signal and the Tx impulse response function
            tx_impulse_waveform = self.tx.impulse();
            
            % get temporal buffer as twice the longest signal
            t_buf = 2* max(...
                abs(tx_impulse_waveform.tend - tx_impulse_waveform.t0), ...
                abs(self.sequence.pulse.tend - self.sequence.pulse.t0)...
                );

            % get a default ending time based on the longest possible 
            % single-way reflection
            tend = t_buf + 2 * norm([kgrid.x_size, kgrid.y_size, kgrid.z_size]) / c_min;
            
            % define the time axis
            kgrid.setTime(ceil(tend/dt), dt);            
            
        end        
        
        function [ksource, tsig, sig0] = getKWaveSource(self, kgrid, kgrid_origin, el_sub_div, c0)
            % GETKWAVESOURCE - Create a k-wave compatible source structures
            % 
            % [ksource, tsig, sig0] = GETKWAVESOURCE(self, kgrid, kgrid_origin, el_sub_div, c0)
            % creates k-Wave compatible source structures.
            % 
            % See also ULTRASOUNDSYSTEM/GETKWAVEGRID
            % ULTRASOUNDSYSTEM/GETKWAVERECEIVE

            %
            if nargin < 4 || isempty(el_sub_div), el_sub_div = [1,1]; end
            if nargin < 5 || isempty(c0)        , c0 = 1540; end
            
            nPulse = self.sequence.numPulse; % number of pulses
            if(isnan(nPulse)), nPulse = self.tx.numel; end % for FSA apertures

            % helper function
            vec = @(x) x(:);
            
            % beamforming
            sigt0 = self.sequence.t0Offset(); % 1 x S
            steering_delays = -sigt0 - self.sequence.delays(self.tx); % (N x S)
            el_apodization  = self.sequence.apodization(self.tx); % ([1|N] x [1|S])
            el_apodization = el_apodization + zeros([self.tx.numel nPulse]); % (N x S)
            
            % get the frequency upsampling ratio
            % kwave_upsampling_ratio = round(inv(kgrid.dt * self.fs));

            % define the transmitted signal
            % get the source signal and the Tx impulse response function
            tx_impulse_waveform = self.tx.impulse;
            chirp_waveform = arrayfun(@(c)copy(c), self.sequence.pulse); % copy the waveforms ([1|N] x [1|S])
            
            % choose the longest signal in time for the time-domain
            imp_dur = tx_impulse_waveform.tend - tx_impulse_waveform.t0;
            sig_dur = reshape([chirp_waveform.tend] - [chirp_waveform.t0], size(chirp_waveform));
            
            % get the largest subelement-surface to pixel step size
            grid_delay = hypot(hypot(kgrid.dx, kgrid.dy), kgrid.dz) / c0;
            
            % check for delta functions - functions less than a max delay
            sig_is_delta = abs(sig_dur) < grid_delay;
            imp_is_delta = abs(imp_dur) < grid_delay;
            
            % get minimum and maximum steering angle delays
            [min_steer, max_steer] = deal(min(steering_delays(:)), max(steering_delays(:)));
            
            % get max cumulative delay to ensure that signal is sampled
            % starting from before the first waveform begins until after
            % the last waveform ends
            tdur = imp_dur + sig_dur + grid_delay;
            
            % get the signal time domain 
            t0   = tx_impulse_waveform.t0   + min([chirp_waveform.t0])   - tdur + min_steer;
            tend = tx_impulse_waveform.tend + max([chirp_waveform.tend]) + tdur + max_steer;
            [n0, nend] = deal(floor(t0 / kgrid.dt), ceil(tend / kgrid.dt));
            t  = (  n0   : 1 :   nend  )' * kgrid.dt; % includes zero if t0 <= 0 <= tend
            tc = (2*n0-1 : 1 : 2*nend-1)' * kgrid.dt; % includes zero if t0 <= 0 <= tend
            tsig = tc + shiftdim(sigt0(:),-2); % offset output time zero (T x 1 x [1|S])

            % tests
            assert(min(t) <= 0 && max(t) >= 0, 'The sampled signal does not begin and end after 0.'); % test that we are actually sampling through t == 0
            assert(any(t==0), 'The sampled signal does not pass through 0.'); % test we pass through 0.
            
            % define the sampling functional: 
            % always the same length as t
            % t always includes 0
            function sig = shiftSignal(delay,ind)
                % check for delta functions - functions less than a max delay
                sig_is_delta = abs(sig_dur(ind)) < grid_delay; 
                imp_is_delta = abs(imp_dur) < grid_delay;
                
                if(~imp_is_delta) % impulse not a delta: delay the impulse
                    sig_samp = vec(tx_impulse_waveform(ind).sample(t - delay));
                    imp_samp = vec(chirp_waveform.sample(t));
                    sig = conv(sig_samp, imp_samp, 'full');
                elseif(~sig_is_delta) % signal is not a delta: delay on signal
                    sig_samp = vec(tx_impulse_waveform(ind).sample(t));
                    imp_samp = vec(chirp_waveform.sample(t - delay));
                    sig = conv(sig_samp, imp_samp, 'full');
                else % both are deltas use linear resampling (sketchy)
                    sig_samp = vec(tx_impulse_waveform(ind).sample(t));
                    imp_samp = vec(chirp_waveform.sample(t));
                    sig = interp1(t,conv(sig_samp, imp_samp, 'full'), t - delay, 'linear', 0);
                end
                % sig = real(sig) .* inv(norm(vec(abs(sig))));
                % normalize samples?
                nsig = norm(vec(abs(sig)));
                isig = nsig > 0;
                sig(isig) = real(sig(isig)) .* inv(nsig(isig));
            end
            
            % output the non-delayed signals
            for ind = numel(chirp_waveform):-1:1, sig0(:,ind) = shiftSignal(0, ind); end
            sig0 = reshape(sig0, [size(sig0,1), size(chirp_waveform)]);
            
            % get the sensor, sensor masks, and element orientations
            [ksensor, ~, sens_map] = self.tx.getKWaveSensor(kgrid, kgrid_origin, el_sub_div);
            
            % initialize the transmit sources
            parfor puls = 1:nPulse
                % get the apodization elements for this pulse
                pulse_apodization = el_apodization(:, puls);
                el_non_zero = pulse_apodization ~= 0;
                
                % get the mask for this pulse
                pulse_mask = false(size(ksensor{1}.mask)); %#ok<PFBNS>
                for k = 1:numel(el_non_zero), if el_non_zero(k), pulse_mask = pulse_mask | ksensor{k}.mask; end, end
                
                % initialize output source
                ksource{puls} = struct(...
                    'u_mask', pulse_mask ...
                    );
            end
            
            % send constant inputs to the workers
            if(isvalid(gcp('nocreate')))
                sendDataToWorkers = @parallel.pool.Constant;
            else
                sendDataToWorkers = @(s)struct('Value', s);
            end
            [t_vec] = dealfun(sendDataToWorkers, vec(t));
            
            % dereference
            [nTx, T] = deal(self.tx.numel, numel(tc));
            
            % define the transmit sources
            fprintf('\nDefining sources ...\n');
            ticBytes(gcp('nocreate'));
            parfor (puls = 1:nPulse)
                % report performance
                fprintf('\nProcessing %i samples for all %i elements of pulse %i ... \n', T, nTx, puls); tt = tic;
                
                % dereference
                t_ = t_vec.Value; %#ok<PFBNS>
                
                % initialize
                ind_msk_pulse   = find(ksource{puls}.u_mask);
                ksource_vec_vel = zeros([3, numel(ind_msk_pulse), T], 'like', t_);
                
                for el = 1:nTx
                    % steering delay for the element
                    tau_steer = steering_delays(el,puls);

                    % get sensitivity map for this element
                    sens_map_el = sens_map(:,el); %#ok<PFBNS>
                    
                    % get sensitivity mapping for all subelements
                    %{
                    for i = nSub:-1:1
                        [ind_msk, dis, el_dir, amp] = deal(...
                            vec(sens_map_el(i).mask_indices).', ...
                            vec(sens_map_el(i).dist).', ...
                            vec(sens_map_el(i).element_dir),...
                            vec(sens_map_el(i).amp) ...
                            );
                        
                        % vector amplitude
                        vec_amp_i{i} = amp .* el_dir + zeros(size(dis));
                        
                        % corresponding indices
                        ind_msk_i{i} = ind_msk;
                        
                        % get distance from element center to grid points
                        % + the steering delay for the element
                        delays_i{i} = dis ./ c0 + tau_steer;
                    end
                    %}
                    
                    [vec_amp_i, ind_msk_i, delays_i] = arrayfun(@(sme)deal(...
                        vec(sme.amp).' .* (sme.element_dir) + zeros(size(sme.dist(:)')), ...
                        vec(sme.mask_indices).', ...
                        vec(sme.dist).' ./ c0 + tau_steer ...
                        ), sens_map_el, 'UniformOutput', false);
                    
                    % combine
                    [vec_amp, ind_msk, delays] = dealfun(@cell2mat, vec_amp_i(:)', ind_msk_i(:)', delays_i(:)');
                    
                    % filter out zero amplitudes
                    ind_non_zero = any(logical(vec_amp),1); % amplitude zero in all dims
                    if(~any(ind_non_zero)), continue, end % move to next element if everything is zero
                    [vec_amp, ind_msk, delays] = dealfun(@(x)x(:,ind_non_zero), vec_amp, ind_msk, delays);
                    
                    % get mask indices
                    [~, ind_umsk] = ismember(ind_msk, ind_msk_pulse);
                    if(~any(ind_umsk)), continue, end % move to next element if no valid elements
                    
                    % resample with the delay                    
                    %%%
                    if(~imp_is_delta) % impulse not a delta: delay the impulse
                        imp_samp = (tx_impulse_waveform.sample(t_ - delays));  %#ok<PFBNS>
                        sig_samp = (chirp_waveform.sample(     t_));           %#ok<PFBNS>
                        sig = convn(imp_samp, sig_samp, 'full');
                    elseif(~sig_is_delta) % signal is not a delta: delay on signal
                        imp_samp = (tx_impulse_waveform.sample(t_));
                        sig_samp = (chirp_waveform.sample(     t_ - delays));
                        sig = convn(sig_samp, imp_samp, 'full');
                    else % both are deltas use linear resampling (sketchy)
                        warning('UltrasoundSystem:getKWaveSource:convolvingDeltas', ...
                            ['Using linear interpolation to convolve two delta functions with a shift. ', ...
                            'Use a small rect function to avoid this. ']);
                        sig_samp = (tx_impulse_waveform.sample(t_));
                        imp_samp = (chirp_waveform.sample(     t_));
                        t_samp = t_vec.Value - delays;
                        conv_samp = griddedInterpolant(t_, conv(sig_samp, imp_samp, 'full'), 'linear', 'none');
                        sig = reshape(conv_samp(vec(t_samp).'), size(t_samp));
                    end
                    
                    % normalize samples?
                    nsig = norm(vec(abs(sig)));
                    isig = nsig > 0;
                    sig(isig) = real(sig(isig)) .* inv(nsig(isig));
                    %%%

                    % add contributing signal from each subsample to each identical pixels
                    ind_umsk_unique = flip(unique(ind_umsk(ind_umsk ~= 0)));
                    for i=ind_umsk_unique(:)'
                        ksource_vec_vel(:,i,:) =  ksource_vec_vel(:,i,:) ...
                            + sum(vec_amp(:,ind_umsk == i) .* permute(sig(:,ind_umsk == i),[3 2 1]), 2);
                    end
                                        
                    %{
                    % initialize output variable
                    % ksource_vec_vel_sz = [3, numel(ind_msk_pulse), numel(t)];
                    kvv = zeros(ksource_vec_vel_sz);

                    % add contributing signal from each subsample to each identical pixels
                    ind_umsk_unique = flip(unique(ind_umsk));
                    for i=ind_umsk_unique, kvv(:,i,:) = sum(vec_amp(:,ind_umsk == i) .* permute(sig(:,ind_umsk == i),[3 2 1]), 2); end
                    
                    % add contributions from this element to all elements
                    ksource_vec_vel = ksource_vec_vel + kvv;
                    %}
                    
                    fprintf('.');
                end
                
                % move x/y/z vector to the corresponding fields and remove
                % vector version
                if any(ksource_vec_vel(1,:) > 0), ksource{puls}.ux = shiftdim(ksource_vec_vel(1,:,:),1); end
                if any(ksource_vec_vel(2,:) > 0), ksource{puls}.uy = shiftdim(ksource_vec_vel(2,:,:),1); end
                if any(ksource_vec_vel(3,:) > 0), ksource{puls}.uz = shiftdim(ksource_vec_vel(3,:,:),1); end

                % report
                fprintf('\nFinished processing pulse %i (%0.5g seconds).', puls, toc(tt))
            end            
            tocBytes(gcp('nocreate'));
        end
        
        function [resp, t_resp] = getKWaveReceive(self, kgrid, ksensor, sens_map, sensor_data, c0, varargin)
            % GETKWAVERECEIVE - Get the received response from the simulation
            %
            % [resp, t_resp] = GETKWAVERECEIVE(self, kgrid, ksensor, sens_map, sensor_data, c0)
            % takes the kWaveGrid kgrid, the kWaveSensor ksensor, the
            % transducer sensitivity map sens_map, a reference sound speed 
            % c0, and the received data sensor_data and creates the element
            % response resp and it's time axes t_resp.
            %
            % [...] = GETKWAVERECEIVE(..., Name, Value, ...) additionally
            % specifies options via name/value pairs.
            %
            % Inputs:
            %   device - GPU selection. 0 for no GPU, -1 to use the
            %            currently selected GPU, or n to select and reset 
            %            the GPU with id n in MATLAB.
            %   interp - interpolation method. Interpolation is provided by
            %   griddedInterpolant. This behaviour is likely to change.
            %
            % This signature (inputs/outputs) is likely to change.
            % 
            % See also TRANSDUCER/GETKWAVESENSOR GRIDDEDINTERPOLANT

            % set defaults
            kwargs.device = []; % use native type = -1 * logical(gpuDeviceCount); % use available gpu
            kwargs.interp = 'linear';

            % helper function
            vec = @(x) x(:);
            
            % set options
            for i = 1:numel(varargin)
                kwargs.(varargin{i}) = varargin{i+1};
            end

            % set device
            if isempty(kwargs.device), % do nothing
            elseif kwargs.device > gpuDeviceCount
                warning('Attempted to select device %i but only %i devices are available. Using the default device')
            else
                switch kwargs.device
                    case -1
                        % force on gpu
                        sensor_data = gpuArray(sensor_data);
                    case 0
                        % don't place on gpu
                    otherwise
                        % select gpu, then place
                        g = gpuDevice();
                        if g.Index ~= kwargs.device
                            sensor_data = gather(sensor_data);
                            g = gpuDevice(kwargs.device);
                            sensor_data = gpuArray(sensor_data);
                        end
                end
            end
            
            % define the rx impulse response function
            max_delay = hypot(hypot(kgrid.dx, kgrid.dy), kgrid.dz) / c0; % maximum grid delay
            rx_imp = self.rx.impulse(); % rx temporatl impulse response function
            n0  = floor((rx_imp.t0   - max_delay) / kgrid.dt);
            nend = ceil((rx_imp.tend + max_delay) / kgrid.dt);
            t_rx = vec(n0 : 1 : nend) * kgrid.dt; % includes zero if t0 <= 0 <= tend
            rx_imp_is_delta = (max(t_rx) - min(t_rx)) < max_delay; % check for delta functions - functions less than a max delay
            [nst, nfin] = deal(n0 + 0 + 1, nend + kgrid.Nt - 1 - 1); % time axis after convolution
            % [nst, nfin] = deal(min(n0, 0), max(nend, kgrid.Nt-1));
            % rx_pad = [-nst, nfin - (kgrid.Nt-1)];
            % T_data = kgrid.Nt; % size of kgrid data
            
            % set time to precision and device of sensor data
            t_rx = real(cast(t_rx, 'like', sensor_data));
            
            % get output signal sizing
            t_resp = vec(nst : 1 : nfin) * kgrid.dt; % includes zero if t0 <= 0 <= tend
            T_resp = numel(t_resp); % size of kgrid data after convolution
                        
            % get index mapping
            ind_msk_all = sort(find(ksensor.mask));
            
            % explicitly preallocate receiver outputs - 
            % implicitly preallocated when using a parfor loop
            % resp = zeros([T_resp, self.rx.numel], 'like', sensor_data);
            
            % splice
            [device, interp] = deal(kwargs.device, kwargs.interp);

            % for each (sub)element
            parfor (el = 1:self.rx.numel)
                
                % rename distributed data
                resp_samp = sensor_data;
                sens_map_el = sens_map(:,el);
                
                % get indices, delays, amplitudes for all
                % sub-elements
                [amp_i, ind_msk_i, delays_i] = arrayfun(@(sme)deal(...
                    vec(sme.amp) + zeros(size(sme.dist(:)')), ...
                    vec(sme.mask_indices).', ...
                    vec(sme.dist).' ./ c0 ...
                    ), sens_map_el, 'UniformOutput', false);
                
                % combine
                [amp, ind_msk, delays] = dealfun(@cell2mat, amp_i(:)', ind_msk_i(:)', delays_i(:)');
                
                % get sensor data indices for the subelement
                [~, ind_sen] = ismember(ind_msk, ind_msk_all);
                nSubEl = numel(ind_sen);

                % get sensor data (T_data x nSubEl)
                resp_samp = cat(1, ...
                    ... zeros([rx_pad(1), numel(ind_sen)]), ...
                    amp .* (resp_samp(ind_sen,:).') ...
                    ... zeros([rx_pad(2), numel(ind_sen)])...
                    ); % %#ok<PFBNS>
                
                % apply delay via convolution with a
                % shifted impulse response function
                % TODO: check that the signs here are okay
                if(~rx_imp_is_delta) % rx impulse not a delta: delay the impulse
                    imp_samp = rx_imp.sample(t_rx - delays); %#ok<PFBNS> % (T_data x nSubEl)
                    resp_samp = convd(resp_samp, imp_samp, 1, 'full', 'device', device); % (T_resp x nSubEl)
                else % impulse is a delta: use linear resampling (sketchy)
                    warning('UltrasoundSystem:kWaveSensor2ChannelData:convolvingDeltas', ...
                        ['Using linear interpolation to convolve two delta functions with a shift. ', ...
                        'Use a small rect function to avoid this. ']);
                    imp_samp = vec(rx_imp.sample(t_rx)); % (T_data x nSubEl)
                    conv_samp = convn(resp_samp, imp_samp, 'full'); % (T_resp x nSubEl) % GPU OOM?
                    conv_sampler = griddedInterpolant(t_resp, zeros([T_resp,1], 'like', t_resp), 'linear', 'none');
                    resp_samp = zeros([T_resp, nSubEl], 'like', conv_samp);
                    for i = nSubEl:-1:1
                        conv_sampler.Values(:) = conv_samp(:,i);
                        resp_samp(:,i) = conv_sampler(t_resp - delays(i));
                    end
                end
                
                % integrate over subelements
                resp(:, el) = trapz_(resp_samp,2);
            end
        end
    
        function [chd, cgrid] = kspaceFirstOrderND(self, target, element_subdivisions, varargin)
            % KSPACEFIRSTORDERND - Simulate channel data via k-Wave
            % 
            % chd = KSPACEFIRSTORDERND(self, target) simulates the Target
            % target and returns a ChannelData object chd via k-Wave.
            %
            % chd = KSPACEFIRSTORDERND(self, target, element_subdivisions)
            % uses the 1x2 array of element_subdivisions to subdivide the
            % elements prior to placing them on the simulation grid.
            %
            % [chd, cgrid] = kspaceFirstOrderND(...) also returns a
            % structure with the coordinate transforms from the kWaveGrid. 
            % This behaviour is likely to be deprecated
            %
            % chd = KSPACEFIRSTORDERND(self, target, element_subdivisions, 
            % Name, Value, ... ) specifies name value pairs.
            %
            % Inputs:
            %     PML_min - minimum (one-sided) PML size
            %     PML_max - maximum (one-sided) PML size
            %     CFL_max - maximum cfl number (for stability)
            %     buffer - x/y/z computational buffer (3 x 2)
            %     resolution_ratio - grid resolution as proportion of wavelength
            %     dims - dimensionality of the simulation (one of {2*,3})
            %     parcluster - parallel cluster for running simulations (use 0 for no cluster)
            %     bounds - minimum boundaries of the sim (3 x 2)
            %
            % Other Name/Value pairs that are valid for kWave's
            % kspaceFirstOrderND functions are valid here with the
            % exception of the PML definition
            %
            %
            % See also ULTRASOUNDSYSTEM/FULLWAVESIM


            % setup a default keyword arguments structure
            kwargs = struct( ...
                'PML_min', 20, ... minimum (one-sided) PML size
                'PML_max', 56, ... maximum (one-sided) PML size
                'CFL_max', 0.25, ... maximum cfl number (for stability)
                'buffer', [0, 0; 0, 0; 0, 0], ... % x/y/z computational buffer
                'resolution_ratio', 0.25, ... grid resolution as proportion of wavelength
                'dims', 2, ... number of dimensions of the simulation
                'parcluster', 0, ... parallel cluster for running simulations (use 0 for no cluster)
                'bounds', [0, 0; 0, 0; 0, 0], ... boundaries of the sim (minimum)
                'DataCast', 'gpuArray-single',...
                'DataRecast', false, ...
                'LogScale', true, ...
                'MovieName', 'kwave-sim', ...
                'PlotPML', false, ...
                'PlotSim', false, ...
                'RecordMovie', false, ...
                'Smooth', true ...
                );

            if nargin < 3 || isempty(element_subdivisions), 
                element_subdivisions = self.getLambdaSubDiv(kwargs.resolution_ratio, target.c0); 
            end


            % store NV pair arguments into kwargs
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            %% define the kgrid and medium sizing
            % check that the grid is larger than the transducer;
            % expand if this is not the case
            bc = [kwargs.bounds(:,1) >= sub(self.xdc.bounds,1,2), sub(self.xdc.bounds,2,2) >= kwargs.bounds(1,2)]; % boundary check
            if any(bc)
                warning('Expanding grid to encompass the transducer.');
                bd = self.xdc.bounds(); % transducer boundaries
                kwargs.bounds(find(bc)) = bd(bc); % set as the minumum boundaires for k-wave
            end

            % get arguments for the grid object
            kgrid_args = kwargs; % all input arguments
            gfields = {'PML_min', 'PML_max', 'CFL_max', 'CFL_max', 'resolution_ratio', 'dims', 'bounds'}; % grid args
            kgrid_args = rmfield(kgrid_args , setdiff(fieldnames(kgrid_args), intersect(fieldnames(kgrid_args), gfields))); % isolate grid args
            kgrid_args = struct2nvpair(kgrid_args);

            % get a kWaveGrid object ( also outputs the medium,
            % because it needs to set the CFL number )
            % TODO: use a Scan definition to define the region for the
            % simulation
            [kgrid, kgrid_PML_size, kgrid_origin, kgrid_size, kgrid_step, kmedium] ...
                = self.getkWaveGrid(target, kgrid_args{:});

            % get the translated / permuted grid sizing
            cgrid = struct(...
                'origin', kgrid_origin([2,3,1]) + [kgrid.y_vec(1); kgrid.z_vec(1); kgrid.x_vec(1);], ...
                'step'  , kgrid_step([2,3,1]), ...
                'size'  , kgrid_size([2,3,1]) ...
                );

            % get subdivision length
            el_sub_div = element_subdivisions;

            % define the sensor on the grid
            [ksensor_rx, ksensor_ind, sens_map] = self.rx.getKWaveSensor(kgrid, kgrid_origin, el_sub_div);

            % create the full sensor mask as the combination of all of the transducer
            % elements
            ksensor.mask = false;
            for el = 1:self.rx.numel, ksensor.mask = ksensor.mask | ksensor_rx{el}.mask; end

            % define the source signal(s) for each transmit in the
            % pulse sequence
            [ksource, t_sig, sig] = self.getKWaveSource(kgrid, kgrid_origin, el_sub_div, target.c0);

            % check the size of the temporal response and
            % maximum number of subelements. An array of this size
            % will be created
            nSubel = max(sum(arrayfun(@(s) numel(s.mask_indices), sens_map), 1), [], 2);
            warning(sprintf('Using up to %i subelements for %i points in time (%0.2fGB / %0.2fGB single/double (complex) temp variables)', nSubel, kgrid.Nt, prod(kgrid.Nt * nSubel) * 2 * [4, 8] / 2^30))


            %% Submit a k-wave simulation for each transmit
            tt_kwave = tic;

            % lam_min_grid_spac = min(kmedium.sound_speed,[],'all') / (self.us.xdc.fc + 0.5*self.us.xdc.bw);
            % lam_max_temp_freq = min(kmedium.sound_speed,[],'all') / max([kgrid.dx, kgrid.dy, kgrid.dz]) / 2;
            fprintf('There are %i points in the grid which is %0.3f megabytes per grid variable.\n', ...
                kgrid.total_grid_points, kgrid.total_grid_points*4/2^20);

            kwfields = { ...
                'DataCast', 'DataRecast', 'LogScale', 'MovieName', ... 
                'PlotPML', 'RecordMovie', 'Smooth', 'PlotSim' ...
                };

            % strip all other arguments from input
            kwave_args = kwargs; % all input arguments
            kwave_args = rmfield(kwave_args , setdiff(fieldnames(kwave_args), intersect(fieldnames(kwave_args), kwfields))); % isolate grid args
            
            % set global arguments
            kwave_args.PMLInside = false;
            kwave_args.PMLSize = kgrid_PML_size(1:kgrid.dim);

            
            % set per pulse properities
            kwave_args = repmat(kwave_args, [self.sequence.numPulse, 1]);
            for puls = 1:self.sequence.numPulse
                if(isfield(kwave_args(puls), 'MovieName'))
                    kwave_args(puls).MovieName = char(kwave_args(puls).MovieName + sprintf("_pulse%2.0i",puls));
                end
            end

            % select function
            switch kwargs.dims
                case 2, kspaceFirstOrderND_ = @kspaceFirstOrder2D;
                case 3, kspaceFirstOrderND_ = @kspaceFirstOrder3D;
            end

            % get all arguments
            kwave_args_ = arrayfun(...
                @(pulse) struct2nvpair(kwave_args(pulse)), ...
                1:self.sequence.numPulse, ...
                'UniformOutput', false ...
                );

            out = cell(self.sequence.numPulse, 2);
            [c0, Np] = deal(target.c0, self.sequence.numPulse);
            parfor (puls = 1:self.sequence.numPulse, kwargs.parcluster)
                fprintf('\nComputing pulse %i of %i\n', puls, Np);
                tt_pulse = tic;

                % prealloc
                outp = cell(1,2);
                % gpuDevice(1+mod(pulse-1, ngpu));

                % simulate
                sensor_data = kspaceFirstOrderND_(kgrid, kmedium, ksource{puls}, ksensor, kwave_args_{puls}{:}); %#ok<PFBNS>

                % Process the simulation data
                [outp{1}, outp{2}] = self.getKWaveReceive(kgrid, ksensor, sens_map, sensor_data, c0); %#ok<PFBNS>

                % enforce on the CPU - save to output var
                outp = cellfun(@gather, outp, 'UniformOutput', false);
                [out{puls,:}] = deal(outp{:});

                % report timing
                fprintf('\nFinished pulse %i of %i\n', puls, Np);
                toc(tt_pulse)
            end

            % get voltages (T x N x M)
            voltages = cast(cell2mat(shiftdim(out(:,1), -2)), 'single');

            % get timing (T x 1 x M) -> (T x 1)
            t_resp   = mode(cell2mat(shiftdim(out(:,2), -2)), 3);

            % TODO: make reports optional
            fprintf(string(self.sequence.type) + " k-Wave simulation completed in %0.3f seconds.\n", toc(tt_kwave));

            % get time axis
            time = t_resp + min(t_sig);

            % get numerically correct frequency
            fs_ = 1 ./ kgrid.dt;
            if fs_ > self.fs, %#ok<BDSCI> % fs_ is scalar
                fs_ = self.fs * round(fs_ / self.fs); % use the rounded upsampling factor
            else
                fs_ = self.fs / round(self.fs / fs_); % use the rounded downsampling factor
            end

            % Create a channel data object
            chd = ChannelData('t0', sub(time,1,1), 'data', voltages, 'fs', fs_);
        end
    end

    % k-Wave calls (new)
    methods
        function [chd, job, readfun] = kspaceFirstOrder(self, target, sscan, varargin)
            % KSPACEFIRSTORDERND - Simulate channel data via k-Wave
            % 
            % chd = KSPACEFIRSTORDERND(self, target) simulates the Target
            % target and returns a ChannelData object chd via k-Wave.
            %
            % chd = KSPACEFIRSTORDERND(self, target, sscan) operates using
            % the simulation region defined by the ScanCartesian sscan. The
            % step sizes in all dimensions must be identical for results to
            % be valid.
            %
            % [chd, job, readfun] = kspaceFirstOrderND(...) also returns
            % a parallel.Job job and a function readfun to read the output 
            % of the Job into a ChannelData object. If these outputs are
            % requested, chd is empty and the job is not submitted. The job
            % can be submitted with the submit function. When the job has
            % completed, it can be read into a ChannelData object with
            % readfun.
            %
            % [...] = KSPACEFIRSTORDERND(self, target, sscan, Name, Value, ...)
            % specifies name value pairs.
            %
            % Inputs:
            %     T - time limit of the simulation 
            %     PML - upper and lower bound on the PML size
            %     CFL_max - maximum cfl number (for stability)
            %     parcluster - parallel cluster for running simulations (use 0 for no cluster)
            %     ElemMapMethod - computational method for mapping elements to the grid. 
            %           Must be one of 
            %          {'nearest*, 'karray-direct', 'karray-depend'}. The
            %          'nearest' method uses the nearest pixel. The 
            %           'karray-direct' method uses the karray method but 
            %           avoids recomputing intermediate results. The 
            %           'karray-depend' method always uses the kWaveArray 
            %           methods, but can be slower. 
            %
            % Other Name/Value pairs that are valid for kWaveArray's 
            % constructor or kWave's kspaceFirstOrderND functions are 
            % valid here except for PML definition arguments.
            %
            % See also ULTRASOUNDSYSTEM/FULLWAVESIM


            % setup a default keyword arguments structure
            kwargs = struct( ...
                'T', [], ... simulation time (s)
                'PML', [20 56], ... (one-sided) PML size range
                'CFL_max', 0.25, ... maximum cfl number (for stability)
                'parcluster', 0, ... parallel cluster for running simulations (use 0 for no cluster)
                'ElemMapMethod', 'nearest', ... one of {'nearest'*,'linear','karray-direct', 'karray-depend'}
                'el_sub_div', self.getLambdaSubDiv(0.1, target.c0), ... element subdivisiosn (width x height)
                'UpsamplingRate', 10, ...
                'BLITolerance', 0.05, ...
                'BLIType', 'sinc', ... stencil - exact or sinc
                'DataCast', 'gpuArray-single',...
                'DataRecast', false, ...
                'LogScale', false, ...
                'MovieName', 'kwave-sim', ...
                'PlotPML', false, ...
                'PlotSim', false, ...
                'PlotScale', 'auto', ...
                ...'DisplayMask', 'off', ...
                'RecordMovie', false, ...
                'Smooth', true ...
                );

            % store NV pair arguments into kwargs
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % only supported with tx == rx for now
            assert(self.tx == self.rx, 'Transmitter and receiver must be identical.')

            % parse
            if isinf(sscan.dx), kwargs.el_sub_div(1) = 1; end % don't use sub-elements laterally for 1D sims
            if isinf(sscan.dy), kwargs.el_sub_div(2) = 1; end % don't use sub-elements in elevation for 2D sims

            % intialize empty outputs
            [chd, job, readfun] = deal([]);

            % get the kWaveGrid
            % TODO: check that the stpe sizes are all equal - this is
            % required by kWaveArray
            [kgrid, Npml] = getkWaveGrid(sscan, 'PML', kwargs.PML);

            % get the kWave medium struct
            kmedium = getMediumKWave(target, sscan);

            % make an iso-impedance medium
            kmedium_iso = kmedium;
            kmedium_iso.density = (target.c0 * target.rho0) ./ kmedium.sound_speed;

            % get the minimum necessary time step from the CFL
            dt_cfl_max = kwargs.CFL_max * min([sscan.dz, sscan.dx, sscan.dy]) / max(kmedium.sound_speed,[],'all');
            
            % use a time step that aligns with the sampling frequency
            dt_us = inv(self.fs);
            cfl_ratio = dt_cfl_max / dt_us;
        
            if cfl_ratio < 1
                warning('Upsampling the kWave time step for an acceptable CFL.');
                time_step_ratio = ceil(inv(cfl_ratio));            
            else
                % TODO: do (or don't) downsample? Could be set for temporal reasons
                time_step_ratio = 1;
            end
            kgrid.dt = dt_us / time_step_ratio;
            fs_ = self.fs * time_step_ratio;
            
            % kgrid to sscan coordinate translation mapping
            kgrid_origin = cellfun(@median, {sscan.z, sscan.x, sscan.y});

            % get the source signal
            txsig = conv(self.tx.impulse, self.sequence.pulse, 4*fs_); % transmit waveform, 4x intermediate convolution sampling
            apod =  self.sequence.apodization(self.tx); % transmit apodization (M x V)
            del  = -self.sequence.delays(self.tx); % transmit delays (M x V)
            txn0 = floor((txsig.t0   + min(del(:))) * fs_); % minimum time sample - must pass through 0
            txne = ceil ((txsig.tend + max(del(:))) * fs_); % maximum time sample - must pass through 0
            t_tx = shiftdim((txn0 : txne)' / fs_, -2); % transmit signal time indices (1x1xT')
            txsamp = permute(apod .* txsig.sample(t_tx - del), [3,1,2]); % transmit waveform (T' x M x V)
            
            % define the sensor on the grid 
            % generate {psig, mask, elem_weights}
            switch kwargs.ElemMapMethod
                case 'nearest'
                    pg = sscan.getImagingGrid(); % grid points {x, y, z}
                    pg = cell2mat(cellfun(@(x)shiftdim(x,-1), pg(:), 'UniformOutput',false)); % -> 3 x Z x X x Y
                    pn = self.rx.positions; % element positions
                    [Nx, Ny, Nz] = dealfun(@(n) n + (n==0), kgrid.Nx, kgrid.Ny, kgrid.Nz); % kwave sizing
                    mask = false(Nx, Ny, Nz); % grid size
                    assert(all(size(mask,1:3) == sscan.size), 'kWave mask and sacn size do not correspond.');
                    for n = self.rx.numel:-1:1, % get nearest pixel for each element
                        ind(n) = argmin(vecnorm(pn(:,n) - pg,2,1),[],'all', 'linear');
                    end
                    mask(ind) = true;
                    psig = permute(txsamp,[2,1,3]); % -> (J' x T' x V) with M == J'
                    elem_weights = eye(self.rx.numel);

                case 'linear'
                    % get ambient sound speed
                    c0map = target.c0;

                    % get a mapping of delays and weights to all
                    % (sub-)elements (J' x M)
                    [mask, el_weight, el_dist, el_ind] = self.rx.elem2grid(sscan, kwargs.el_sub_div);% perm(X x Y x Z), (J' x M)
                    el_map_grd = ((1:nnz(mask))' == el_ind(:)'); % matrix mapping (J' x J'')

                    % apply to transmit signal: for each element
                    [del, apod, t_tx] = dealfun(@(x)shiftdim(x, -1), del, apod, t_tx); % (1 x M x V), (1 x 1 x 1 x T')
                    tau = t_tx - del - el_dist/c0map; % J''' x M x V x T'
                    psig = permute(el_weight .* apod .* txsig.sample(tau), [1,2,4,3]); % per sub-element transmit waveform (J''' x M x T' x V)
                    psig = reshape(psig, [prod(size(psig,1:2)) size(psig,3:4)]); % per element transmit waveform (J'' x T' x V)
                    psig = pagemtimes(double(el_map_grd), psig); % per grid-point transmit waveform (J' x T' x V)


                case {'karray-direct', 'karray-depend'}
                    % [ksensor_rx, ksensor_ind, sens_map] = self.rx.getKWaveSensor(kgrid, kgrid_origin, el_sub_div);
                    karray_opts = {'BLITolerance', kwargs.BLITolerance, 'UpsamplingRate', kwargs.UpsamplingRate, 'BLIType', kwargs.BLIType};
                    karray = kWaveArray(self.rx, kgrid.dim, kgrid_origin, karray_opts{:});
                    mask = karray.getArrayBinaryMask(kgrid);

                    % assign source for each transmission (J' x T' x V)
                    % TODO: abstract this logic to come from transducer directly so
                    % it can be used in fullwave or some other FDTD method
                    switch kwargs.ElemMapMethod
                        case 'karray-direct'
                            % get the offgrid source weights
                            elem_weights = arrayfun(@(i){sparse(vec(karray.getElementGridWeights(kgrid, i)))}, 1:self.xdc.numel);  % (J x {M})
                            elem_weights = cat(2, elem_weights{:}); % (J x M)
                            elem_weights = elem_weights(mask(:),:); % (J' x M)
                            psig = pagemtimes(full(elem_weights), 'none', txsamp, 'transpose'); % (J' x M) x (T' x M x V) -> (J' x T' x V)

                            % get the offgrid source sizes
                            elem_meas = arrayfun(@(i)karray.elements{i}.measure, 1:self.xdc.numel);
                            elem_dim  = arrayfun(@(i)karray.elements{i}.dim    , 1:self.xdc.numel);
                            elem_norm = elem_meas ./ (kgrid.dx) .^ elem_dim; % normalization factor
                            elem_weights = elem_weights ./ elem_norm; 

                        case 'karray-depend', % compute one at a time and apply casting rules
                            psig = cellfun(@(x) ...
                                {cast(karray.getDistributedSourceSignal(kgrid, x.'), 'like', x)}, ...
                                num2cell(real(txsamp), [1,2]) ...
                                );
                            psig = cat(3, psig); % (J' x T' x V)
                    end
            end

            % define the source and sensor
            ksensor.mask = mask; % pixels to record
            ksensor.record = {'u'}; % record the original particle velocity
            % ksensor.record = {'u','u_non_staggered', 'p'}; % record everything
            for v = self.sequence.numPulse:-1:1, ksource(v).ux = real(sub(psig,v,3)); end % set transmit pulses
            [ksource.u_mask] = deal(mask); % set transmit aperture mask

            % set the total simulation time: default to a single round trip at ambient speed
            if isempty(kwargs.T)
                kwargs.T = 2 * (vecnorm(range([sscan.xb; sscan.yb; sscan.zb], 2),2,1) ./ target.c0);
            end
            Nt = 1 + floor((kwargs.T / kgrid.dt) + max(range(t_tx,1))); % number of steps in time
            kgrid.setTime(Nt, kgrid.dt);

            % get the receive impulse response function
            rx_imp = self.rx.impulse;
            t_rx = rx_imp.getSampleTimes(fs_);
            rx_sig = gather(real(rx_imp.sample(t_rx(:))));

            % simulation start time
            t0 = gather(t_tx(1) + t_rx(1));

            %% Submit a k-wave simulation for each transmit
            tt_kwave = tic;

            fprintf('There are %i points in the grid which is %0.3f megabytes per grid variable.\n', ...
                kgrid.total_grid_points, kgrid.total_grid_points*4/2^20);

            % strip all other arguments from input
            nonkwfields = {'T', 'PML','CFL_max', 'parcluster', 'ElemMapMethod', 'el_sub_div', ... % not kWave args
                 'BLIType', 'BLITolerance','UpsamplingRate',  ... % not kspaceFirstOrder args
                }; 
            kwave_args = rmfield(kwargs, nonkwfields); % forward all others
            
            % set global arguments: these are always overridden
            kwave_args.PMLInside = false;
            kwave_args.PMLSize = Npml(1:kgrid.dim);
            
            % make a new movie name for each pulse
            kwave_args = repmat(kwave_args, [self.sequence.numPulse, 1]);
            mv_nm = {kwave_args.MovieName} + sprintf("_pulse%2.0i",1:self.sequence.numPulse);
            [kwave_args.MovieName] = deal(mv_nm{:});

            % select function
            switch kgrid.dim
                case 1, kspaceFirstOrderND_ = @kspaceFirstOrder1D;
                case 2, kspaceFirstOrderND_ = @kspaceFirstOrder2D;
                case 3, kspaceFirstOrderND_ = @kspaceFirstOrder3D;
                otherwise, error("Unsupported dimension size (" + kgrid.dim  +").");
            end

            % get all arguments
            kwave_args_ = arrayfun(...
                @(pulse) struct2nvpair(kwave_args(pulse)), ...
                1:self.sequence.numPulse, ... 
                'UniformOutput', false ...
                );

            % processing step: get sensor data, enforce on CPU (T x N)
            switch kwargs.ElemMapMethod
                case 'karray-depend'
                    proc_fun = @(x) gather(convn( karray.combineSensorData(kgrid, x.ux).', rx_sig, 'full'));
                case {'karray-direct'}
                    proc_fun = @(x) gather(convn(x.ux.' * full(elem_weights)), rx_sig, 'full'); % (J' x T)' x (J' x N) -> T x N
                case {'nearest'}
                    proc_fun = @(x) gather(x.ux.'); % -> (T x N)
                case 'linear'
                    % create the advanced impulse response function with
                    % which to convolve the output
                    vec = @(x)x(:);
                    N = self.rx.numel;
                    rx_sig = gather(real(rx_imp.sample(t_rx(:)' + el_dist(:)/c0map))); % J'' x T''
                    el_map_el = ((1:N) == vec(ones(size(el_ind)) .* (1:N)))'; % map from convolved samples to elements
                    proc_fun = @(x) gather(...
                         (el_map_el * (el_weight(:) .* convd(el_map_grd' * x.ux, rx_sig, 2, 'full'))).' ... % [(N x J'') x [[(J'' x J') x (J' x T')] x (T' x T | J'')]]' -> (T x N) 
                         ...
                        ); % [(N x J'') x (J'' x T)]' -> T x N

                otherwise, warning('Unrecognized mapping option - mapping to grid pixels by default.');
                    proc_fun = @(x) gather(x.ux.');
            end

            % TODO: provide more options to toggle between cluster job and
            % parfor behaviour - for now, only parcluster == 0 -> parfor
            switch kwargs.parcluster
                case 0 % process directly in a parfor loop if 0
                     out = cell(self.sequence.numPulse, 1); % init
                    [Np] = deal(self.sequence.numPulse); % splice

                    % for puls = self.sequence.numPulse:-1:1
                    parfor (puls = 1:self.sequence.numPulse, kwargs.parcluster)
                        % TODO: make this part of some 'info' logger or something
                        fprintf('\nComputing pulse %i of %i\n', puls, Np);
                        tt_pulse = tic;

                        % simulate
                        sensor_data     = kspaceFirstOrderND_(kgrid, kmedium    , ksource(puls), ksensor, kwave_args_{puls}{:}); %#ok<PFBNS>
                        sensor_data_iso = kspaceFirstOrderND_(kgrid, kmedium_iso, ksource(puls), ksensor, kwave_args_{puls}{:}); 

                        % Process the simulation data
                        out{puls} = proc_fun(sensor_data) - proc_fun(sensor_data_iso); %#ok<PFBNS> data is small 

                        % report timing % TODO: make this part of some 'info' logger or something
                        fprintf('\nFinished pulse %i of %i\n', puls, Np);
                        toc(tt_pulse)
                    end

                    % create ChannelData objects
                    chd = ChannelData('data', cat(3, out{:}), 't0', t0, 'fs', fs_);

                    % TODO: make reports optional
                    fprintf(string(self.sequence.type) + " k-Wave simulation completed in %0.3f seconds.\n", toc(tt_kwave));

                otherwise
                    % make a job on the cluster
                    clu = kwargs.parcluster;
                    job = createJob(clu, 'AutoAddClientPath', true, 'AutoAttachFiles',true);

                    % arguments for each simulation
                    parfor (puls = 1:self.sequence.numPulse, 0)
                        kargs_sim{puls} = [...
                            {kgrid, kmedium    , ksource(puls), ksensor}, ...
                            kwave_args_{puls}(:)'...
                            ];
                        kargs_sim_iso{puls} = [...
                            {kgrid, kmedium_iso, ksource(puls), ksensor}, ...
                            kwave_args_{puls}(:)'...
                            ];
                    end
                    
                    % add simulation and processing task
                    job.createTask(@(varargin) proc_fun(kspaceFirstOrderND_(varargin{:})), 1, [kargs_sim_iso, kargs_sim]);

                    % function to read into a ChannelData object
                    % reshape and convole with receive impulse-response
                    % we can make this a lambda because we only
                    % have one output per task
                    V = self.sequence.numPulse;
                    readfun = @(job) ... 
                        ChannelData('t0', t0, 'fs', fs_, 'data', ... 
                        diff(cell2mat(arrayfun(@(t)t.OutputArguments, reshape(job.Tasks, [1,1,V,2]))),1,4) ... 
                        );

                    % if no job output was requested, run the job and
                    % create the ChannelData
                    if nargout < 2, 
                        submit(job); 
                        wait(job);
                        chd = readfun(job);

                        % TODO: make reports optional
                        fprintf(string(self.sequence.type) + " k-Wave simulation completed in %0.3f seconds.\n", toc(tt_kwave));
                    end
            end
        end
    end
    
    % Beamforming
    methods(Access=public)
        function b = DAS(self, chd, c0, varargin)
            % DAS - Delay and sum
            %
            % b = DAS(us, chd, c0) performs delay-and-sum beamforming on 
            % the ChannelData chd using an assumed sound speed of c0. The
            % ChannelData must conform to the delays given by the Sequence 
            % model in us.sequence. The output is defined on the Scan 
            % object us.scan.
            % 
            % b = DAS(us, chd, c0) specifies the receive aperture
            % reduction function (defaults to summation)
            %
            % 
            % Inputs:
            %  chd      - A ChannelData object  
            %  c0       - A sound speed, or any object with a .c0 property
            %
            % Name/Value pair arguments
            %  prec -  compute precision of the positions 
            %            {'single'* | 'double'}
            %  device - GPU device index: -1 for default gpu, 0 for cpu, n
            %           to select (and reset!) a gpuDevice.
            %  interp - select a method of interpolation for the transmit
            %           synthesis
            %  apod   - define the apodization. The apodization is defined
            %           as a 5D compatible array across 
            %           I1 x I2 x I3 x N x M where N is the number of 
            %           receivers, M is the number of transmits, and I[1-3]
            %           are the dimensions according to us.scan. Typically,
            %           apod should be singleton in at least two or three
            %           dimensions.
            %
            %  keep_rx - whether to keep the receive dimension rather than
            %           sum. The default is false.
            %
            %  keep_tx - whether to keep the transmit dimension rather than
            %           sum. The default is false. If keep_tx is true, then
            %           keep_rx must also be true.
            %
            %
            % outputs:
            %   - X\Y\Z:    3D coordinates of the output
            %   - B:        B-mode image
            
            % TODO: switch to using a struct kwargs to facilitate later
            % using the arguments blocks for R2020b+
            % default name-value pair arguments
            prec = 'single';
            device = int64(-1 * logical(gpuDeviceCount)); % {0 | -1} -> {CPU | GPU}
            apod = 1;
            interp_args = {};
            sumtx = true;
            sumrx = true; 
            
            % name-value pairs
            for n = int64(1):2:numel(varargin)
                switch varargin{n}
                    case 'prec'
                        prec = varargin{n+1};
                    case 'device'
                        device = varargin{n+1};
                    case 'interp'
                        interp_args = varargin(n:n+1);
                    case 'apod', apod = varargin{n+1}; 
                    case 'keep_tx', sumtx = ~varargin{n+1};
                    case 'keep_rx', sumrx = ~varargin{n+1};
                    otherwise
                        error('Unrecognized name-value pair');
                end
            end

            % if we want to keep tx, we must keep rx too because of
            % available DAS functions
            if ~sumtx && sumrx
                error('Unable to keep transmit dimension but not receive dimension. Try bfDAS.'); 
            end

            % accept c0 as a value, or an object with a c0 property
            if ~isnumeric(c0) && ~istall(c0), c0 = c0.c0; end % accept tall numbers too for now.

            % make sure t0 is a scalar
            if ~isscalar(chd.t0), chd = rectifyt0(chd); end
            
            % get positions of the imaging plane 
            [X, Y, Z, image_size] = self.scan.getImagingGrid();

            % reshape into I x N x M
            apod_args = {'apod', apod};
            
            % convert to x/y/z in 1st dimension
            P_im = permute(cat(4, X, Y, Z),[4,1,2,3]); % 3 x I1 x I2 x I3 = 3 x I
            
            % get positions of the aperture(s)
            P_tx = self.tx.positions(); % cast(self.tx.positions(), 'like', time(end)); % 3 x M
            P_rx = self.rx.positions(); % cast(self.rx.positions(), 'like', time(end)); % 3 x N
            
            % get the beamformer arguments
            dat_args = {chd.data, chd.t0, chd.fs, c0, 'device', device, 'position-precision', prec}; % data args
            ext_args = [interp_args, apod_args]; % extra args
            
            switch self.sequence.type
                case 'FSA'
                    pos_args = {P_im, P_rx, P_tx, [0;0;1]};
                case 'PW'
                    pos_args = {P_im, P_rx, [0;0;0], self.sequence.focus}; % TODO: use origin property in tx sequence
                    ext_args{end+1} = 'plane-waves'; 
                case 'VS'
                    pos_args = {P_im, P_rx, self.sequence.focus, [0;0;1]};
            end

            % beamform and collapse the aperture
            if      sumtx &&  sumrx, fun = 'DAS'; 
            elseif  sumtx && ~sumrx, fun = 'SYN'; 
            elseif ~sumtx && ~sumrx, fun = 'BF';
            end

            % beamform the data (I1 x I2 x I3 x N x M x F x ...)
            b = beamform(fun, pos_args{:}, dat_args{:}, ext_args{:});

            % move data dimension, back down raise aperture dimensions (I1 x I2 x I3 x F x ... x N x M)
            b = permute(b, [1:3,6:ndims(b),4:5]);
        end
        
        function chd = focusTx(self, chd0, seq, varargin)
            % FOCUSTX - Synthesize transmits
            %
            % chd = FOCUSTX(self, chd0) focuses the FSA ChannelData chd0 by
            % linearly synthesizing transmits (e.g. delay and sum across 
            % transmits)
            %
            % chd = FOCUSTX(self, chd0, seq) uses the Sequence seq to focus
            % the data instead of self.sequence.
            %
            % chd = FOCUSTX(..., Name, Value, ...) uses name-value pairs to
            % specify 
            %
            % Name/value pair arguments
            %  length -               Length of DFT / time vector
            %  interp -               Interpolation method
            %
            %  Interpolation methods are passed to the ChannelData/sample
            %  function. An additional method 'freq' is provided to delay
            %  the data in the frequency domain.
            %
            % Outputs:
            %   chd -                   New ChannelData object
            %
            % See also CHANNELDATA/SAMPLE
            
            % defaults
            kwargs.length = [];
            kwargs.interp = 'cubic';
            
            % focus using self's sequence 
            if nargin < 3 || isempty(seq), seq = self.sequence; end

            % parse optional inputs
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end
            L = kwargs.length;

            % Copy semantics
            chd = copy(chd0);

            % nothing to do for FSA acquisitions
            switch seq.type, case 'FSA', return; end 

            % dist/time to receiver
            tau_focal = - seq.delays(self.tx); % M x M'
            apod      =   seq.apodization(self.tx); % [1|M] x [1|M']

            % resample only within the window where we currently have data.
            nmin = floor(min(tau_focal,[],'all') .* chd.fs); % minimum sample time
            nmax =  ceil(max(tau_focal,[],'all') .* chd.fs); % maximum sample time
            chd.t0    =    chd.t0 + nmin / chd.fs; % shift time axes forwards to meet minimum sample time
            tau_focal = tau_focal - nmin / chd.fs; % shift delays backwards to meet time axes
            chd = zeropad(chd,0,(nmax - nmin)); % expand time axes to capture all data
            
            % legacy: pick new signal length
            if kwargs.interp == "freq",
                if isempty(L) % default
                    L = chd.T;
                elseif ischar(L) % strategy
                    switch L
                        case 'min'
                            L = chd.T;
                        case 'pow2'
                            L = 2^(nextpow2(chd.T));
                        otherwise
                            error('Unrecognized strategy for DFT length');

                    end
                elseif isscalar(L) % size
                    if L < chd.T, warning('Signal length may be too short!'); end
                else
                    error("L must be a scalar or one of {'min' | 'pow2'}");
                end
                chd = zeropad(chd,0,max(0, L - chd.T)); % expand only
            end
            %}
            
            % align dimensions
            D = 1+ndims(chd.data); % get a free dimension for M'
            tau_focal = swapdim(tau_focal, [1,2], [chd.mdim, D]); % move data
            apod      = swapdim(apod     , [1,2], [chd.mdim, D]); % move data

            % choose an interpolation method: frequency is implemented here
            %{
            switch kwargs.interp
                case 'freq'
                    % frequency vector
                    l = shiftdim(0 : 1 : (L - 1), 1); % L x 1 x 1
                    x = chd.data; % extract the data

                    % apply phase shifts and sum
                    wL = apod .* exp(-2i * pi * chd.fs  .* tau_focal ./ L); % steering vector (1 x 1 x M x ... x M')
                    x = fft(x, L, chd.tdim); % put data in freq domain (L x N x M x 1 x F x ...)
                    y = tenmul(wL.^l, x, chd.mdim, chd.tdim); % apply phase shift and sum over transmits (L x N x 1 x M' x F x ...)
                    z = ifft(y, L, chd.tdim, 'nonsymmetric'); % back to time (T' x N x 1 x M' x F x ...)
                    z = circshift(z, nwrap, chd.tdim); % handle FFT wrap-around
                    
                    % swap M' and M dimensions to return in same format
                    z = swapdim(z, chd.mdim, D);

                otherwise
                    [MP, MPa] = deal(size(tau_focal,D), size(apod,D)); % sizing: number of transmits
                    for m = MP:-1:1 % per transmit (implicit pre-allocation)
                        % sample, apodize, and sum over transmits
                        % (T' x N x {M'} x F x ...)
                        a = sub(apod,min(m,MPa),4); % apodization for this transmit
                        z{m} = chd.sample(chd.time - sub(tau_focal,m,D), kwargs.interp, a, chd.mdim); 
                    end
                    z = cat(chd.mdim, z{:}); % combine (T' x N x M' x F x ...)
            end
            %}
            
            % sample and store
            z = chd.sample(chd.time - tau_focal, kwargs.interp, apod, chd.mdim); % sample (perm(T' x N x 1) x F x ... x M')
            z = swapdim(z, chd.mdim, D); % replace transmit dimension (perm(T' x N x M') x F x ...)
            chd.data = z;% store output channel data % (perm(T' x N x M') x F x ...)
        end
        
        function [B, X, Y, Z] = DASEikonal(self, chd, medium, cgrid, rcvfun, varargin)
            % DASEIKONAL - Delay and sum FSA data with an assumed sound speed
            %
            % Inputs:
            %   chd -                   ChannelData object of FSA data
            %   medium -                Medium object
            %   rcvfun -                Receive aperture accumulation function
            %                           (defaults to summation)
            %   cgrid -                 A structure defining the grid on
            %                           which to calculate the speed of
            %                           sound using the Medium object.
            %     cgrid.origin (3 x 1)  The origin of the grid in x/y/z
            %                           defined as the point with the
            %                           smallest value in each dimension
            %     cgrid.step   (3 x 1)  The distance between grid points in
            %                           x/y/z. Inf/NaN for a sliced
            %                           dimension.
            %     cgrid.size   (3 x 1)  The number of grid points in each
            %                           dimension. 1 for a sliced dimension
            %                           
            %   Name/Value pair arguments:
            %   prec -                 compute precision of the positions
            %                           {'single'* | 'double'}
            %   device -               GPU device index: -1 for default
            %                           gpu, 0 for cpu, n for device index
            %
            % Outputs:
            %   B -        B-mode image
            %   X\Y\Z -    3D coordinates of the output
            %
            % See also DAS
            
            
            warning('This function is currently unsupported');
            % default name-value pair arguments
            prec = 'single';
            device = int64(-1 * logical(gpuDeviceCount)); % {0 | -1} -> {CPU | GPU}

            % name-value pairs
            for n = int64(1):2:numel(varargin)
                switch varargin{n}
                    case 'prec'
                        prec = varargin{n+1};
                    case 'device'
                        device = varargin{n+1};
                    otherwise
                        error('Unrecognized name-value pair');
                end
            end
            
            % get positions of the imaging plane
            [X, Y, Z, image_size] = self.scan.getImagingGrid();
            
            % convert to x/y/z in 1st dimension
            shp = @(n) n(:)';
            P_im = cat(1, shp(X), shp(Y), shp(Z)); % 3 x I

            % get sample points of the speed of sound grid
            for d = 3:-1:1
                if ~isfinite(cgrid.step(d)), dp = 0; else, dp = cgrid.step(d); end
                P_c_vec{d} = cgrid.origin(d) + dp*colon(1,cgrid.size(d)); 
            end
            shf = @(n)shiftdim(n,-1);
            [Xc, Yc, Zc] = ndgrid(P_c_vec{:});
            P_cm = cat(1, shf(Xc), shf(Yc), shf(Zc)); % 3 x Ic
            
            % get speed of sound map on the c-grid
            [c, ~] = medium.getPropertyMap(P_cm); % X x Y x Z
            
            % get positions of the aperture(s)
            P_tx = self.tx.positions(); % 3 x M
            P_rx = self.rx.positions(); % 3 x N
            
            % get the beamformer arguments
            dat_args = {cgrid, chd.data, c, chd.t0, chd.fs, 'device', 0, 'position-precision', prec}; %#ok<PROPLC> % data args
            ext_args = {}; % extra args
            
            switch self.sequence.type
                case 'FSA'
                    pos_args = {P_im, P_rx, P_tx};
                case 'PW'
                    pos_args = {P_im, P_rx, self.sequence.focus}; 
                    ext_args{end+1} = 'plane-waves'; 
                    error('Eikonal beamforming not implemented for plane wave transmits.');
                case 'VS'
                    pos_args = {P_im, P_rx, self.sequence.focus};
                    error('Eikonal beamforming not implemented for virtual source transmits.');
            end
            
            % beamform we the speed of sound data and collapse the aperture
            if nargin < 6 || isempty(rcvfun)
                % beamform the data (I x 1)
                B = cbeamform('DAS', pos_args{:}, dat_args{:}, ext_args{:}); 
            else
                % beamform the data (I x N)
                B = cbeamform('SYN', pos_args{:}, dat_args{:}, ext_args{:});
                
                % apply recieve aperture function in dim 2 and shape up
                B = rcvfun(B, 2); % I x 1
            end
            
            % make the shape consistent
            B = reshape(B, image_size);        
        end
    
        function b = bfAdjoint(self, chd, c0, varargin)
            % BFADJOINT - Adjoint method beamformer
            %
            % b = BFADJOINT(self, chd) beamforms the ChannelData chd using
            % an adjoint matrix method.
            % 
            % b = BFADJOINT(self, chd, c0) uses an assumed sound speed c0
            % for the green's functions. The default is 1540 (m/s).
            %             
            % See also ULTRASOUNDSYSTEM/DAS ULTRASOUNDSYSTEM/FOCUSTX
            % ULTRASOUNDSYSTEM/BFEIKONAL

            % TODO: include apodization, device, other keyword arguments

            warning('This function is currently unsupported.'); 
            
            % options
            kwargs.fthresh = -40; % threshold for including frequencies
            kwargs.apod = 1; % apodization
            kwargs.Nfft = chd.T; % FFT-length

            % parse inputs
            if nargin < 3 || isempty(c0), c0 = 1540; end

            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % move the data to the frequency domain, being careful to
            % preserve the time axis
            K = kwargs.Nfft; % DFT length
            f = shiftdim(chd.fs * (0 : K - 1)' / K, 1-chd.tdim); % frequency axis
            df = chd.fs * 1 / K; % frequency step size
            x = fft(chd.data,K,chd.tdim); % K x N x M x ...
            x = x .* exp(-2i*pi*f.*chd.t0); % phase shift to re-align time axis

            % choose frequencies to evaluate
            xmax = max(x, [], chd.tdim); % maximum value per trace
            f_val = mod2db(x) - mod2db(xmax) >= kwargs.fthresh; % threshold 
            f_val = f_val & f < chd.fs / 2; % positive frequencies only
            f_val = any(f_val, setdiff(1:ndims(x), chd.tdim)); % evaluate only freqs across aperture/frames that is above threshold
            
            % get the pixel positions
            D = max(3, gather(ndims(chd.data))); % >= 4
            Pi = self.scan.getImagingGrid();
            Pi = cellfun(@(x) {shiftdim(x, -D)}, Pi); % place I past max data dims 5-7
            Pi = cat(1,Pi{:}); % 3 x 1 x 1 x 1 x [I]

            % get the receive apodization, spliced if it can be applied
            apod = kwargs.apod;
            [a_n, a_m, a_mn] = deal(1);
            if all(size(apod, 1:3) == 1) % image is scalar, apply to data
                a_mn = shiftdim(apod, 3); % N x V
            elseif size(apod, 4) == 1 % receive is scalar, apply over tx
                ord = [D+(1:4), 2]; % send dims 1-5 here
                ord = [ord, setdiff(1:max(ord, 5), ord)]; % complete set of dimensions
                a_m = ipermute(apod, ord); % 1 x V x 1 x 1 x [I] x 1
            elseif size(apod, 5) == 1 % transmit is scalar, apply over rx
                ord = [D+(1:3), 2]; % send dims 1-4 here
                ord = [ord, setdiff(1:max(ord, 5), ord)]; % complete set of dimensions
                a_n = ipermute(apod, ord); % 1 x N x 1 x 1 x [I]
            else % none of the following is true: this request is excluded for now
                error('Unable to apply apodization due to size constraints. Apodization must be scalar in the transmit dimension, receive dimension, or all image dimensions.')
            end
            
            % get the delays for the transmit/receive green's matrix
            % kernels
            tau_tx = vecnorm(self.tx.positions() - Pi,2,1) / c0; % 1 x M x 1 x 1 x [I]
            tau_rx = vecnorm(self.rx.positions() - Pi,2,1) / c0; % 1 x N x 1 x 1 x [I]

            % get the transmit steering vector weights and delays
            del_tx  = self.sequence.delays(self.tx); % M x V
            apod_tx = self.sequence.apodization(self.tx); % M x V

            % transform to frequency step kernels
            w_rx    = exp(-2i*pi*df.*tau_rx); %  receive greens function
            w_tx    = exp(-2i*pi*df.*tau_tx); % transmit greens function
            w_steer = exp(-2i*pi*df.*del_tx); % transmit steering delays

            % TODO: cast data type for efficency?
            [w_tx, w_rx, w_steer, apod_tx] = dealfun(@(w) cast(w, 'like', real(x)), w_tx, w_rx, w_steer, apod_tx);
            b = repmat(cast(0, 'like', x), [1, size(Pi,2:ndims(Pi))]);

            % beamform one frequency at a time to save memory space
            h = waitbar(0,'Beamforming ...');
            % parfor (k = 1:K)
            for k = gather(find(f_val)') % skip unimportant frequencies

                % report progress
                if isvalid(h), waitbar(k/K/2, h, char("Beamforming: " + (gather(f(k))/1e6) + " MHz")); end

                % select frequency
                % xk = shiftdim(sub(x,k,1),1); % data, in freq. domain (N x V x ...)
                xk = permute(sub(x,k,chd.tdim), [chd.ndim, chd.mdim, chd.tdim, 4:D]);
                
                % TODO: adapt for data out of order
                % TODO: adapt for tall types

                % compute the greens functions on transmit/receive
                G_tx = w_tx.^(k-1); % 1 x M x 1 x 1 x [I]
                G_rx = w_rx.^(k-1); % 1 x N x 1 x 1 x [I]

                % compute the inverse steering vector on transmit
                T_tx = apod_tx .* w_steer.^(k-1); % M x V
                A_tx = pagemtimes(G_tx, T_tx); % 1 x V x 1 x 1 x [I]
                Ainv_tx = pagetranspose(A_tx); % V x 1 x 1 x 1 x [I] % make a column vector
                Ainv_tx = Ainv_tx ./ vecnorm(Ainv_tx, 2, 1); % normalize the power
                
                % apodize, delay, and sum the data for this frequency
                % only 1 of the a_* will contain the apodization
                yn = a_m .* pagemtimes(a_n .* conj(G_rx), a_mn .* xk); % 1 x V x 1 x 1 x [I]
                y  = pagemtimes(yn, conj(Ainv_tx)); % 1 x 1 x 1 x 1 x [I]
                
                % integrate over all frequencies
                b = b + y; % 1 x 1 x 1 x 1 x [I]

            end           
            if isvalid(h), close(h); end

            % move to image dimensions
            % b = shiftdim(b,D); % [I] x ...
            b = swapdim(b, 1:3, D+(1:3));
        end    

        function b = bfEikonal(self, chd, medium, cscan, varargin)
            % BFEIKONAL - Delay-and-sum beamformer with Eikonal delays
            %
            % b = BFEIKONAL(self, chd, medium, cscan) creates a b-mode 
            % image b from the ChannelData chd and Medium medium using the 
            % delays given by the solution to the eikonal equation defined 
            % on the ScanCartesian cscan. The transmitter and receiver must
            % fall within the cscan. The step size in each dimension must
            % be identical. The equation is solved via the fast marching 
            % method. 
            % 
            % b = BFEIKONAL(..., Name,Value, ...) defines additional
            % parameters via Name/Value pairs
            %
            % Inputs:
            %   
            % keep_rx -     setting this to true returns an extra dimension
            %               containing the data for each receive prior to
            %               summation. The default is false.
            %
            % keep_tx -     setting this to true returns an extra dimension
            %               containing the data for each transmit prior to
            %               summation. The default is false.
            %
            % interp  -     specifies the method for interpolation. Support 
            %               is provided by the ChannelData/sample method. 
            %               The default is 'linear'.
            % 
            %   
            % apod    -     specifies and ND-array A for apodization. A must 
            %               be broadcastable to size  I1 x I2 x I3 x N x M 
            %               where I1 x I2 x I3 is the size of the image, N 
            %               is the number of receivers, and M is the number
            %               of transmits. The default is 1.
            % 
            %   
            % bsize   -     sets the ChannelData block size to at most B 
            %               transmits at a time in order to limit memory 
            %               usage. If memory is a concern, lowering B may 
            %               help. If there is ample memory available, a 
            %               higher value of B may yield better performance.
            % 
            %   
            % parcluster  - specifies a parcluster object for 
            %               parallelization. The default is the current 
            %               parallel pool returned by gcp.
            %               Setting clu = 0 avoids using a parallel cluster
            %               or pool. 
            %               A parallel.ThreadPool will tend to perform 
            %               better than a parallel.ProcessPool because
            %               threads are allowed to share memory.
            %               Use 0 when operating on a GPU or if memory 
            %               usage explodes on a parallel.ProcessPool.
            % 
            % See also DAS BFDAS BFADJOINT

            % if not given a new Scan, use the UltrasoundSystem Scan
            if nargin < 3 || isempty(cscan), cscan = self.scan; end

            % ensure that the cscan is on Cartesian coordinates
            assert(isa(cscan, 'ScanCartesian'), "The medium Scan must be a ScanCartesian (given a " + class(cscan) + ").");

            % defaults
            kwargs.interp = 'linear';
            kwargs.parcluster = gcp('nocreate');
            kwargs.apod = 1;
            kwargs.keep_rx = false;
            kwargs.keep_tx = false;
            kwargs.bsize = 16;

            % set input options
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % get summation options
            sumtx = ~kwargs.keep_tx;
            sumrx = ~kwargs.keep_rx;

            % get cluster
            if isempty(kwargs.parcluster), kwargs.parcluster = 0; end % empty -> 0
            clu = kwargs.parcluster; % compute cluster/pool/threads
            % travel times cannot use threadPool because mex call
            if isa(clu, 'parallel.ThreadPool'), clu = 0; end

            % get worker transfer function
            if(clu == 0)
                 constfun = @(x) struct('Value', x);
            else,constfun = @parallel.pool.Constant;
            end

            % get the grid definition
            dp = [cscan.dx, cscan.dy, cscan.dz]; % step size each dim
            grd = {cscan.x, cscan.y, cscan.z}; % grid definition
            nsdims = (~isinf(dp)); % non-singleton dimensions
            dp = uniquetol(dp(nsdims)); % spatial step size
            ord = cscan.getPermuteOrder(); % order of grid
            sel = ord(ismember(ord, find(nsdims))); % selection of grid indices
            grd = grd(sel); % re-order and trim the grid
            
            % check the grid spacing is identical in non-singleton 
            % dimensions
            assert(numel(dp) == 1, ...
                'The cgrid must have equally sized steps in all non-singleton dimensions.'...
                );

            % check the data is FSA and matches transmit / receive
            % transducers
            assert(self.tx.numel == chd.M, 'Number of transmits must match number of transmitter elements.')
            assert(self.rx.numel == chd.N, 'Number of receives must match number of receiver elements.')

            % get the transmit, receive
            Pv = self.tx.positions();
            Pr = self.rx.positions();
            
            % get the sound speed in the Medium
            c = props(medium, cscan, 'c');

            % convert positions to sound speed grid coordinates (1-based)
            og = [cscan.x(1); cscan.y(1); cscan.z(1)];
            [Pvc, Prc] = dealfun(@(x) sub((x - og) ./ dp, sel, 1) + 1, Pv, Pr); % ([2|3] x [N|M])

            % enforce type and send data to workers maybe
            cnorm = constfun(double(gather(c ./ dp)));

            % get one-way delays within the field then generate samplers, using
            % reduced dimensions ([M|N] x {Cx x Cy x Cz})
            gi_opts = {'cubic', 'none'};
            tt = tic; fprintf('\nComputing Eikonal time delays ... \n');
            parfor (m = 1:chd.M, clu)
                % fprintf('tx %i\n', m);
                [tau_map_tx] = msfm(squeeze(cnorm.Value), double(Pvc(:,m))); %#ok<PFBNS> % travel time to each point
                tx_samp{m} = griddedInterpolant(grd, tau_map_tx,gi_opts{:}); %#ok<PFBNS> % make interpolator on cscan
            end
            parfor (n = 1:chd.N, clu)
                % fprintf('rx %i\n', n);
                [tau_map_rx] = msfm(squeeze(cnorm.Value), double(Prc(:,n))); %#ok<PFBNS> % travel time to each point
                rx_samp{n} = griddedInterpolant(grd, tau_map_rx,gi_opts{:}); %#ok<PFBNS> % make interpolator on cscan
            end
            fprintf('\nEikonal time delays completed in %0.3f seconds.\n', toc(tt));

            % get the imaging grid
            gi = self.scan.getImagingGrid(); % {I1 x I2 x I3} each
            gi = gi(sel); % select and trim dimensions 

            % splice args
            interp_method = kwargs.interp; 
            apod = kwargs.apod;

            % get sample times for each tx/rx
            tau_tx = cellfun(@(f) f(gi{:}), tx_samp, 'UniformOutput', false); % all receive delays
            tau_rx = cellfun(@(f) f(gi{:}), rx_samp, 'UniformOutput', false); % all receive delays
            tau_rx = cat(4, tau_rx{:}); % use all at a time
            tau_tx = cat(5, tau_tx{:}); % reconstruct matrix
            D = max(5, ndims(chd.data)); % data dimensions
            ord = [D+(1:3), chd.ndim, chd.mdim]; % apodization permutation order
            ord = [ord, setdiff(1:max(ord), ord)]; % full order (all dimes)

            % splice data, apod, delays per transmit
            [chds, ix] = splice(chd, chd.mdim, kwargs.bsize); % split into groups of data
            tau_tx = cellfun(@(ix){sub(tau_tx,ix,5)}, ix); % split per transmit block
            if size(apod,5) == 1, apod = {apod}; % store as a single cell
            else, apod = cellfun(@(ix){sub(apod,ix,5)}, ix); % split per transmit block
            end
            Mp = numel(chds); % number of transmit blocks
            
            % identify dimensiosn for summation
            sdim = [];
            if sumtx, sdim = [sdim, chd.mdim]; end % tx dimensions
            if sumrx, sdim = [sdim, chd.ndim]; end % rx dimensions
            
            b = 0; % initialize
            hw = waitbar(0,'Beamforming ...'); % create a wait bar
            parfor (m = 1:Mp, 0) % for each transmit (no cluster because parpool may cause memory issues here)
                % make the eikonal delays and apodization align with the channel data
                tau = tau_rx + tau_tx{m}; % I1 x I2 x I3 x N x M
                if isscalar(apod), a = apod{1}; else, a = apod{m}; end % recieve apodization per transmit (I1 x I2 x I3 x N x M)
                a = ipermute(a, ord); % move time delays / apodization into matching dimensions  
                tau = ipermute(tau, ord); % both: perm(1 x N x M) [x 1 x ... ] x I1 x I2 x I3
                    
                % sample and unpack the data (perm(1 x N x M) [x F x ... ] x I1 x I2 x I3)
                ym = sample(chds(m), tau, interp_method, a, sdim); % sample and sum over rx/tx maybe
                
                % sum or accumulate over transmit blocks
                if sumtx, b = b + ym; else, bm{m} = ym; end 
                if isvalid(hw), waitbar((Mp-m+1)/Mp, hw); end % update if not closed: parfor loops go backwards
            end
            if ~sumtx, b = cat(chd.mdim,bm{:}); end % combine if not summing tx
            if isvalid(hw), close(hw); end % close waitbar if not already closed

            % move image output into lower dimension
            if istall(b), b = gather(b); end % we have to gather tall arrays to place anything in dim 1
            b = swapdim(b, 1:3, D+(1:3)); % move image dimensions down (I1 x I2 x I3 [x F x ... ] x perm(1 x N x M))
        end
    
        function b = bfDAS(self, chd, c0, varargin)
            % BFDAS - Delay-and-sum beamformer
            %
            % b = BFDAS(self, chd, c0) creates a b-mode image b from the 
            % ChannelData chd and sound speed c0 (m/s). 
            % 
            % b = BFDAS(..., Name, Value, ...) passes additional Name/Value
            % pair arguments
            % 
            % 
            % Inputs:
            %   
            % keep_rx -     setting this to true returns an extra dimension
            %               containing the data for each receive prior to
            %               summation. The default is false.
            %
            % keep_tx -     setting this to true returns an extra dimension
            %               containing the data for each transmit prior to
            %               summation. The default is false.
            %
            % interp  -     specifies the method for interpolation. Support 
            %               is provided by the ChannelData/sample method. 
            %               The default is 'linear'.
            % 
            %   
            % apod    -     specifies and ND-array A for apodization. A must 
            %               be broadcastable to size  I1 x I2 x I3 x N x M 
            %               where I1 x I2 x I3 is the size of the image, N 
            %               is the number of receivers, and M is the number
            %               of transmits. The default is 1.
            % 
            %   
            % bsize   -     sets the ChannelData block size to at most B 
            %               transmits at a time in order to limit memory 
            %               usage. If memory is a concern, lowering B may 
            %               help. If there is ample memory available, a 
            %               higher value of B may yield better performance.
            % 
            %   
            % parcluster  - specifies a parcluster object for 
            %               parallelization. The default is the current 
            %               parallel pool returned by gcp.
            %               Setting clu = 0 avoids using a parallel cluster
            %               or pool. 
            %               A parallel.ThreadPool will tend to perform 
            %               better than a parallel.ProcessPool because
            %               threads are allowed to share memory.
            %               Use 0 when operating on a GPU or if memory 
            %               usage explodes on a parallel.ProcessPool.
            % 
            % See also DAS BFADJOINT CHANNELDATA/SAMPLE PARCLUSTER

            % defaults
            kwargs.interp = 'linear';
            kwargs.parcluster = gcp('nocreate');
            kwargs.apod = 1;
            kwargs.keep_rx = false;
            kwargs.keep_tx = false;
            kwargs.bsize = 16;

            % set input options
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % parse inputs
            if ~isnumeric(c0), c0 = c0.c0; end

            % get cluster
            clu = kwargs.parcluster;
            if isempty(clu) || isa(chd.data, 'gpuArray'), clu = 0; end % run on CPU by default

            % get image pixels, outside of range of data
            Pi = self.scan.getImagingGrid();
            Pi = cellfun(@(x) shiftdim(x, -1), Pi, 'UniformOutput',false);
            Pi = cat(1, Pi{:}); % 3 x I1 x I2 x I3

            % get the transmit, receive
            Pr = self.rx.positions(); % receiver positions

            % get virtual source or plane wave geometries
            switch self.sequence.type
                case 'FSA'
                    [Pv, Nv] = deal(self.tx.positions(), [0;0;1]);
                case 'VS'
                    [Pv, Nv] = deal(self.sequence.focus, [0;0;1]);
                case 'PW'
                    [Pv, Nv] = deal([0;0;0], self.sequence.focus); % TODO: use origin property in tx sequence
            end
            Pr = swapdim(Pr,2,5); % move N to dim 5
            [Pv, Nv] = deal(swapdim(Pv,2,6), swapdim(Nv,2,6)); % move M to dim 6

            % get receive delays
            dr = vecnorm(Pi - Pr,2,1); % 1 x I1 x I2 x I3 x N x 1
                
            % transmit sensing vector
            dv = Pi - Pv; % 3 x I1 x I2 x I3 x 1 x M
            switch self.sequence.type, 
                case {'VS', 'FSA'}, dv = vecnorm(dv, 2, 1) .* sign(sum(dv .* Nv,1));
                case {'PW'}, dv = sum(dv .* Nv, 1);
            end % 1 x I1 x I2 x I3 x 1 x M

            % bring to I1 x I2 x I3 x 1 x M
            [dv, dr] = deal(shiftdim(dv,1), shiftdim(dr,1));

            % make same type as the data
            [dv, dr] = dealfun(@(x)cast(x, 'like', real(chd.data)), dv, dr);

            % splice args
            apod = kwargs.apod; % apodization (I1 x I2 x I3 x N x M)
            interp_method = kwargs.interp; 
            cinv = 1 ./ c0;
            sumtx = ~kwargs.keep_tx;
            sumrx = ~kwargs.keep_rx;

            % move image dimensions beyond the data
            D = max([ndims(chd.data), ndims(dv), ndims(dr), 5]); % highest dimension of data
            [dv, dr, apod] = dealfun(@(x) swapdim(x, 1:3, D+(1:3)), dv, dr, apod);
            sdim = []; % dimensions to sum after apodization
            if sumrx, sdim = [sdim, chd.ndim]; end 
            if sumtx, sdim = [sdim, chd.mdim]; end 

            % TODO: fully support tall data by doing these permutations
            % in a different order to avoid permuting the tall
            % dimension perhaps by
            % * splicing in the tall dimension?
            % * moving the image dimensions beyond the data dimensions and
            %   operating there
            assert(ismember('T', chd.ord(1:3)), 'The time dimension must be in one of the first 3 dimensions.');
            tord = [chd.tdim, chd.ndim, chd.mdim]; % order to move to ChannelData dims
            
            % splice per transmit according to the block size
            [chds, is] = splice(chd, chd.mdim, kwargs.bsize); % splice transmit in chunks of size bsize
            dvm = cellfun(@(i){sub(dv,i,5)}, is); % splice per transmit
            if size(apod,5) ~= 1, 
                am = cellfun(@(i){sub(apod,i,5)}, is); % cell per transmit
            else
                am = {apod}; % single cell 
            end

            b = 0; % initialize
            % hw = waitbar(0,'Beamforming ...'); % create a wait bar
            for m = 1:numel(chds)
            % parfor (m = 1:numel(chds), clu) % for each transmit
                tau = cinv .* (dvm{m} + dr); % get sample times (1 x 1 x 1 x N x 1 x ... x I1 x I2 x I3)
                if isscalar(am), a = am{1}; else, a = am{m}; end % (1 x 1 x 1 x N x 1 x ... x I1 x I2 x I3)
                % a = sub(apod, min(m, Ma), 5); % recieve apodization 
                % move to permutation of (1 x N x M) - N/M aligned with ChannelData
                tau = swapdim(tau, [4,5], [chds(m).ndim, chds(m).mdim]); 
                a   = swapdim(a  , [4,5], [chds(m).ndim, chds(m).mdim]); 

                % sample, apodize, and sum over rx if requested
                z = sample(chds(m), tau, interp_method, a, sdim); % (perm(1 x N x M) x F x ... x I1 x I2 x I3)

                % swap the imaging and time/aperture dimensions
                z = ipermute(z, [tord, 4:gather(ndims(z))]); % 1 x N x M x F x ...
                z = swapdim(z, 1:3, D+(1:3)); % (I1 x I2 x I3 x F x ... x N x M)

                if sumtx, bm{m} = []; else; bm{m} = z; end % store tx to combine later (if requested)
                if sumtx, b = b + z; end % sum tx here (if requested)
                % if isvalid(hw), waitbar(m/M, hw); end % update if not closed
            end
            % if isvalid(hw), close(hw); end % close if not closed
            if ~sumtx, b = cat(D+3, bm{:}); end % combine all transmits
        end
    end
    
    % dependent methods
    methods
        function f = get.fc(self)
            if self.rx == self.tx
                f = self.rx.fc;
            else
                f = [self.tx.fc, self.rx.fc];
            end
        end
        function set.fc(self, f)
            if self.rx == self.tx
                self.rx.setCentralFrequency(f);
            else
                self.tx.setCentralFrequency(f(1));
                self.rx.setCentralFrequency(f(end));
            end
        end
        function x = get.xdc(self)
            if self.tx == self.rx
                x = self.rx;
            else
                error("Call to 'xdc' ambiguous; transmitter and receiver are not the same.");
            end
        end
        function set.xdc(self, xdc), [self.tx, self.rx] = deal(xdc); end
        function c = get.pulse(self), c = self.sequence.pulse; end
        function set.pulse(self, c), self.sequence.pulse = c; end
    end
    
    % recompilation helper functions
    methods
        function recompile(self), recompileMex(self); if gpuDeviceCount, recompileCUDA(self); end, end
        % RECOMPILE - Recompile mex and CUDA files
        %
        % RECOMPILE(self) recompiles all mex binaries and CUDA files and 
        % stores them in self.tmp_folder. If there are no MATLAB compatible
        % GPUs, it does not attempt to recompile CUDA files.
        %
        % See also ULTRASOUNDSYSTEM/RECOMPILECUDA ULTRASOUNDSYSTEM/RECOMPILEMEX
        function recompileMex(self)
            % RECOMPILEMEX - Recompile mex files
            %
            % RECOMPILEMEX(self) recompiles all mex binaries and stores
            % them in self.tmp_folder.
            %
            % See also ULTRASOUNDSYSTEM/RECOMPILECUDA
            % ULTRASOUNDSYSTEM/RECOMPILE
            
            
            defs = UltrasoundSystem.getMexFileDefs();

            % compile each definition
            for d = defs
                % make full command
                com = cat(1,...
                    '-outdir', self.tmp_folder, ... place binaries in system's tmp folder
                    join("-I" + d.IncludePath), ...
                    join("-L" + d.Libraries), ...
                    join("-D" + d.DefinedMacros),...
                    fullfile(UltrasoundSystem.getSrcFolder(), 'FMM', 'functions', d.Source) ...
                    );
                
                try
                    com = cellstr(com);
                    s = mex(com{:});
                    if s, warning("Error recompiling code!"); else, disp("Success recompiling " + d.Source); end
                catch
                    warning("Unable to recompile code!");
                end
            end
   
        end
        function recompileCUDA(self)
            % RECOMPILECUDA - Recompile CUDA ptx files
            %
            % RECOMPILECUDA(self) recompiles all CUDA files and stores
            % them in self.tmp_folder.
            %
            % See also ULTRASOUNDSYSTEM/RECOMPILEBFCONST
            % ULTRASOUNDSYSTEM/RECOMPILE ULTRASOUNDSYSTEM/RECOMPILEMEX
            
            % src file folder
            src_folder = UltrasoundSystem.getSrcFolder();
            
            % get all source code definitions
            defs = UltrasoundSystem.getCUDAFileDefs();
            
            % compile each
            for d = defs
                % make full command
                com = join(cat(1,...
                    "nvcc ", ...
                    "--ptx " + fullfile(src_folder, d.Source), ...
                    "-arch=native ", ... % for half types: TODO move to compile option
                    "-o " + fullfile(self.tmp_folder, strrep(d.Source, 'cu', 'ptx')), ...
                    join("--" + d.CompileOptions),...
                    join("-I" + d.IncludePath), ...
                    join("-L" + d.Libraries), ...
                    join("-W" + d.Warnings), ...
                    join("-D" + d.DefinedMacros)...
                    ));
                
                try
                    s = system(com);
                    if s, warning("Error recompiling code!"); else, disp("Success recompiling " + d.Source); end
                catch
                    warning("Unable to recompile code!");
                end
            end
        end
        function recompileBFCONST(self, chd)
            % RECOMPILEBFCONST - Constant size compilation for beamforming
            %
            % RECOMPILEBFCONST(self) recompiles the beamforming cuda 
            % executables for the current scan, transducer, and 
            % transmit sequence. 
            % 
            % Using a fixed data size triggers compiler-level optimizations
            % that can improve performance for iterative calls.
            %
            % After calling this function, calls with a different size of 
            % data will have unexpected results. Use
            % UltrasoundSystem/recompileCUDA to reset the binaries to
            % handle variable sizes.
            % 
            % RECOMPILEBFCONST(self, chd) additionally uses a fixed size
            % ChannelData object.
            % 
            % See also ULTRASOUNDSYSTEM/RECOMPILECUDA 
            % ULTRASOUNDSYSTEM/RECOMPILE 
            
            % src file folder
            src_folder = UltrasoundSystem.getSrcFolder();
            
            % get the other sizes for beamform.m
            VS = ~(self.sequence.type == "PW"); % whither virtual source
            Isz = self.scan.size; % size of the scan
            N = self.xdc.numel; % number of receiver elements
            M = self.sequence.numPulse; % number of transmits
            assert(~isnan(M) && M > 0); % M must be defined
            
            % get all source code definitions
            defs = UltrasoundSystem.genCUDAdef_beamform();
            
            % add the defined macros
            defs.DefinedMacros = cat(1, ...
                defs.DefinedMacros, ... keep the current defs
                "QUPS_" + {... prepend 'QUPS_'
                "VS="+VS,... virtual source model
                "N="+N,... elements
                "M="+M,... transmits
                "I1="+Isz(1),... pixel dim 1
                "I2="+Isz(2),... pixel dim 2
                "I3="+Isz(3) ... pixel dim 3
                }');
            
            % if T is provided, include it
            if nargin >= 2
                defs.DefinedMacros = [defs.DefinedMacros; ...
                    {"T="+chd.T}; ... number of time samples
                    ];
            end
            
            % compile each
            for d = defs
                % make full command
                com = join(cat(1,...
                    "nvcc --ptx ", ...
                    "-arch=native ", ... % for half types: TODO move to compile option
                    fullfile(src_folder, d.Source) + " ", ...
                    "-o " + fullfile(self.tmp_folder, strrep(d.Source, 'cu', 'ptx')), ...
                    join("--" + d.CompileOptions),...
                    join("-I" + d.IncludePath), ...
                    join("-L" + d.Libraries), ...
                    join("-W" + d.Warnings), ...
                    join("-D" + d.DefinedMacros)...
                    ));
                
                try
                    s = system(com);
                    if s, warning("Error recompiling code!"); else, disp("Success recompiling " + d.Source); end
                catch
                    warning("Unable to recompile code!");
                end
            end
        end
    end
    
    % source file recompilation definitions
    methods(Static,Hidden)
        function defs = getCUDAFileDefs()
            % no halp :(
            
            % get all source code definitions
            defs = [...
                UltrasoundSystem.genCUDAdef_beamform(),...
                UltrasoundSystem.genCUDAdef_interpd(),...
                UltrasoundSystem.genCUDAdef_convd(),...
                UltrasoundSystem.genCUDAdef_greens(),...
                ];
        end

        function defs = getMexFileDefs()
            % no halp :(
            
            % get all source code definitions
            defs = [...
                UltrasoundSystem.genMexdef_msfm(), ... % both msfm files
                ];
        end

        function f = getSrcFolder()
            f = fileparts(mfilename('fullpath'));
        end
        
        function d = genCUDAdef_beamform()
            % no halp :(

            % filename
            d.Source = {...
                'bf.cu', ...
                }';

            d.IncludePath = {}; % include folders
            d.Libraries = {}; % libraries

            d.CompileOptions = {...  compile options
                'use_fast_math',...
                };

            d.Warnings = {... warnings
                'no-deprecated-gpu-targets'...
                };
            
            d.DefinedMacros = {};
        end

        function d = genCUDAdef_interpd()
            % no halp :(

            % filename
            d.Source = {...
                'interpd.cu', ...
                }';

            d.IncludePath = {}; % include folders
            d.Libraries = {}; % libraries

            d.CompileOptions = {...  compile options
                'use_fast_math',...
                };

            d.Warnings = {... warnings
                'no-deprecated-gpu-targets'...
                };
            
            d.DefinedMacros = {};
        end

        function d = genCUDAdef_greens()
            % no halp :(

            % filename
            d.Source = {...
                'greens.cu', ...
                }';

            d.IncludePath = {}; % include folders
            d.Libraries = {}; % libraries

            d.CompileOptions = {...  compile options
                'use_fast_math',...
                };

            d.Warnings = {... warnings
                'no-deprecated-gpu-targets'...
                };

            d.DefinedMacros = {};
        end

        function d = genCUDAdef_convd()
            % no halp :(

            % filename
            d.Source = {...
                'conv_cuda.cu', ...
                }';

            d.IncludePath = {}; % include folders
            d.Libraries = {}; % libraries

            d.CompileOptions = {...  compile options
                'use_fast_math',...
                };

            d.Warnings = {... warnings
                'no-deprecated-gpu-targets'...
                };
            
            d.DefinedMacros = {};
        end

        function d = genMexdef_msfm()
            %

            % generate for msfms2d
            d.Source = 'msfm2d.c';

            d.IncludePath = {... include folders
                fullfile(UltrasoundSystem.getSrcFolder(), 'FMM', 'functions'), ... helper_math
                };

            d.Libraries = {...
                ...
                };

            d.DefinedMacros = {...
                };

            d = repmat(d, [1, 2]);

            % generate for msfms3d - same files, different structure
            d(2).Source = 'msfm3d.c';
        end
    end
end

