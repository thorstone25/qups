classdef UltrasoundSystem < handle
    
    % objects
    properties
        tx
        rx
        sequence
        scan
    end
    
    % parameters
    properties
        fs = 40e6           % system sampling frequency (Hz)
    end
    
    properties(Dependent)
        fc                  % central operating frequency (Hz)
        xdc                 % if a single transducer
        pulse               % transmit signal, if only one
    end
    
    properties(Hidden,SetAccess=protected)
        tmp_folder
    end
        
    % get/set & constructor
    methods
        % constructor
        function self = UltrasoundSystem(varargin)
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
            defs = self.getCUDAFileDefs();
            fls = arrayfun(@(d) string(strrep(d.Source, 'cu', 'ptx')), defs);
            s = arrayfun(@(fl) copyfile(which(fl), fullfile(self.tmp_folder, fl)), fls);
            if any(~s), self.recompileCUDA(); end % attempt to recompile code
            
            % copy code or recompile it
            % TODO: generalize to mex extension for other machines (with 
            % isunix or iswindows or ismac)
            defs = self.getMexFileDefs();
            fls = arrayfun(@(d) string(strrep(d.Source, 'c', 'mexa64')), defs);
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
            % sub_div = GETLAMBDASUBDIV(self, p, c, ap)
            % returns the element subdivision vector corresponding to a
            % proportion p of the wavelength given a medium or speed of
            % sound. In other words, for 
            
            if(nargin < 4), ap = 'rx'; end
            if isa(c, 'Medium'), c = c.c0; end
            
            % get wavelength
            lam = c / self.fc;
            
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
        function [chd, wv] = comp_RS_FSA(self, target, element_subdivisions, varargin)
            % COMP_RS_FSA - Compute rayleigh-Sommerfield via full-synthetic-aperture
            % 
            % chd = COMP_RS_FSA(self, target)
            % computes the full synthetic aperture data using a simple
            % Green's function kernel applied to all sub elements and all
            % point scatterers. It then applies the transmit sequence to
            % form the channel data.
            %
            % chd = COMP_RS_FSA(self, target, element_subdivisions) uses
            % multiple sub-apertures in the integration. This argument is
            % passed to FieldII to construct the subdivisions. FieldII must
            % be on the path.
            %
            % chd = COMP_RS_FSA(..., Name, Value, ...) provides additional
            % name/value pair arguments.
            %
            % [chd, wv] = COMP_RS_FSA(...) additionally returns the final
            % waveform convolved across the point scatterers and aperture.
            %
            % Inputs:
            %   - target:               a Target object
            %   - element_subdivisions: 1 x 2 vector of subdivisions per
            %                           transducer element (optional)
            %
            % Name-Value Pairs
            %   - device    integer representing gpu selection. -1 selects
            %   the current GPU. 0 selects a cpu. n where n > 0 selects and
            %   resets a GPU
            %
            %   - method    interpolation method {'interpn'* | 'griddedInterpolant'}
            %
            % Outputs:
            %   - time (T x 1):             time sample values (s)
            %   - voltages (T x M x N):     voltage samples (pure)
            %
            % where T -> time, M -> transmitters, N -> receivers
            %
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_ALL
            
            % get Tx/Rx apertures subelement positions
            % (3 x N x E)
            if nargin < 3, element_subdivisions = self.getLambdaSubDiv(0.1, target); end
            device = -logical(gpuDeviceCount());
            method = 'interpn';
            interp_method = 'linear';
            
            for i = 1:2:numel(varargin)
                switch lower(varargin{i})
                    case 'device', device = varargin{i+1};
                    case 'method', method = varargin{i+1};
                    case 'interp', interp_method = varargin{i+1};
                end
            end
            
            % get the centers of all the sub-elements
            if all(element_subdivisions == 1) % no sub-elements
                ptc_tx = self.tx.positions();
                ptc_rx = self.rx.positions();
            else % use FieldII's definitions
                ptc_tx = self.tx.getFieldIIBaryCenters(element_subdivisions);
                ptc_rx = self.rx.getFieldIIBaryCenters(element_subdivisions);
            end
            
            % get the Tx/Rx impulse response function ( T' x 1)
            tx_impulse_waveform = self.tx.impulse;
            t_tx = tx_impulse_waveform.getSampleTimes(self.fs);
            tx_impulse_samples = reshape(tx_impulse_waveform.fun(t_tx), 1, []);
            
            rx_impulse_waveform = self.rx.impulse;
            t_rx = rx_impulse_waveform.getSampleTimes(self.fs);
            rx_impulse_samples = reshape(rx_impulse_waveform.fun(t_rx), 1, []);
            
            % set the waveform
            %%% TODO: handle multiple chirps in the sequence ... or not?
            signal_waveform = self.sequence.pulse.copy(); % leave the original intact
            t_sig = signal_waveform.getSampleTimes(self.fs);
            signal_samples  = reshape(signal_waveform.fun(t_sig), 1, []);

            % get maximum necessary time sample (use manhattan distance and
            % sound speed to give an upper bound)
            maxdist = @(p) max(vecnorm(p,2,1), [], 'all');
            mindist = @(p) min(vecnorm(p,2,1), [], 'all');
            taumax = (2 * maxdist(target.pos) + maxdist(ptc_tx) + maxdist(ptc_rx)) ./ min(target.c0);
            taumin = (2 * mindist(target.pos) + mindist(ptc_tx) + mindist(ptc_rx)) ./ max(target.c0);
            tstart = t_tx(1) + t_rx(1) + t_sig(1); 
            
            % get the convolved kernel
            kern = conv(conv(tx_impulse_samples, rx_impulse_samples, 'full'), signal_samples, 'full');
            tk = tstart + (0 : 1 : numel(kern) - 1 )' ./ self.fs;
            wv = Waveform('samples', kern, 't', tk);

            % get minimum/maximum sample times
            tmin = taumin + wv.t0;
            tmax = taumax + wv.tend;

            % create time vector (T x 1)
            % this formulation is guaranteed to pass through t == 0
            t = (floor(tmin * self.fs) : ceil(tmax * self.fs))';% ./ self.fs;

            % pre-allocate output
            [T, N, M, E] = deal(numel(t), self.rx.numel, self.tx.numel, prod(element_subdivisions));
            x   = zeros([T N M]);

            % splice
            c0  = target.c0;
            pos = target.pos; % 3 x S
            amp = target.amp; % 1 x S
            fs_ = self.fs;

            % TODO: set data types | cast to GPU? | perform on GPU?
            if device > 0, gpuDevice(device); end
            if device, 
                [pos, ptc_rx, ptc_tx, t, tk, kern] = dealfun(@gpuArray, pos, ptc_rx, ptc_tx, t, tk, kern);
            end

            % create 1D interpolator function
            switch method
                case 'griddedInterpolant'
                    % terp = griddedInterpolant(gather(tk), gather(kern), 'linear','none');
                    terp = griddedInterpolant(gather(round(tk*fs_)), gather(kern), interp_method,'none');
                    f = @(tau) (terp(gather(t - tau)));
                case 'interpn'
                    f = @(tau) interpn(round(fs_*tk), kern, t - tau, interp_method, 0);
                case 'interpd'
                    f = @(tau) interpd(kern(:), (t - tau) - (fs_*tk(1)), 1, interp_method);
                otherwise
                    f = @(tau) interpn(fs_*tk, kern, t - tau, interp_method, 0);
            end

            % cast dimensions to compute in parallel
            ptc_rx = permute(ptc_rx, [1,6,2,4,3,5]); % 3 x 1 x N x 1 x En x 1
            ptc_tx = permute(ptc_tx, [1,6,4,2,5,3]); % 3 x 1 x 1 x M x  1 x Em

            % for each tx/rx pair, extra subelements ...
            % TODO: parallelize where possible
            w = waitbar(0);
            for s = 1:target.numScat
            for em = 1:E
                for en = 1:E
                    % unpack indexing (MATLAB optimization)
                    % [en, em, n, m] = ind2sub([E, E, N, M], i);

                    % compute time delays
                    % TODO: do this via ray-path propagation through a
                    % medium
                    % 1 x S x N x M x 1 x 1
                    r_rx = vecnorm(sub(pos,s,2) - sub(ptc_rx, en, 5));
                    r_tx = vecnorm(sub(pos,s,2) - sub(ptc_tx, em, 6)); 
                    tau_rx = r_rx ./ c0;
                    tau_tx = r_tx ./ c0;

                    % compute the contribution
                    % T x S x N x M x 1 x 1
                    att = sub(amp,s,2);% .* (1 ./ r_rx) .* (1 ./ r_tx); % propagation attenuation
                    s_ = att .* f((tau_rx + tau_tx) * fs_); % temporal delay (no phase shift)

                    % add contribution (T x N X M)
                    x(:) = x + permute(nan2zero(sum(s_, 2, 'omitnan')), [1,3,4,2]);

                    % update waitbar
                    if isvalid(w), waitbar(sub2ind([E,E],en,em) / (E*E), w); end
                end
            end
            end
            if isvalid(w), delete(w); end

            % make a channel data object
            chd_ = ChannelData('t0', sub(t ./ self.fs,1,1), 'fs', self.fs, 'data', x);

            % synthesize linearly
            [chd] = self.focusTx(chd_, self.sequence, 'interp', interp_method);

            % send data (back) to CPU
            chd = gather(chd);
        end
    end

    % SIMUS calls
    methods
        function chd = simus(self, target, varargin)
            % no halp :(
            
            % defaults
            % TODO: forward arguments to params or opt as appropriate
            kwargs = struct(...
                'device', 0, ...
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
            [X, Y, Z, A] = deal(sub(target.pos,1,1), sub(target.pos,2,1), sub(target.pos,3,1), target.amp);
            if isempty(kwargs.dims) 
                if all(Y == 0, 'all'), kwargs.dims = 2; Y = []; % don't simulate in Y if it is all zeros 
                else, kwargs.dims = 3; end
            end
            if kwargs.dims == 2 && any(Y ~= 0)
                warning("QUPS:UltrasoundSystem:simus:casting", "Projecting all points onto Y == 0 for a 2D simulation.");
            end

            % get all other param struct values (implicitly force same
            % transducer)
            p = {getSIMUSParam(target), getSIMUSParam(self.xdc)};
            p = cellfun(@struct2nvpair, p, 'UniformOutput', false);
            p = cat(1, p{:});
            p = struct(p{:});

            % set transmit sequence ... the only way we can
            % TODO: forward arguments to transmit parameters
            p.fs    = self.fs;
            p.TXnow = kwargs.periods; % number of wavelengths
            p.TXapodization = zeros([self.xdc.numel,1]); % set tx apodization
            p.RXdelay = zeros([self.xdc.numel,1]); % receive delays (none)

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

            % call the sim: FSA approach
            M = self.xdc.numel; % splice
            parfor (m = 1:M, kwargs.parcluster)
                p_ = p; % splice
                p_.TXapodization(m) = 1; % transmit only on element m
                rf{m} = simus(X,Y,Z,A,zeros([M,1]),p_,opt); % rf for each transmit (no delays)
            end

            % extend time axis so that it's identical for all transmits
            T = (cellfun(@(x)size(x,1), rf));
            if numel(unique(T)) > 1
                rf = cellfun(@(x) cat(1, x, zeros([max(T) - size(x,1), size(x,2:ndims(x))], 'like', x)), rf, 'UniformOutput', false);
            end
            rf = cat(3, rf{:}); % T x N x M

            % create the output QUPS ChannelData object
            chd = ChannelData('data', rf, 't0', 0, 'fs', self.fs);

            % synthesize linearly
            chd = self.focusTx(chd, self.sequence, 'interp', kwargs.interp);
        end
    end

    % Field II calls
    methods(Access=public)
        function [chd] = calc_scat_all(self, target, element_subdivisions, varargin)
            % [chd] = calc_scat_all(self, target, element_subdivisions)
            % 
            % function to compute the full synthetic aperture data using
            % Field II.
            %
            % Inputs:
            %   - target:               a Target object
            %   - element_subdivisions: 2 x 1 vector of subdivisions per
            %                           transducer element (optional)
            %
            % Outputs:
            %   - time (T x 1):             time sample values (s)
            %   - voltages (T x M x N):     voltage samples (pure)
            % where T -> time, M -> transmitters, N -> receivers
            
            % defaults
            kwargs = struct('interp', 'linear');

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
            
            % set sound speed/sampling rate
            set_field('fs', self.fs);
            set_field('c',target.c0);
            
            % get Tx/Rx apertures
            if nargin < 3, element_subdivisions = self.getLambdaSubDiv(0.1, target); end
            p_focal = mean(target.pos,2);
            Tx = self.tx.getFieldIIAperture(p_focal.', element_subdivisions);
            Rx = self.rx.getFieldIIAperture(p_focal.', element_subdivisions);
            
            % set the Tx/Rx impulse response function
            tx_impulse_waveform = self.tx.impulse();
            t_tx = tx_impulse_waveform.getSampleTimes(self.fs);
            tx_impulse_samples = reshape(tx_impulse_waveform.fun(t_tx), 1, []);
            xdc_impulse(Tx, real(tx_impulse_samples));
            
            rx_impulse_waveform = self.rx.impulse();
            t_rx = rx_impulse_waveform.getSampleTimes(self.fs);
            rx_impulse_samples = reshape(rx_impulse_waveform.fun(t_rx), 1, []);
            xdc_impulse(Rx, real(rx_impulse_samples));
            
            % set the waveform
            %%% TODO: handle multiple chirps in the sequence
            signal_waveform = self.sequence.pulse.copy(); % leave the original intact
            signal_samples  = reshape(signal_waveform.toSampled(self.fs), 1, []);
            xdc_excitation(Tx, real(signal_samples));
            
            % call the sim
            down_sampling_factor = 1;
            [voltages, start_time] = calc_scat_all(Tx, Rx, ...
                target.pos.', target.amp.', down_sampling_factor);
            
            % number of samples in time
            T = size(voltages,1);
            
            % create the time vector
            waveform_start = signal_waveform.t0 + t_tx(1) + t_rx(1);
            time = (0:T-1)' * (1 / self.fs) + (start_time + waveform_start);
            
            % reduce precision to save RAM
            time = single(time);
            voltages = single(voltages);
            
            % reshape to T x N x M
            M = self.tx.numel;
            N = self.rx.numel;
            voltages = reshape(voltages, [T N M]);

            % create the output QUPS ChannelData object
            chd_ = ChannelData('data', voltages, 't0', time(1), 'fs', self.fs);

            % synthesize linearly
            [chd] = self.focusTx(chd_, self.sequence, 'interp', kwargs.interp);

            % cleanup
            if field_started, evalc('field_end'); end
        end        
    end
    
    % k-Wave calls
    methods(Access=public)
        function [kgrid, PML_size, kgrid_origin, kgrid_size, kgrid_step, kmedium] = getkWaveGrid(self, target, varargin)
            % [kgrid, PML_size, origin, grid_size, step] = getkWaveGrid(self, target)
            %
            % function to create a kWaveGrid object for a given Target.
            %
            % Inputs:
            %   - target:       a Target object
            %
            % Name-Value Pair Inputs:
            %   - dims:         number of dimensions (default = 2)
            %
            %   - PML_min:      minimum PML size (default = 4)
            %
            %   - PML_max:      minimum PML size (default = 48)
            %
            %   - CLF_max:      maximum CFL: the kgrid will use a physical
            %                   and temporal spacing beneath this value
            %                   that is a ratio of the sampling frequency
            %                   (default = 0.25)
            %
            %   - buffer (3 x 2):   3 x 2 array specifying computational 
            %                   buffer to provide spacing from the other 
            %                   objects including the transmitter, 
            %                   receiver, and target objects. This is the 
            %                   amount of room that the grid is expanded 
            %                   by. Given in meters as 
            %                   [xlo, xhi; ylo, yhi, zlo, zhi] 
            %                   (default = repmat(5e-3, [3,2]) )
            %
            %   - resolution-ratio: sets the ratio of grid spacing to
            %                   spatial wavelength assuming a given a 
            %                   reference speed of sound from the target
            %                   (default = 0.25)
            %
            %   - reference-sound-speed: overrides the reference sound
            %                   speed for the resolution ratio. 
            %                   (default = target.c0)
            %
            % Outputs:
            %   - kgrid:            a kWaveGrid object
            %
            %   - PML_size (3 x 1): the chosen PML size that minimizes the
            %                       maximum prime number within PML-min and
            %                       PML-max
            %
            %   - origin (3 x 1):   cartesian grid origin to recover the
            %                       original coordinates when using
            %                       kgrid.{x|y|z} or similar functions
            %
            %   - grid_size (3 x 1): size of the computational grid, 
            %                       without the PML. The size is 1 for
            %                       sliced dimensions
            %
            %   - step (3 x 1):     descritization size in each dimension
            %                       in meters. The size is inf for sliced
            %                       dimensions
            
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
            lam = c0 ./ (max(self.xdc.fc) + self.xdc.bw*[-1,1]/2);
            [dx, dy, dz] = deal(lam_res * min(lam));
            dp = [dx;dy;dz]; % vector difference
            
            % finds the smallest PML buffer from PML_min_size to search_max that
            % minimizes the factoring size
            PML_buf = @(n) (argmin(arrayfun(@(a)max(factor(n+2*a)), PML_min_size:PML_max_size)) + PML_min_size - 1);
            
            % get total min/max bounds for the transducers and the target
            pb_t = target.getBounds(); % get min/max bounds of the target
            pb_tx = self.tx.bounds(); % get min/max bounds of the tx
            pb_rx = self.rx.bounds(); % get min/max bounds of the tx
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
            % ULTRASOUNDSYSTEM/GETKWAVESOURCE - Create a k-wave compatible source structures
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
            %

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
                if(~rx_imp_is_delta) % rx impulse not a delta: delay the impulse
                    imp_samp = rx_imp.sample(t_rx - delays); %#ok<PFBNS> % (T_data x nSubEl)
                    resp_samp = convd(resp_samp, imp_samp, 1, 'full', 'device', device); % (T_resp x nSubEl)
                else % impulse is a delta: use linear resampling (sketchy)
                    warning('UltrasoundSystem:kWaveSensor2ChannelData:convolvingDeltas', ...
                        ['Using linear interpolation to convolve two delta functions with a shift. ', ...
                        'Use a small rect function to avoid this. ']);
                    imp_samp = vec(rx_imp.sample(t_rx)); % (T_data x nSubEl)
                    conv_samp = convn(resp_samp, imp_samp, 'full'); % (T_resp x nSubEl) % GPU OOM?
                    conv_sampler = griddedInterpolant(t_resp, zeros([T_resp,1], 'like', t_resp), interp, 'none');
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
    
        function [chd, cgrid]= kspaceFirstOrderND(self, target, element_subdivisions, varargin)
            %

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
    
    % Beamforming
    methods(Access=public)
        function [B, X, Y, Z] = DAS(self, chd, medium, rcvfun, varargin)
            % Delay and sum
            %
            % [B, X, Y, Z] = DAS(self, chd, medium, rcvfun, varargin)
            % 
            % function to perform delay and sum beamforming on SAR data
            % inputs:
            %   - time (T x 1) or {2}:  time sample values or a start time, 
            %                           frequency pair (s) | (s, Hz)
            %   - resp (T x M x N):     voltage samples (pure)
            %   - medium:               Medium object
            %   - rcvfun:               Receive aperture accumulation
            %                           (defaults to a summation)
            %   Name-value pair arguments
            %   - prec:                 compute precision of the positions
            %                           {'single'* | 'double'}
            %   - device:               GPU device index: -1 for default
            %                           gpu, 0 for cpu, 1-4 for device index
            % where T -> time, M -> transmitters, N -> receivers
            %
            % outputs:
            %   - X\Y\Z:    3D coordinates of the output
            %   - B:        B-mode image
            
            % default name-value pair arguments
            prec = 'single';
            device = int64(-1 * logical(gpuDeviceCount)); % {0 | -1} -> {CPU | GPU}
            apod = 1;
            interp_args = {};
            
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
                    otherwise
                        error('Unrecognized name-value pair');
                end
            end
            
            % get positions of the imaging plane % 
            [X, Y, Z, image_size] = self.scan.getImagingGrid();

            % get the apodization
            % apod = self.scan.apodScanline(self.sequence);

            % reshape into I x M x N?
            apod_args = {'apod', apod};
            
            % convert to x/y/z in 1st dimension
            P_im = permute(cat(4, X, Y, Z),[4,1,2,3]); % 3 x I1 x I2 x I3 = 3 x I
            
            % get positions of the aperture(s)
            P_tx = self.tx.positions(); % cast(self.tx.positions(), 'like', time(end)); % 3 x M
            P_rx = self.rx.positions(); % cast(self.rx.positions(), 'like', time(end)); % 3 x N
            
            % get the beamformer arguments
            dat_args = {chd.data, chd.t0, chd.fs, medium.c0, 'device', device, 'position-precision', prec}; % data args
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
            if nargin < 5 || isempty(rcvfun)
                % beamform the data (I1 x I2 x I3 x 1)
                B = beamform('DAS', pos_args{:}, dat_args{:}, ext_args{:});
            else
                % beamform the data (I1 x I2 x I3 x N)
                B = beamform('SYN', pos_args{:}, dat_args{:}, ext_args{:});
                
                % apply recieve aperture function in dim 1 and shape up
                B = rcvfun(B, 4); % I1 x I2 x I3 x 1
            end
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
            %  sum -                  Set false to preserve transmits
            %  interp -               Interpolation method
            %
            %  Interpolation methods are passed to the ChannelData/sample
            %  function. An additional method 'freq' is provided to delay
            %  the data in the frequency domain.
            % 
            % where T -> time, M -> transmitters, N -> receivers, M' ->
            % synthetic transmits
            %
            % Outputs:
            %   chd -                   New ChannelData object
            %
            % See also CHANNELDATA/SAMPLE
            
            % defaults
            L = [];
            summation = true;
            interp_method = 'linear';
            
            % focus using self's sequence 
            if nargin < 3 || isempty(seq), seq = self.sequence; end

            % optional inputs
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'length', L = varargin{i+1};
                    case 'sum', summation = varargin{i+1};
                    case 'interp', interp_method = varargin{i+1};
                    otherwise, error('Unrecognized input.'); 
                end     
            end

            % nothing to do for FSA acquisitions
            switch seq.type, case 'FSA', chd = copy(chd0); return; end 

            % dist/time to receiver
            tau_focal = - shiftdim(seq.delays(self.tx),      -2); % 1 x 1 x M x M'
            apod      =   shiftdim(seq.apodization(self.tx), -2); % 1 x 1 x [1|M] x [1|M']
            
            % get time vector description
            tstart = chd0.t0;
            dt = 1 ./ chd0.fs;
            
            % find minimum signal length
            nwrap = ceil(max(abs(tau_focal(:))) ./ dt); % amount the signal extends beyond original
            tend = max(chd0.t0) + dt*(chd0.T - 1 + nwrap); % maximum signal time
            S = ceil((tend - tstart) ./ dt) + 1; % minimum signal length
            
            % pick new signal length
            if isempty(L)
                L = S + 1;
            elseif ischar(L)
                switch L
                    case 'min'
                        L = S + 1;
                    case 'pow2'
                        try 
                            L = 2^(nextpow2(S));
                        catch
                            L = gpuArray(2^nextpow2(gather(S)));
                        end
                    otherwise
                        error('Unrecognized choice for DFT length');
                        
                end
            elseif isscalar(L)
                if L < S, warning('Signal length may be too short!'); end
            else
                error("L must be a scalar or one of {min | pow2}");
            end
            
            % extended time vector
            l = shiftdim(0 : 1 : (L - 1), 1); % L x 1 x 1

            % choose domain to process data in
            switch interp_method
                case 'freq'
                    wL = apod .* exp(-2i * pi ./ dt ./ L .* l .* tau_focal); % steering vector (L x 1 x M x M')
                    x = fft(chd0.data, L, 1); % data in freq domain (L x N x M x 1)
                    if summation, y = tenmul(wL,   x, 3, 1); % apply phase shift and sum over transmits (L x N x 1 x M')
                    else,         y =        wL .* x;       % apply phase shift, no sum over transmits (L x N x M x M') <- huge!!!
                    end
                    z = ifft(y, L, 1, 'nonsymmetric'); % back to time (L x N x [1|M] x M')
                    z = circshift(z, nwrap, 1); % handle FFT wrap-around

                otherwise
                    time = cast(dt .* l + tstart, 'like', chd0.t0); % T' x 1 x 1
                    for n = chd0.N:-1:1 % per receive channel (implicit pre-allocation)
                        chd_ = copy(chd0); % don't modify original data
                        chd_.data = sub(chd_.data,n,2); % choose only this channel
                        y = apod .* chd_.sample(time - tau_focal - dt*nwrap, interp_method); % sample and apodize the data (T x 1 x M x M')
                        if summation, y = sum(y, 3); end
                        z(:,n,:,:) = y; % store data (T x 1 x [1|M] x M')
                    end
            end

            % create output channel data object
            chd = copy(chd0); % cpy all properties
            chd.t0 = tstart - dt*nwrap; % shift time axes
            chd.data = permute(z,[1,2,4,3]); ... % convert to L x N x M' x [1|M]
        end
        
        function [B, X, Y, Z] = DASEikonal(self, chd, medium, cgrid, rcvfun, varargin)
            % function to perform delay and sum beamforming on SAR data
            % inputs:
            %   - time (T x 1) or {2}:  time sample values or a start time, 
            %                           frequency pair (s) | (s, Hz)
            %   - resp (T x M x N):     voltage samples (pure)
            %   - medium:               Medium object
            %   - rcvfun:               Receive aperture accumulation
            %                           (defaults to a summation)
            %   - cgrid:                A structure defining the grid on
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
            %   Name-value pair arguments
            %   - prec:                 compute precision of the positions
            %                           {'single'* | 'double'}
            %   - device:               GPU device index: -1 for default
            %                           gpu, 0 for cpu, 1-4 for device index
            % where T -> time, M -> transmitters, N -> receivers
            %
            % outputs:
            %   - X\Y\Z:    3D coordinates of the output
            %   - B:        B-mode image
            
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
    
        function b = matrixbf(self, chd, c0)
            % MATRIXBF Matrix beamformer
            %
            % b = matrixbf(self, chd)
            % 
            % 
            % See also ULTRASOUNDSYSTEM/DAS ULTRASOUNDSYSTEM/FOCUSTX

            % options
            kwargs.fthresh = -40; % threshold for including frequencies

            % move the data to the frequency domain, being careful to
            % preserve the time axis
            K = round(chd.T); % DFT length
            f = chd.fs * (0 : K - 1)' / K; % frequency axis
            df = chd.fs * 1 / K; % frequency step size
            x = fft(chd.data,K,1); % K x N x M x ...
            x = x .* exp(-2i*pi*f.*chd.t0); % re-align t0 axis

            % choose frequencies to evaluate
            xmax = max(x, [], 1); % maximum value per trace
            f_val = mag2db(abs(x)) - mag2db(abs(xmax)) >= kwargs.fthresh; % threshold 
            f_val = (any(f_val, 2:ndims(x))); % evaluate any freqs across aperture/frames that is above threshold
            
            % get the pixel positions
            Pi = cell(3,1);
            [Pi{:}] = self.scan.getImagingGrid();
            Pi = cellfun(@(x) {shiftdim(x, -4)}, Pi); % place I in dims 5-7
            Pi = cat(1,Pi{:}); % 3 x 1 x 1 x 1 x [I]
            
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
            w_steer = exp(-2i*pi*df.*del_tx); % transmits steering delays

            % TODO: cast data type for efficency?
            [w_tx, w_rx, w_steer, apod_tx] = dealfun(@(w) cast(w, 'like', real(x)), w_tx, w_rx, w_steer, apod_tx);
            b = zeros([1, size(Pi,2:ndims(Pi))], 'like', x);

            % beamform one frequency at a time to save memory space
            h = waitbar(0,'Beamforming ...');
            % parfor (k = 1:K)
            for k = gather(find(f_val)')
                % skip unimportant frequencies
                % if ~f_val(k), continue; end

                % report progress
                if isvalid(h), waitbar(k/K, h, char("Beamforming: " + (gather(f(k))/1e6) + " MHz")); end

                % select frequency
                xk = shiftdim(x(k,:,:,:,:,:),1); % data, in freq. domain (N x V x ...)
                % TODO: adapt for multiple frames

                % compute the greens functions on transmit/receive
                G_tx = w_tx.^(k-1); % 1 x M x 1 x 1 x [I]
                G_rx = w_rx.^(k-1); % 1 x M x 1 x 1 x [I]

                % compute the inverse steering vector on transmit
                T_tx = apod_tx .* w_steer.^(k-1); % M x V
                A_tx = pagemtimes(G_tx, T_tx); % 1 x V x 1 x 1 x [I]
                Ainv_tx = pagetranspose(A_tx); % V x 1 x 1 x 1 x [I] % make a column vector
                
                % delay and sum the data for this frequency
                yn = pagemtimes(conj(G_rx), xk); % 1 x V x 1 x 1 x [I]
                y  = pagemtimes(yn, conj(Ainv_tx)); % 1 x 1 x 1 x 1 x [I]
                
                % integrate over all frequencies
                b = b + y; % 1 x 1 x 1 x 1 x [I]

            end           
            if isvalid(h), close(h); end

            % revert to normal dimensions
            b = shiftdim(b,4); % [I] x ...
        end    
    end
    
    % Receive Aperture beamforming methods: operate along dimension 2
    methods(Static)
        function z = rcvDefault(x, dim)
            if nargin < 2, dim = 2; end
            z = UltrasoundSystem.rcvSum(x, dim);
        end
        
        function z = rcvSum(x, dim)
           if nargin < 2, dim = 2; end
           z = sum(x, dim);
        end
        
        function z = rcvSLSCAvg(x, dim, maxlag)
            %
            
            % defaults
            if nargin < 2, dim = 2; end
            if nargin < 3, maxlag = 10; end
            
            % normalize magnitude per sample (averaging)
            x = x ./ abs(x);
            x(isnan(x)) = 0; % 0/0 -> 0
            
            % get weighting filter across receiver pairs
            L = size(x,dim);
            [M, N] = ndgrid(1:L, 1:L);
            H = abs(M - N); % lag
            S = (0 < H & H <= maxlag); % valid lags for adding
            O = 1 ./ (L - H); % debiasing weights
            W = S .* O; % final weights per pair
            
            % place cross-receiver across cells
            xc = num2cell(x, setdiff(1:ndims(x), dim));
            
            % place weights as cross receiver over cells, receiver in dim
            W = num2cell(shiftdim(W, -(dim-2)), dim);
            
            % correlation sum across receivers per cross-receiver kernel 
            vn = @(w,xc) sum(w .* conj(xc) .* x, dim);
                          
            % compute product and sum per cross-receiver
            z = 0;
            parfor i = 1:numel(W)
                z = z + vn(W{i}, xc{i});
            end
        end
        
        function z = rcvSLSCEns(x, dim, maxlag)
            %

            % defaults
            if nargin < 2, dim = 2; end
            if nargin < 3, maxlag = 10; end
            
            % get weighting filter across receiver pairs
            L = size(x,dim);
            [M, N] = ndgrid(1:L, 1:L);
            H = abs(M - N);
            W = (0 < H & H <= maxlag);
            
            % place cross-receiver across cells
            xc = num2cell(x, setdiff(1:ndims(x), dim));
            
            % place weights as cross receiver over cells, receiver in "dim" 
            W = num2cell(shiftdim(W, -(dim-2)), dim);
            
            % correlation across receivers per cross-receiver kernel 
            vn = @(w,xc) deal(...
                sum(w .* conj(xc) .* x,  dim),...
                sum(w .* abs(x).^2,      dim),...
                sum(w .* abs(xc).^2 * L, dim));
                
            % compute and sum for each cross-receiver
            [z, a, b] = deal(0);
            for i = 1:numel(W)
                [zr, ar, br] = vn(W{i}, xc{i});
                z = z + zr;
                a = a + ar;
                b = b + br;
            end
            
            % get final image
            z = z ./ (a .* b);
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
    
    % helper functions
    methods
        function recompile(self), recompileMex(self); recompileCUDA(self); end
        function recompileMex(self)
            % RECOMPILEMEX(self)
            %
            %
                     
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
        function recompileCUDA(self, cuda_folder)
            % RECOMPILECUDA(self, cuda_folder)
            %
            %
            
            % CUDA folder: defaults to default install location
            if nargin < 2, cuda_folder = UltrasoundSystem.getDefaultCUDAFolder(); end
            
            % src file folder
            src_folder = UltrasoundSystem.getSrcFolder();
            
            % get all source code definitions
            defs = UltrasoundSystem.getCUDAFileDefs();
            
            % compile each
            for d = defs
                % make full command
                com = join(cat(1,...
                    fullfile(cuda_folder,'bin/nvcc'), ...
                    ['--ptx ' fullfile(src_folder, [d.Source])], ...
                    ['-o ' fullfile(self.tmp_folder, strrep(d.Source, 'cu', 'ptx'))], ...
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
    methods(Static)
        function defs = getCUDAFileDefs(cuda_folder)
            % no halp :(
            
            % CUDA folder
            if nargin < 1, cuda_folder = UltrasoundSystem.getDefaultCUDAFolder(); end

            % get all source code definitions
            defs = [...
                UltrasoundSystem.genCUDAdef_beamform(cuda_folder),...
                ];
        end

        function defs = getMexFileDefs()
            % get all source code definitions
            defs = [...
                UltrasoundSystem.genMexdef_msfm(), ... % both msfm files
                ];
        end

        function f = getSrcFolder()
            f = fileparts(mfilename('fullpath'));
        end
        
        function f = getDefaultCUDAFolder()
            %

            % TODO: search for CUDA via environmental variables
            if isunix()
                f = '/usr/local/cuda';
            else
                error("No known default CUDA folder");
            end
        end

        function d = genCUDAdef_beamform(cuda_folder)
            % no halp :(

            % CUDA folder
            if nargin < 1, cuda_folder = UltrasoundSystem.getDefaultCUDAFolder(); end

            % filename
            d.Source = 'bf.cu';

            d.IncludePath = {... include folders
                fullfile(cuda_folder, 'samples/common/inc'), ...
                };

            d.Libraries = {...
                ...
                };

            d.CompileOptions = {...  compile options
                "use_fast_math",...
                };

            d.Warnings = {... warnings
                "no-deprecated-gpu-targets"...
                };

            d.DefinedMacros = {... macro definitions
                ..."T=5", "M=3", "N=1", "I=6"...
                ...["T","M","N","I"] + "=" + string([1e3, 2e3, 4e3, 3e3])
                };
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

