% ULTRASOUNDSYSTEM - Complete ultrasound system class
%
% The ULTRASOUNDSYSTEM class is a synthesis class containing the properties
% describing a medical ultrasound system and providing methods to simulate
% channel data or beamform channel data. The complete system is described
% by the transmit and receive Transducer, the transmit Sequence, and the
% Scan defining the region for simulation or beamforming.
%
% Multiple simulators are supported, but external simulators must be
% installed separately. They include:
% 
% * greens via QUPS
% * simus via MUST
% * calc_scat_all, calc_scat_multi via FieldII
% * kspaceFirstOrder via K-wave
% * fullwaveSim via Fullwave
% 
% Multiple beamformers are provided which include
% 
% * bfDAS - a naive delay-and-sum beamformer
% * DAS - a more performant naive delay-and-sum beamformer
% * bfEikonal - a sound speed delay-and-sum beamformer using the eikonal equation
% * bfAdjoint - a frequency domain beamformer
% 
% See also CHANNELDATA TRANSDUCER SEQUENCE SCAN SCATTERERS MEDIUM 

classdef UltrasoundSystem < matlab.mixin.Copyable
    
    % objects
    properties
        tx Transducer = TransducerArray() % Transducer object (transmit)
        rx Transducer = TransducerArray() % Transducer object (receive)
        sequence Sequence = Sequence()    % Sequence object
        scan Scan = ScanCartesian() % Scan object
    end
    
    % parameters
    properties
        % FS - Simulation sampling frequency 
        %
        % ULTRASOUNDSYSTEM.FS sets the sampling frequency for simulation
        % routines.
        %
        % Example:
        %
        % us = UltrasoundSystem('xdc', TransducerArray());
        % 
        % 
        % See also GREENS SIMUS CALC_SCAT_MULTI KSPACEFIRSTORDER
        fs (1,1) {mustBeNumeric} = 40e6   % simulation sampling frequency
    end
    
    properties(Dependent)
        fc  {mustBeNumeric} % central operating frequency (from the Transducer)
        xdc Transducer      % Transducer object (if receive and transmit are identical)
        pulse Waveform      % Waveform object (from the Sequence)
    end
    
    properties(Hidden,SetAccess=protected)
        tmp_folder (1,1) string = mktempdir() % temporary folder for compiled binaries
    end
        
    % constructor/destructor
    methods
        % constructor
        function self = UltrasoundSystem(kwargs, opts)
            % ULTRASOUNDSYSTEM - Construct an UltrasoundSystem object
            %
            % us = ULTRASOUNDSYSTEM() constructs an UltrasoundSystem
            % object.
            % 
            % us = ULTRASOUNDSYSTEM(Name,Value,...) constructs an
            % UltrasoundSystem object using Name/Value pairs. 
            % 
            % us = ULTRASOUNDSYSTEM(...,'tx', tx) sets the transmit 
            % Transducer to be the Transducer tx.
            %
            % us = ULTRASOUNDSYSTEM(...,'rx', rx) sets the receive 
            % Transducer to be the Transducer rx.
            %
            % us = ULTRASOUNDSYSTEM(...,'xdc', xdc) sets the transmit and 
            % receive Transducer us.tx and us.rx to be the Transducer xdc.
            % When set like this, the property us.xdc is available.
            %
            % us = ULTRASOUNDSYSTEM(...,'sequence', seq) sets the Sequence
            % to be the Sequence seq.
            %
            % us = ULTRASOUNDSYSTEM(...,'scan', scan) sets the Scan to be
            % the Scan scan.
            %
            % us = ULTRASOUNDSYSTEM(...,'fs', fs) sets the simulation 
            % sampling frequency to fs.
            %
            % us = ULTRASOUNDSYSTEM(...,'recompile', false) avoids
            % attempting to recompile mex and CUDA files when an
            % UltrasoundSystem object is created.
            %
            % See also TRANSDUCER SEQUENCE SCAN
            arguments
                kwargs.tx Transducer
                kwargs.rx Transducer
                kwargs.xdc Transducer
                kwargs.sequence Sequence
                kwargs.scan Scan
                kwargs.fs (1,1) {mustBeNumeric}
                opts.recompile (1,1) logical = true
            end
            
            % initailize via name-Value pairs
            f = string(fieldnames(kwargs))'; % name-value fields (1 x F)

            % initialize
            for s = f, self.(s) = kwargs.(s); end
            
            % create a default Sequence if none provided
            if ~isfield(kwargs, 'sequence')
                % set the default pulse
                fc = self.tx.fc; % extract - otherwise, we reference the changing value of the object
                excitation = @(t)exp(-2j*pi*fc*t); % functional form
                P = 2; % 2 periods

                % set the default pulse sequence
                wv = Waveform( ...
                    't0', - P / 2 / self.tx.fc, ...
                    'tend', P / 2 / self.tx.fc, ...
                    'fun', excitation);

                self.sequence = Sequence(...
                    'type','FSA',...
                    'focus', [0;0;0], ...
                    'pulse', wv, ...
                    'numPulse', self.tx.numel ...
                    );
            end

            % shadow with the (newly created) temp folder for binaries and 
            % const-compiled code
            addpath(self.tmp_folder);

            % copy code or recompile it
            if gpuDeviceCount % only do CUDA stuff if there's a MATLAB-compatible GPU
                defs = self.getCUDAFileDefs();
                fls = arrayfun(@(d) string(strrep(d.Source, '.cu', '.ptx')), defs);
                e = logical(arrayfun(@exist, fullfile(self.tmp_folder, fls))); % already exists?
                s = arrayfun(@(fl) copyfile(which(fl), fullfile(self.tmp_folder, fl)), fls(~e)); % move there if not?
                if opts.recompile && any(~s), self.recompileCUDA(); end % attempt to recompile code
            end

            % copy code or recompile it
            defs = self.getMexFileDefs();
            fls = arrayfun(@(d) string(strrep(d.Source, 'c', mexext())), defs);
            e = logical(arrayfun(@exist, fullfile(self.tmp_folder, fls))); % already exists?
            s = arrayfun(@(fl) copyfile(which(fl), fullfile(self.tmp_folder, fl)), fls(~e)); % move there if not?
            if opts.recompile && any(~s), self.recompileMex(); end % attempt to recompile code

        end

        % destructor
        function delete(self)
            % DELETE - Destroy an UltrasoundSystem ... programatically.
            %
            % On object destruction, any temporary directories are removed.
            %
            % See also HANDLE
            arguments, self (1,1) UltrasoundSystem, end

            % if we made a temp folder, clean it up
            if ~isempty(self.tmp_folder) && exist(self.tmp_folder, 'dir')
                rmpath(self.tmp_folder) % remove from the path
                list = dir(self.tmp_folder); % all files in the folder
                nms = string({list.name}); % get the file names
                % check that it's only ptx and mex files we are deleting
                assert(all(...
                    endsWith(nms, [".ptx", string(mexext())]) ...
                    | nms == ".." | nms == "." ...
                    ), ...
                    "Call for deletion of " + self.tmp_folder + " failed due to unexpected files present." ...
                    ); 
                
                % rmdir(self.tmp_folder, 's'); % recursive deletion - dangerous
                
                % safe: remove specific files first, then (attempt) the folder
                rmdir(fullfile(self.tmp_folder, "*" + ".ptx")); % remove any ptx files
                rmdir(fullfile(self.tmp_folder, "*" + string(mexext()))); % remove any mex files
                rmdir(self.tmp_folder);
            end
        end
    end
    
    % overloading methods
    methods(Access=protected)
        % copy 
        function other = copyElement(self)
            arguments, self (1,1) UltrasoundSystem, end
            n = cellstr(["tx", "rx", "sequence", "scan", "fs"]); % copy props
            v = cellfun(@(n){self.(n)}, n);
            nv = cat(1,n,v);
            other = UltrasoundSystem(nv{:}, 'recompile', false);            
        end
    end
    methods
        function sub_div = getLambdaSubDiv(self, p, c, ap)
            % GETLAMBDASUBDIV - Get subelement divisions w.r.t. wavelength
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

    % display
    methods
        function h = plot(self, ax, im_args)
            % PLOT - Plot the geometry of the UltrasoundSystem 
            %
            % h = PLOT(self) plots the UltrasoundSystem self by plotting 
            % the transmit and receive transducer(s), the imaging pixels,
            % and representation of the transmit sequence on the same axes.
            %
            % h = PLOT(self, ax) plots on the axes ax. The default is the
            % current axes returned by gca.
            %
            % h = PLOT(..., Name, Value, ...) passes following arguments to
            % the call to plot for each of the objects i.e. to the plot
            % function for self.scan, self.sequence, and self.tx (and
            % self.rx if necessary).
            %
            % Example:
            % % Create a default system using a focused transmit
            % xdc = TransducerArray();
            % us = UltrasoundSystem('xdc', xdc);
            % us.sequence.type = 'VS';
            % us.sequence.focus = [0;0;30e-3] + ...
            % linspace(-xdc.numel/4, xdc.numel/4, xdc.numel/2+1) .* [xdc.pitch;0;0];
            %
            % % plot it
            % figure;
            % h = plot(us, 'LineWidth', 2);
            %
            % See also PLOT ULTRASOUNDSYSTEM/XDC ULTRASOUNDSYSTEM/SCAN ULTRASOUNDSYSTEM/SEQUENCE

            arguments
                self UltrasoundSystem
                ax (1,1) matlab.graphics.axis.Axes = gca;
            end
            arguments
                im_args.?matlab.graphics.primitive.Line
            end

            hold(ax, 'on');
            set(ax, 'ydir', 'reverse');
            title(ax, 'Geometry');
            plargs = struct2nvpair(im_args);
            hps = plot(self.scan, 'm.', 'DisplayName', 'Image', plargs{:}); % the imaging points
            if self.tx == self.rx
                hxdc = plot(self.xdc, ax, 'r+', 'DisplayName', 'Elements', plargs{:}); % elements
            else
                htx = plot(self.tx, ax, 'b+', 'DisplayName', 'Transmit Elements', plargs{:}); % tx elements
                hrx = plot(self.rx, ax, 'r+', 'DisplayName', 'Receive Elements', plargs{:}); % rx elements
                hxdc = [htx, hrx];
            end
            switch self.sequence.type % show the transmit sequence
                case 'PW', hseq = plot(self.sequence, ax, 3e-2, 'k.', 'DisplayName', 'Tx Sequence', plargs{:}); % scale the vectors for the plot
                otherwise, hseq = plot(self.sequence, ax, 'k.', 'DisplayName', 'Tx Sequence', plargs{:}); % plot focal points, if they exist
            end
            h = [hxdc, hps, hseq];
            legend(ax, h);
            grid(ax, 'on'); 
            grid(ax, 'minor')
        end
    end

    % property modification
    methods
        function self = scale(self, kwargs)
            % SCALE - Scale the units of the system
            %
            % self = SCALE(self, 'dist', dscale) scales the values of 
            % distance by dscale.
            % 
            % self = SCALE(self, 'time', tscale) scales the values of
            % time by tscale, and values of (temporal) frequency by the
            % inverse of tscale.
            %
            % Example:
            % % Define a system in meters, seconds
            % xdc = TransducerArray('pitch', 1e-3, 'fc', 1e6); % defined in meters, seconds
            % us = UltrasoundSystem('xdc', xdc); 
            %
            % % Scale the values to millimiters, microseconds
            % us = scale(us, 'dist', 1e3, 'time', 1e6);
            % 
            %
            arguments
                self (1,1) UltrasoundSystem
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end

            self = copy(self);
            args = struct2nvpair(kwargs);
            
            % Scan: convert distance only
            if isfield(kwargs, 'dist')
                self.scan = scale(self.scan, 'dist', kwargs.dist);
            end
            % Sampling: convert time only
            if isfield(kwargs, 'time')
                self.fs = self.fs / kwargs.time;
            end

            % Transducer: convert distance and/or time/freq
            if self.tx == self.rx,
                self.xdc = scale(self.xdc, args{:}); % convert in distance and time
            else
                self.tx = scale(self.tx, args{:}); % convert in distance and time
                self.rx = scale(self.rx, args{:}); % convert in distance and time
            end

            % Sequence
            self.sequence = scale(self.sequence, args{:}); % convert in distance and time
        end
    end

    % Modified Green's function based direct computations
    methods
        function [chd, wv] = greens(self, scat, element_subdivisions, kwargs)
            % GREENS - Simulate ChannelData via a shifted Green's function.
            % 
            % chd = GREENS(self, scat) simulates the response of the 
            % Scatterers scat from the UltrasoundSystem self and returns 
            % the corresponding ChannelData chd. It first computes the full
            % synthetic aperture data using a simple Green's function 
            % kernel applied to all sub-elements and all point scatterers 
            % and then applies the transmit Sequence via focusTx.
            %
            % chd = GREENS(self, scat, element_subdivisions) uses the 
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
            % [...] = GREENS(..., 'interp', method) specifies the method for
            % interpolation. Support is provided by the ChannelData/sample 
            % method. The default is 'cubic'.
            % 
            % [...] = GREENS(..., 'device', index) uses the gpuDevice with ID
            % index. If index == -1, the current device is used. Index == 0
            % avoids using a gpuDevice. The default is -1 if a gpu is
            % available.
            %
            % [...] = GREENS(..., 'bsize', B) uses a block size of B when
            % vectorizing computations. A larger block size will run
            % faster, but use more memory. The default is 1.
            %
            % [...] = GREENS(..., 'device', 0, 'tall', true, ...) uses a 
            % tall type for intermediate computations. This may help 
            % prevent out-of-memory errors.
            %
            % Example:
            % 
            % % Simulate some data
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % chd = greens(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            %
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_MULTI ULTRASOUNDSYSTEM/FOCUSTX CHANNELDATA/SAMPLE
            
            arguments
                self (1,1) UltrasoundSystem
                scat Scatterers
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
            end
            arguments
                kwargs.device (1,1) {mustBeInteger} = -logical(gpuDeviceCount());
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.tall (1,1) logical = false;
                kwargs.bsize (1,1) {mustBeInteger, mustBePositive} = 1;
                kwargs.verbose (1,1) logical = false;
            end
            
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
            taumax = arrayfun(@(scat)(2 * maxdist(scat.pos) + maxdist(ptc_tx) + maxdist(ptc_rx)) ./ min(scat.c0), shiftdim(scat(:),-3));
            taumin = arrayfun(@(scat)(2 * mindist(scat.pos) + mindist(ptc_tx) + mindist(ptc_rx)) ./ max(scat.c0), shiftdim(scat(:),-3));
            
            % Directly convolve the Waveform objects to get the final
            % convolved kernel
            wv = conv(self.rx.impulse, ...
                conv(self.tx.impulse, self.sequence.pulse, self.fs), ...
                self.fs); % transmit waveform, convolved at US frequency
            wv.fs = self.fs;
            kern = wv.samples;

            F = numel(scat);
            if kwargs.verbose, hw = waitbar(0); end

            for f = F:-1:1 % for each Scatterers
            % get minimum/maximum sample times
            tmin = taumin(f) + wv.t0;
            tmax = taumax(f) + wv.tend;

            % create time vector (T x 1)
            % this formulation is guaranteed to pass through t == 0
            n0 = floor(tmin * self.fs);
            ne = ceil(tmax * self.fs);
            t = (n0 : ne)';

            % pre-allocate output
            [T, N, M, E] = deal(numel(t), self.rx.numel, self.tx.numel, prod(element_subdivisions));
            x   = complex(zeros([1 T N M]));

            % splice
            c0  = scat(f).c0;
            pos = scat(f).pos; % 3 x S
            amp = scat(f).amp; % 1 x S
            fs_ = self.fs;
            if kwargs.device && exist('greens.ptx', 'file') ... % use the GPU kernel
                    && (ismember(kwargs.interp, ["nearest", "linear", "cubic", "lanczos3"]))
                % function to determine type
                isftype = @(x,T) strcmp(class(x), T) || any(arrayfun(@(c)isa(x,c),["tall", "gpuArray"])) && strcmp(classUnderlying(x), T);

                % determine the data type
                if     isftype(kern, 'double'), suffix = "" ; cfun = @double;
                elseif isftype(kern, 'single'), suffix = "f"; cfun = @single;
                elseif isftype(kern, 'halfT' ), suffix = "h"; cfun = @(x) alias(halfT(x));
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
                    x, pos, amp, ptc_rx, ptc_tx, kern, t(1)/fs_, wv.t0, fs_, 1/c0 ...
                    );

                % re-map sizing
                [QI, QS, QT, QN, QM] = deal(scat(f).numScat, length(t), length(kern), N, M);

                % compute the minimum and maximum time delay for each scatterer
                ps = gpuArray(ps);
                [rminrx, rmaxrx, rmintx, rmaxtx] = deal(+inf, -inf, +inf, -inf);
                for p = self.tx.positions, rminrx = min(rminrx, vecnorm(ps - p, 2, 1) ./ c0); end % minimum scat time
                for p = self.tx.positions, rmaxrx = max(rmaxrx, vecnorm(ps - p, 2, 1) ./ c0); end % maximum scat time
                for p = self.rx.positions, rmintx = min(rmintx, vecnorm(ps - p, 2, 1) ./ c0); end % minimum scat time
                for p = self.rx.positions, rmaxtx = max(rmaxtx, vecnorm(ps - p, 2, 1) ./ c0); end % maximum scat time
                [rmin, rmax] = deal(rmintx + rminrx, rmaxtx + rmaxrx);

                % sort points by their maximum delay
                [~, i] = sort(rmax); % get sorting
                [ps, as, rmin, rmax] = dealfun(@(x)sub(x,i,2), ps,as,rmin,rmax); % apply to all position variables

                % get the index bounds for the output time axis
                sb = ([rmin; rmax] + t0x - t0k) * fs_ + [0; QT];

                % grab the kernel reference
                k = parallel.gpu.CUDAKernel('greens.ptx', 'greens.cu', 'greens' + suffix);
                k.setConstantMemory( 'QUPS_S', uint64(QS) ); % always set S
                try k.setConstantMemory('QUPS_T', uint64(QT), ...
                    'QUPS_N', uint64(QN), 'QUPS_M', uint64(QM) ... 
                    ); end % already set by const compiler
                try k.setConstantMemory('QUPS_I', uint64(QI)); end % might be already set by const compiler
                k.ThreadBlockSize = min(k.MaxThreadsPerBlock,QS); 
                k.GridSize = [ceil(QS ./ k.ThreadBlockSize(1)), N, M];

                % call the kernel
                x = k.feval(x, ps, as, pn, pv, kn, sb, [t0k, t0x, fs_, cinv_], [E,E], flagnum);
                
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
                svec = num2cell((1:kwargs.bsize)' + (0:kwargs.bsize:scat(f).numScat-1), 1);
                svec{end} = svec{end}(svec{end} <= scat(f).numScat);
                S = numel(svec); % number of scatterer blocks

                for sv = 1:S, s = svec{sv}; % vector of indices
                    for em = 1:E
                        for en = 1:E
                            % compute time delays
                            % TODO: do this via ray-path propagation through a
                            % medium
                            % S x 1 x N x M x 1 x 1
                            r_rx = vecnorm(sub(pos.',s,1) - sub(ptc_rx, en, 5),2,2);
                            r_tx = vecnorm(sub(pos.',s,1) - sub(ptc_tx, em, 6),2,2);
                            tau_rx = (r_rx ./ c0); % S x 1 x N x 1 x 1 x 1
                            tau_tx = (r_tx ./ c0); % S x 1 x 1 x M x 1 x 1

                            % compute the attenuation (S x 1 x [1|N] x [1|M] x 1 x 1)
                            att = sub(amp.',s,1);% .* (1 ./ r_rx) .* (1 ./ r_tx); % propagation attenuation

                            % get 0-based sample time delay
                            % switch time and scatterer dimension
                            % S x T x N x M x 1 x 1
                            t0 = gather(wv.t0); % 1 x 1

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
                                it = it & T-1 >= tvec - max(tau(:));
                                it = it &   0 <= tvec - min(tau(:));

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
            x = reshape(x, size(x,2:ndims(x))); % same as shiftdim(x,1), but without memory copies on GPU
            chd(f) = ChannelData('t0', sub(t,1,1) ./ fs_, 'fs', fs_, 'data', x);

            % truncate the data if possible
            iszero = all(chd(f).data == 0, 2:ndims(chd(f).data)); % true if 0 for all tx/rx/targs
            n0 = gather(find(cumsum(~iszero, 'forward'), 1, 'first'));
            T_ = gather(find(cumsum(~iszero, 'reverse'), 1, 'last' ));
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
        function conf = fullwaveConf(self, medium, sscan, kwargs)
            % FULLWAVECONF - Generate a Fullwave simulation configuration
            %
            % conf = FULLWAVECONF(self, medium) creates a simulation
            % configuration struct to be used with fullwaveJob to simulate 
            % the response from the Medium medium using the 
            % UltrasoundSystem self.
            % 
            % conf = FULLWAVECONF(self, medium, sscan) uses the  
            % ScanCartesian sscan as the simulation region. The default is
            % self.scan.
            %
            % conf = FULLWAVECONF(..., 'f0', f0) uses a reference frequency
            % of f0 to configure the simulation. The default is self.tx.fc.
            %
            % conf = FULLWAVECONF(..., 'CFL_max', cfl) scales the 
            % sampling frequency until the CFL is at most cfl. A rule of
            % thumb for a stable simulation is 0.3. A lesser CFL will be
            % more numerically stable whereas a greater CFL will be faster
            % to compute, but will be less numerically accurate and may
            % lead to instability.
            % 
            % conf = FULLWAVECONF(..., 'txdel', method) uses the specified 
            % method for applying the transmit delays. The 'discrete'
            % method applies delays per pixel rounded to the nearest time
            % interval. The 'continuous' method samples the transmit
            % waveform continuously at the exact delays, but takes longer
            % to compute. The 'interpolate' method samples the transmit
            % waveform once and then resamples the waveform with the 
            % transmit delays. This method offers more compute performance
            % while sacrificing much less accuracy.
            % 
            % Example:
            % 
            % % Setup a system
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(sscan.size), 1000*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            % rho(Xg == 0 & Zg == 30e-3) = 1000*2; % double the density
            % med = Medium.Sampled(sscan, c, rho);
            % 
            % % Simulate the ChannelData
            % simdir = fullfile(pwd, 'fwsim'); % simulation directory
            % conf = fullwaveConf(us, med, sscan, 'CFL_max', 0.5); % configure the sim
            % job = UltrasoundSystem.fullwaveJob(conf, 'simdir', simdir); % create a job
            % submit(job); % submit the job
            % ... do other things
            % wait(job); % wait for the job to finish processing
            % chd = UltrasoundSystem.readFullwaveSim(simdir); % extract the ChannelData
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(real(chd));
            % 
            % See also ULTRASOUNDSYSTEM/FULLWAVEJOB

            arguments
                self (1,1) UltrasoundSystem
                medium Medium
                sscan ScanCartesian = self.scan
                kwargs.f0 (1,1) {mustBeNumeric} = self.tx.fc; % center frequency of the transmit / simulation
                kwargs.CFL_max (1,1) {mustBeReal, mustBePositive} = 0.3 % maximum CFL
                kwargs.txdel (1,1) string {mustBeMember(kwargs.txdel, ["discrete", "continuous", "interpolate"])} = 'interpolate';
            end            

            %% Configuration variables

            % basic vars
            c0       = medium.c0;       % speed of sound (m/s)
            omega0   = 2*pi*kwargs.f0;  % center radian frequency of transmitted wave
            dur      = diff(sscan.zb)*2.3/c0; % duration of simulation (s) TODO: make this more general, or an input?

            % determine other grid vars
            dX   = min(abs([sscan.dx, sscan.dy, sscan.dz]), [], 'omitnan');  % limit of spatial step size - will be the same in both dimensions
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
            xdcfw = self.xdc.getFullwaveTransducer(sscan);

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
            wv_tx.dt = dT; % set the sampling interval
            t = wv_tx.time; % get discrete sampling times
            nTic = numel(t); % number of transmit samples
            
            % get transmit signal per input-pixel, time index, transmit for nTic time
            % indices 
            tau_tx_pix = shiftdim(tau_tx_pix,  -1); % 1 x nPxIn x nTx, nTx == M
            switch kwargs.txdel
                % apply discrete shifts in time
                case 'discrete',  icmat = cell2mat(arrayfun(@(tau) circshift(wv_tx.sample(t(:)), round(tau ./ dT)), tau_tx_pix, 'UniformOutput', false));
                    % continuous resampling of the waveform
                case 'continuous', icmat = cell2mat(arrayfun(@(tau) {wv_tx.sample(t(:) - tau)}, tau_tx_pix));
                case 'interpolate' % interpolate, upsampling by 10x in time first
                    t_up = wv_tx.getSampleTimes(10*fs_);
                    icmat = interp1(t_up, wv_tx.sample(t_up), t(:) - tau_tx_pix, 'spline', 0); 
            end

            % apply transmit apodization
            icmat = icmat .* shiftdim(tx_apod,-1); % (T' x nInPx x nTx)

            %% Define the Medium
            % get maps in some order
            maps = medium.getFullwaveMap(sscan);
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
        
        function [chd, conf] = fullwaveSim(self, medium, sscan, kwargs)
            % FULLWAVESIM - Simulate channel data via Fullwave
            %
            % chd = FULLWAVESIM(self, medium, sscan) simulates the Medium 
            % medium on the simulation grid sscan and returns a ChannelData
            % object chd. The simulation scan should be large and fine 
            % enough that all elements of the Transducer can be placed.
            %
            % chd = FULLWAVESIM(self, medium) uses self.scan as the
            % simulation scan sscan.
            %
            % chd = FULLWAESIM(..., 'simdir', dir) uses the directory dir
            % to store temporary simulation files. The default is a folder
            % in the working directory.
            %
            % chd = FULLWAESIM(..., 'parcluster', clu) uses the
            % parallel.Cluster clu to compute each pulse. The simulation
            % directory must be accesible by the parallel.Cluster clu.
            % 
            % [chd, conf] = FULLWAVESIM(...) also returns the configuration
            % structure used to launch the simulation.
            %
            % Example:
            % 
            % % Setup a system
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(sscan.size), 1000*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            % rho(Xg == 0 & Zg == 30e-3) = 1000*2; % double the density
            % med = Medium.Sampled(sscan, c, rho);
            % 
            % % Simulate the ChannelData
            % chd = fullwaveSim(us, med, sscan);
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(real(chd));
            % 
            % See also ULTRASOUNDSYSTEM/KSPACEFIRSTORDERND
            
            arguments
                self (1,1) UltrasoundSystem
                medium Medium
                sscan ScanCartesian = self.scan
                kwargs.parcluster (1,1) parallel.Cluster = parcluster() % parallel cluster
                kwargs.simdir (1,1) string = fullfile(pwd, 'fwsim'); % simulation directory
                kwargs.f0 (1,1) {mustBeNumeric} = self.xdc.fc; % center frequency of the transmit / simulation
                kwargs.CFL_max (1,1) {mustBeReal, mustBePositive} = 0.5 % maximum CFL
                kwargs.txdel (1,1) string {mustBeMember(kwargs.txdel, ["disc", "cont", "terp"])} = 'terp'; % delay method
            end            

            % create the configuration
            conf_args = rmfield(kwargs, setdiff(fieldnames(kwargs), {'f0', 'CFL_max', 'txdel'}));
            conf_args_ = struct2nvpair(conf_args);
            conf = fullwaveConf(self, medium, sscan, conf_args_{:});

            % create a job to process it
            job = UltrasoundSystem.fullwaveJob(conf, kwargs.parcluster, 'simdir', kwargs.simdir);

            % submit the job
            submit(job);

            % wait for it to finish
            wait(job);

            % read in the data
            chd = UltrasoundSystem.readFullwaveSim(kwargs.simdir);
        end
    end
    methods(Static)
        function job = fullwaveJob(conf, clu, kwargs)
            % FULLWAVEJOB - Create a Fullwave simulation job.
            %
            % job = ULTRASOUNDSYSTEM.FULLWAVEJOB(conf) creates a parallel.Job 
            % to run the fullwave simulation from conf.
            %
            % job = ULTRASOUNDSYSTEM.FULLWAVEJOB(conf, clu) creates a 
            % parallel.Job on the parallel.Cluster clu. The parallel.Cluster clu 
            % must have access to the simulation directory simdir. The resource 
            % configurations for the parallel.Cluster clu should be setup prior
            % to this call to ensure enough RAM is present.
            %
            % job = ULTRASOUNDSYSTEM.FULLWAVEJOB(..., 'simdir', dir) 
            % uses the directory dir to store the simulation inputs and outputs.
            %
            % Example:
            % 
            % % Setup a system
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(sscan.size), 1000*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            % rho(Xg == 0 & Zg == 30e-3) = 1000*2; % double the density
            % med = Medium.Sampled(sscan, c, rho);
            % 
            % % Simulate the ChannelData
            % simdir = fullfile(pwd, 'fwsim'); % simulation directory
            % conf = fullwaveConf(us, med, sscan, 'CFL_max', 0.5); % configure the sim
            % job = UltrasoundSystem.fullwaveJob(conf, 'simdir', simdir); % create a job
            % submit(job); % submit the job
            % ... do other things
            % wait(job); % wait for the job to finish processing
            % chd = UltrasoundSystem.readFullwaveSim(simdir); % extract the ChannelData
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(real(chd));
            % 
            % See also FULLWAVECONF PARCLUSTER PARALLEL.CLUSTER

            arguments
                conf (1,1) struct
                clu (1,1) parallel.Cluster = parcluster()
                kwargs.simdir (1,1) string = fullfile(pwd, 'fwsim')
            end

            % make simulation directory
            simdir = kwargs.simdir;
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
            job.UserData = struct('simdir', kwargs.simdir);

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
            % Example:
            % 
            % % Setup a system
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(sscan.size), 1000*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            % rho(Xg == 0 & Zg == 30e-3) = 1000*2; % double the density
            % med = Medium.Sampled(sscan, c, rho);
            % 
            % % Simulate the ChannelData
            % simdir = fullfile(pwd, 'fwsim'); % simulation directory
            % conf = fullwaveConf(us, med, sscan, 'CFL_max', 0.5); % configure the sim
            % job = UltrasoundSystem.fullwaveJob(conf, 'simdir', simdir); % create a job
            % submit(job); % submit the job
            % ... do other things
            % wait(job); % wait for the job to finish processing
            % chd = UltrasoundSystem.readFullwaveSim(simdir); % extract the ChannelData
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(real(chd));
            % 
            % See also RUNFULLWAVETX ULTRASOUNDSYSTEM/FULLWAVEJOB
            arguments
                simdir (1,1) string
                conf (1,1) struct
            end

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
        function chd = simus(self, scat, kwargs, simus_kwargs)
            % SIMUS - Simulate channel data via MUST
            %
            % chd = SIMUS(self, scat) simulates the Scatterers scat and
            % returns a ChannelData object chd.
            %
            % When calling this function the transmit sequence pulse is 
            % ignored. SIMUS only supports tone bursts at the central 
            % frequency of the transducer.
            %
            % The transmit and receive transducer must be identical i.e. 
            % self.rx == self.tx must be true.
            %
            % chd = SIMUS(...,'dims', D, ...) selects the number of 
            % dimensions for the simulation. D must be one of {2, 3} or 
            % empty. If D is empty, the dimensions are chosen based on the 
            % point scatterers.
            %
            % chd = SIMUS(...,'periods', T, ...) selects the number of
            % periods of the tone burst. T can be any positive number.
            %
            % chd = SIMUS(...,'interp', method, ...) selects the method of
            % interpolation to use when synthesizing transmits from the
            % full-synthetic-aperture data.
            %
            % chd = SIMUS(...,'parenv', clu, ...) or 
            % chd = SIMUS(...,'parenv', pool, ...) uses the parallel.Cluster 
            % clu or the parallel.Pool pool to run each transmit in 
            % parallel. 
            % 
            % chd = SIMUS(..., 'parenv', 0) avoids using a parallel.Cluster
            % or parallel.Pool. Use 0 when operating on a GPU or if memory 
            % usage explodes on a parallel.ProcessPool.
            %
            % The default is the current parpool returned by gcp. To use a
            % parallel.ThreadPool with SIMUS, see this <a href="matlab:web('https://github.com/thorstone25/qups/issues/2')">issue</a>.
            % 
            % Example:
            % 
            % % Simulate some data
            % us = UltrasoundSystem(); % get a default system
            % us.rx = us.tx; % ensure the receiver and transmitter are identical
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % chd = simus(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            % 
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_ALL FOCUSTX
            
            arguments
                self (1,1) UltrasoundSystem
                scat (1,1) Scatterers
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate')
                simus_kwargs.periods (1,1) {mustBePositive} = 1
                simus_kwargs.dims {mustBeScalarOrEmpty, mustBeMember(simus_kwargs.dims, [2,3])} = []
            end

            % load options
            if isempty(kwargs.parenv), kwargs.parenv = 0; end % select 0 workers if empty

            % TODO: check the transmit/receive/sequence impulse: they 
            % cannot be satisfied if not a Delta or empty
            %if ~ismember("periods", varargin(cellfun(@ischar,varargin) | cellfun(@isstring, varargin))) % was this an input?
                %warning("QUPS:UltrasoundSystem:simus:unsatisfiable", "Transmit sequence determined by 'periods', property.");
            %end
            
            % get the points and the dimensions of the simulation(s)
            [X, Y, Z, A] = arrayfun(@(scat) ...
                deal(sub(scat.pos,1,1), sub(scat.pos,2,1), sub(scat.pos,3,1), scat.amp), ...
                scat, 'UniformOutput',false);
            if isempty(simus_kwargs.dims) 
                if all(cellfun(@(Y)all(Y == 0,'all'),Y), 'all'), simus_kwargs.dims = 2; [Y{:}] = deal([]); % don't simulate in Y if it is all zeros 
                else, kwargs.dims = 3; end
            end
            if simus_kwargs.dims == 2 && cellfun(@(Y)any(Y ~= 0, 'all'),Y)
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
            p.TXnow = simus_kwargs.periods; % number of wavelengths
            p.TXapodization = zeros([self.xdc.numel,1]); % set tx apodization
            p.RXdelay = zeros([self.xdc.numel,1]); % receive delays (none)
            
            % set options per Scatterers
            p = repmat(p, [1,1,1,numel(scat)]); % (1 x 1 x 1 x F)
            pxdc = arrayfun(@(scat) {getSIMUSParam(scat)}, scat); % properties per Scatterers
            for f = 1:numel(scat), 
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

            % select the computing environment
            parenv = kwargs.parenv;
            if isempty(parenv), parenv = 0; end % empty pool -> 0
            isloc = ~isa(parenv, 'parallel.Pool') || ~isa(parenv, 'parallel.Cluster'); % local or parpool
            if isloc, [pclu, parenv] = deal(parenv, 0); else, pclu = 0; end % cluster or local

            % call the sim: FSA approach
            [M, F] = deal(self.xdc.numel, numel(scat)); % splice
            for f = F:-1:1 % per scat
                argf = {X{f},Y{f},Z{f},A{f},zeros([M,1]),p(f),opt}; % args per scat
                parfor (m = 1:M, pclu) % use parallel rules, but execute on main thread
                    args = argf; % copy settings for this frame
                    args{6}.TXapodization(m) = 1; % transmit only on element m
                    if isloc, rf{m,f} = simus(args{:}); % local compute
                    else, out(m,f) = parfeval(parenv, @simus, 1, args{:}); % add cluster job
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
    methods
        function chd = calc_scat_all(self, scat, element_subdivisions, kwargs)
            % CALC_SCAT_ALL - Simulate channel data via FieldII
            %
            % chd = CALC_SCAT_ALL(self, scat) simulates the Scatterers
            % scat and returns a ChannelData object chd.
            % 
            % chd = CALC_SCAT_ALL(self, scat, element_subdivisions)
            % specifies the number of subdivisions in width and height for 
            % each element. The default is [1,1].
            %
            % chd = CALC_SCAT_ALL(..., 'interp', method) specifies the
            % interpolation methods for the transmit synthesis. The method
            % must be supported by focusTx.
            %
            % chd = CALC_SCAT_ALL(..., 'parenv', clu) or 
            % chd = CALC_SCAT_ALL(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to
            % parallelize computations. parallel.ThreadPools are ignored
            % due to mex function restrictions.
            % 
            % chd = CALC_SCAT_ALL(..., 'parenv', 0) avoids using a 
            % parallel.Cluster or parallel.Pool. Use 0 when operating on a 
            % GPU or if memory usage explodes on a parallel.ProcessPool.
            %
            % The default is the current pool returned by gcp.
            % 
            % Example:
            % 
            % % Simulate some data
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % chd = calc_scat_all(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            % 
            % See also ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/CALC_SCAT_MULTI FOCUSTX

            arguments
                self (1,1) UltrasoundSystem
                scat Scatterers
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool"])} = gcp('nocreate')
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
            end
            
            % helper function
            vec = @(x) x(:); % column-vector helper function

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = copy(self.tx.impulse); % transmitter impulse
            wv_rx = copy(self.rx.impulse); % receiver impulse
            wv_pl = copy(self.sequence.pulse);
            
            % get the time axis (which passes through t == 0)
            [wv_tx.fs, wv_rx.fs, wv_pl.fs] = deal(self.fs);
            t_tx = wv_tx.time;
            t_rx = wv_rx.time;
            t_pl = wv_pl.time;

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.samples)')));
            rx_imp = gather(double(real(vec(wv_rx.samples)')));
            tx_pls = gather(double(real(vec(wv_pl.samples)')));

            % choose the cluster to operate on: avoid running on ThreadPools
            parenv = kwargs.parenv;
            if isempty(parenv) || isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool'), parenv = 0; end

            % splice
            [M, F] = deal(self.sequence.numPulse, numel(scat)); % number of transmits/frames
            [fs_, tx_, rx_] = deal(self.fs, self.tx, self.rx); % splice
            [c0, pos, amp] = arrayfun(@(t)deal(t.c0, {t.pos}, {t.amp}), scat); % splice

            % Make position/amplitude and transducers constants across the workers
            if isa(parenv, 'parallel.Pool'), 
                cfun = @parallel.pool.Constant;
            else, 
                cfun = @(x)struct('Value', x);
                [pos, amp] = deal({pos},{amp}); % for struct to work on cell arrays
            end
            [pos_, amp_, tx_, rx_] = dealfun(cfun, pos, amp, tx_, rx_);

            % for each scat (frame)
            parfor (f = 1:F, parenv) % each transmit pulse
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

        function chd = calc_scat_multi(self, scat, element_subdivisions, kwargs)
            % CALC_SCAT_MULTI - Simulate channel data via FieldII
            %
            % chd = CALC_SCAT_MULTI(self, scat) simulates the Scatterers 
            % scat and returns a ChannelData object chd.
            %
            % chd = CALC_SCAT_MULTI(self, scat, element_subdivisions)
            % specifies the number of subdivisions in width and height for
            % each element. The default is [1, 1].
            %
            % chd = CALC_SCAT_MULTI(..., 'parenv', clu) or 
            % chd = CALC_SCAT_MULTI(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to
            % parallelize computations. parallel.ThreadPools are invalid
            % due to mex function restrictions.
            % 
            % chd = CALC_SCAT_MULTI(..., 'parenv', 0) avoids using a 
            % parallel.Cluster or parallel.Pool. Use 0 when operating on a 
            % GPU or if memory usage explodes on a parallel.ProcessPool.
            %
            % The default is the current pool returned by gcp.
            % 
            % Example:
            % 
            % % Simulate some data
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % chd = calc_scat_multi(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            % 
            % See also ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/FOCUSTX

            arguments
                self (1,1) UltrasoundSystem
                scat Scatterers
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool"])} = gcp('nocreate')
            end
            
            % helper function
            vec = @(x) x(:); % column-vector helper function

            % initialize field II
            try
                evalc('field_info');
                field_started = false;
            catch
                field_init(-1);
                field_started = true;
            end

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = copy(self.tx.impulse); % transmitter impulse
            wv_rx = copy(self.rx.impulse); % receiver impulse
            wv_pl = copy(self.sequence.pulse);

            % get the time axis (which passes through t == 0)
            [wv_tx.fs, wv_rx.fs, wv_pl.fs] = deal(self.fs);
            t_tx = wv_tx.time;
            t_rx = wv_rx.time;
            t_pl = wv_pl.time;

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.samples)')));
            rx_imp = gather(double(real(vec(wv_rx.samples)')));
            tx_pls = gather(double(real(vec(wv_pl.samples)')));

            % get the apodization and time delays across the aperture
            apod_tx = self.sequence.apodization(self.tx); % N x M
            tau_tx = -self.sequence.delays(self.tx); % N x M
            tau_offset = min(tau_tx, [], 1); % (1 x M)
            tau_tx = tau_tx - tau_offset; % 0-base the delays for FieldII

            % choose the parallel environment to operate on: avoid running on ThreadPools
            parenv = kwargs.parenv;
            if isempty(parenv) || isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool'), parenv = 0; end

            [M, F] = deal(self.sequence.numPulse, numel(scat)); % number of transmits/frames
            [fs_, tx_, rx_] = deal(self.fs, self.tx, self.rx); % splice
            [c0, pos, amp] = arrayfun(@(t)deal(t.c0, {t.pos}, {t.amp}), scat); % splice
            
            % Make position/amplitude and transducers constants across the workers
            if isa(parenv, 'parallel.Pool'), 
                cfun = @parallel.pool.Constant;
            else, 
                cfun = @(x)struct('Value', x);
                [pos, amp] = deal({pos},{amp}); % for struct to work on cell arrays
            end
            [pos_, amp_, tx_, rx_] = dealfun(cfun, pos, amp, tx_, rx_);

            parfor (m = 1:M, parenv) % each transmit pulse
            for (f = F:-1:1) % each scat frame            
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
                [voltages{m,f}, ts{m,f}] = calc_scat_multi(Tx, Rx, pos_.Value{f}.', amp_.Value{f}.'); %#ok<PFBNS> % constant over workers
            end
            end
            
            % adjust start time based on signal time definitions
            t0 = ... cell2mat(ts) + ... % fieldII start time (1 x 1 x M x F)
                (t_pl(1) + t_tx(1) + t_rx(1)) ... signal delays for impulse/excitation
                + shiftdim(tau_offset,-1) ... 0-basing the delays across the aperture
                ; % 1 x 1 x M
            
            % create the output QUPS ChannelData object 
            chd = cellfun(@(x, t0) ChannelData('data', x, 't0', t0), voltages, ts); % per transmit/frame (M x F) object array
            chd = arrayfun(@(f) join(chd(:,f), 3), 1:F); % join over transmits (1 x F) object array
            chd = join(chd, 4); % join over frames (1 x 1) object array

            % set sampling frequency and transmit times for all
            chd.fs = self.fs;
            chd.t0 = chd.t0 + t0;

            % cleanup
            if field_started, evalc('field_end'); end
        end
    end
    
    % k-Wave calls
    methods
        function [chd, readfun] = kspaceFirstOrder(self, med, sscan, kwargs, karray_args, kwave_args)
            % KSPACEFIRSTORDER - Simulate channel data via k-Wave
            % 
            % chd = KSPACEFIRSTORDER(self, med) simulates the Medium med 
            % and returns a ChannelData object chd via k-Wave.
            %
            % chd = KSPACEFIRSTORDER(self, med, sscan) operates using
            % the simulation region defined by the ScanCartesian sscan. The
            % default is self.scan.
            % 
            % If using an element mapping method provided by kWaveArray,
            % the step sizes in all dimensions must be identical.
            %
            % chd = KSPACEFIRSTORDER(..., 'T', T) runs the simulation
            % until time T. The default is computed heuristically to
            % include a two-way propagation of a signal at the slowest
            % sound speed.
            % 
            % chd = KSPACEFIRSTORDER(..., 'PML', P) or 
            % chd = KSPACEFIRSTORDER(..., 'PML', [PL, PU]) uses a PML of 
            % size P or with a size between PL and PU with the smallest 
            % maximum prime factor. The default is [20 56].
            % 
            % chd = KSPACEFIRSTORDER(..., 'CFL_max', cfl) scales the
            % sampling frequency until the CFL is at most cfl. A rule of
            % thumb for a stable simulation is 0.3. A lesser CFL will be
            % more numerically stable whereas a greater CFL will be faster
            % to compute, but will be less numerically accurate and may
            % lead to instability.
            % 
            % chd = KSPACEFIRSTORDER(..., 'ElemMapMethod', method) 
            % specifies the computational method used for mapping 
            % transducer elements to the computational grid. Must be one of 
            % {'nearest*, 'linear', 'karray-direct', 'karray-depend'}. 
            % 
            % The 'nearest' method uses the nearest pixel. The 'linear' 
            % method uses linear interpolation weights. The 'karray-direct'
            % method uses the karray method but avoids recomputing 
            % intermediate results. The 'karray-depend' method always uses 
            % the kWaveArray methods, but can be slower.
            % 
            % chd = KSPACEFIRSTORDER(..., 'parenv', clu) or 
            % chd = KSPACEFIRSTORDER(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to compute
            % each pulse in parallel via parfor. 
            % 
            % When using a parallel.Cluster clu, the settings for clu must
            % be setup for an independent parallel.Job as separate job is
            % created for each pulse. For example, if using a SLURM cluster
            % with multiple GPUs, clu.SubmitArguments should contain 
            % ' --gpus=1' so that 1 GPU is requested per pulse.
            % 
            % [job, readfun] = KSPACEFIRSTORDER(..., 'parenv', clu) instead 
            % returns a communicating parallel.Job job and a function to
            % read the ChannelData object from the completed job. Use
            % the submit function to begin running the job.
            % 
            % The settings for clu should must be setup for a communicating
            % parallel.Job as a single MPI job is created for all pulses. 
            % For example, if using a SLURM cluster with multiple GPUs,
            % clu.SubmitArguments should contain the number of GPUs desired
            % for execution e.g. ' --gpus=4' and clu.NumWorkers should
            % equal the maximum number of simulataneous CPUs to use e.g. 4.
            % If more CPUs than GPUs are requested, MATLAB will share GPUs
            % between multiple CPUs. 
            % 
            % When the job has completed successfully, the ChannelData
            % object can be extracted using 'chd = readfun(job)'.
            %
            % A MATLAB parallel.Job is saved until it is deleted, so
            % simulations can be recalled later from the job reference.
            % 
            % chd = KSPACEFIRSTORDER(..., 'parenv', 0) avoids using a
            % parallel.Cluster or parallel.Pool. 
            %
            % The default is the current pool returned by gcp.
            % 
            % [...] = KSPACEFIRSTORDER(..., Name, Value, ...)
            % specifies other Name/Value pairs that are valid for 
            % kWaveArray's constructor or kWave's kspaceFirstOrderND. 
            % 
            % kWave's kspaceFirstOrderND `PMLSize` and `PMLInside`
            % arguments are invalid as they are overwritten by the method.
            %
            % Example:
            % 
            % % Setup a system
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(sscan.size), 1000*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            % rho(Xg == 0 & Zg == 30e-3) = 1000*2; % double the density
            % med = Medium.Sampled(sscan, c, rho);
            % 
            % % Simulate the ChannelData
            % if gpuDeviceCount, dtype = 'gpuArray-single'; 
            % else, dtype = 'single'; end
            % chd = kspaceFirstOrder(us, med, sscan, 'DataCast', dtype);
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(real(chd));
            % 
            % See also ULTRASOUNDSYSTEM/FULLWAVESIM PARALLEL.JOB/FETCHOUTPUTS
            arguments % required arguments
                self (1,1) UltrasoundSystem
                med Medium
                sscan (1,1) ScanCartesian = self.scan
            end
            arguments % keyword arguments for this function
                kwargs.T double {mustBeScalarOrEmpty} = [], % simulation time (s)
                kwargs.PML (1,:) double = [20 56], % (one-sided) PML size range
                kwargs.CFL_max (1,1) double {mustBePositive} = 0.25, % maximum cfl number (for stability)
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate'), % parallel environment for running simulations
                kwargs.ElemMapMethod (1,1) string {mustBeMember(kwargs.ElemMapMethod, ["nearest","linear","karray-direct", "karray-depend"])} = 'nearest', % one of {'nearest'*,'linear','karray-direct', 'karray-depend'}
                kwargs.el_sub_div (1,2) double = self.getLambdaSubDiv(0.1, med.c0), % element subdivisions (width x height)
            end
            arguments % kWaveArray arguments - these are passed to kWaveArray
                karray_args.UpsamplingRate (1,1) double =  10, ...
                karray_args.BLITolerance (1,1) double = 0.05, ...
                karray_args.BLIType (1,1) string {mustBeMember(karray_args.BLIType, ["sinc", "exact"])} = 'sinc', ... stencil - exact or sinc
            end
            arguments % kWave 1.1 arguments - these are passed to kWave
                kwave_args.CartInterp (1,1) string {mustBeMember(kwave_args.CartInterp, ["linear", "nearest"])}
                kwave_args.CreateLog (1,1) logical
                kwave_args.DataCast (1,1) string = 'gpuArray-single'
                kwave_args.DataRecast (1,1) logical = false
                kwave_args.DisplayMask % must be 'off' or a mask the size of sensor.mask
                kwave_args.LogScale (1,1) logical
                kwave_args.MeshPlot(1,1) logical
                kwave_args.MovieArgs cell
                kwave_args.MovieName (1,1) string = 'kwave-sim'
                kwave_args.MovieType (1,1) string {mustBeMember(kwave_args.MovieType, ["frame", "image"])}
                kwave_args.PlotFreq (1,1) {mustBeInteger, mustBePositive}
                kwave_args.PlotLayout (1,1) logical
                kwave_args.PlotPML (1,1) logical = false
                kwave_args.PlotScale = 'auto' % must be 'auto' or [min, max]
                kwave_args.PlotSim (1,1) logical = false
                kwave_args.PMLAlpha (1,1) double
                % kwave_args.PMLInside % set by subroutines
                % kwave_args.PMLSize % set by subroutines
                kwave_args.RecordMovie (1,1) logical = false
                kwave_args.SaveToDisk (1,1) string
                kwave_args.Smooth (1,:) logical % can be (1,3) to smooth {p0, c, rho} separately, or (1,1) for all
                kwave_args.StreamToDisk (1,1) {mustBeNumericOrLogical}
            end

            % only supported with tx == rx for now
            assert(self.tx == self.rx, 'Transmitter and receiver must be identical.')

            % parse
            if isinf(sscan.dx), kwargs.el_sub_div(1) = 1; end % don't use sub-elements laterally for 1D sims
            if isinf(sscan.dy), kwargs.el_sub_div(2) = 1; end % don't use sub-elements in elevation for 2D sims

            % start measuring total execution time
            tt_kwave = tic;

            % get the kWaveGrid
            % TODO: check that the stpe sizes are all equal - this is
            % required by kWaveArray
            [kgrid, Npml] = getkWaveGrid(sscan, 'PML', kwargs.PML);

            % get the kWave medium struct
            kmedium = getMediumKWave(med, sscan);

            % make an iso-impedance medium
            kmedium_iso = kmedium;
            kmedium_iso.density = (med.c0 * med.rho0) ./ kmedium.sound_speed;

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
                    c0map = med.c0;

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
                    karray_args.BLIType = char(karray_args.BLIType);
                    karray_opts = struct2nvpair(karray_args);
                    karray = kWaveArray(self.rx, kgrid.dim, kgrid_origin, karray_opts{:});
                    mask = karray.getArrayBinaryMask(kgrid);

                    % assign source for each transmission (J' x T' x V)
                    % TODO: abstract this logic to come from transducer directly so
                    % it can be used in fullwave or some other FDTD method
                    switch kwargs.ElemMapMethod
                        case 'karray-direct'
                            % define locally
                            vec = @(x) x(:);

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
                            psig = cat(3, psig{:}); % (J' x T' x V)
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
                kwargs.T = 2 * (vecnorm(range([sscan.xb; sscan.yb; sscan.zb], 2),2,1) ./ med.c0);
            end
            Nt = 1 + floor((kwargs.T / kgrid.dt) + max(range(t_tx,1))); % number of steps in time
            kgrid.setTime(Nt, kgrid.dt);

            % get the receive impulse response function
            rx_imp = copy(self.rx.impulse);
            rx_imp.fs = fs_;
            t_rx = rx_imp.time;
            rx_sig = gather(real(rx_imp.samples(:)));

            % simulation start time
            t0 = gather(t_tx(1) + t_rx(1));

            fprintf('There are %i points in the grid which is %0.3f megabytes per grid variable.\n', ...
                kgrid.total_grid_points, kgrid.total_grid_points*4/2^20);

            % set global arguments: these are always overridden
            kwave_args.PMLInside = false;
            kwave_args.PMLSize = Npml(1:kgrid.dim);

            % make string inputs character arrays to avoid compatability
            % issues
            for f = string(fieldnames(kwave_args))'
                if isstring(kwave_args.(f))
                    kwave_args.(f) = char(kwave_args.(f));
                end
            end
            
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
                    rx_args = {kgrid, karray};
                    % proc_fun = @(x) gather(convn( karray.combineSensorData(kgrid, x.ux).', rx_sig, 'full'));
                case {'karray-direct'}
                    rx_args = {elem_weights};
                    % proc_fun = @(x) gather(convn(x.ux.' * full(elem_weights)), rx_sig, 'full'); % (J' x T)' x (J' x N) -> T x N
                case {'nearest'}
                    rx_args = {};
                    %proc_fun = @(x) gather(convn(x.ux.', rx_sig, 'full')); % -> (T x N)
                case 'linear'
                    % create the advanced impulse response function with
                    % which to convolve the output
                    vec = @(x)x(:);
                    N = self.rx.numel;
                    rx_sig = gather(real(rx_imp.sample(t_rx(:)' + el_dist(:)/c0map))); % J'' x T''
                    el_map_el = ((1:N) == vec(ones(size(el_ind)) .* (1:N)))'; % map from convolved samples to elements
                    rx_args = {el_weight, el_map_grd, el_map_el};
                    % proc_fun = @(x) gather(...
                    %     (el_map_el * (el_weight(:) .* convd(el_map_grd' * x.ux, rx_sig, 2, 'full'))).' ... % [(N x J'') x [[(J'' x J') x (J' x T')] x (T' x T | J'')]]' -> (T x N)
                    %     ...
                    %     ); % [(N x J'') x (J'' x T)]' -> T x N

                otherwise, warning('Unrecognized mapping option - mapping to grid pixels by default.');
                    rx_args = {};
                    % proc_fun = @(x) gather(convn(x.ux.', rx_sig, 'full'));
            end

            % choose where to compute the data
            % parfor behaviour - for now, only parenv == 0 -> parfor
            elemmethod = kwargs.ElemMapMethod;
            
            % get the arguments to run the sim
            args = {...
                kspaceFirstOrderND_, ...
                kgrid, kmedium, kmedium_iso, ksource, ksensor, kwave_args_, ...
                rx_sig, elemmethod, rx_args, ...
                t0, fs_ ...
                };

            if nargout < 2 || ~isa(kwargs.parenv, 'parallel.Cluster') % no job/cluster requested
                args{end+1} = kwargs.parenv; % parfor arg: execution environment

                chd = UltrasoundSystem.kspaceRunSim(args{:});

                % TODO: make reports optional
                fprintf(string(self.sequence.type) + " k-Wave simulation completed in %0.3f seconds.\n", toc(tt_kwave));
            else
                % set the parfor options argument
                args{end+1} = Inf; % parfor arg: max number of workers

                % make a job on the cluster
                % TODO: use a modified independent launch script to create
                % tasks but without massive memory overhead.
                clu = kwargs.parenv;
                job = createCommunicatingJob(clu, 'AutoAddClientPath', true, 'AutoAttachFiles',true, 'Type', 'Pool');
                job.createTask(@UltrasoundSystem.kspaceRunSim, 1, args, 'CaptureDiary',true);

                % map the outputs
                chd = job;
                readfun = @(job) subsref(job.fetchOutputs(), substruct('{}', {1}));
            end
        end
    end
    methods(Static, Hidden)
        function y = kspaceFirstOrderPostProc(x, rx_sig, method, varargin)
            % processing step: get sensor data, enforce on CPU (T x N)
            switch method
                case 'karray-depend'
                    [kgrid, karray] = deal(varargin{:});
                    y = gather(convn( karray.combineSensorData(kgrid, x.ux).', rx_sig, 'full'));
                case {'karray-direct'}
                    elem_weights = deal(varargin{:});
                    y = gather(convn(x.ux.' * full(elem_weights), rx_sig, 'full')); % (J' x T)' x (J' x N) -> T x N
                case {'nearest'}
                    y = gather(convn(x.ux.', rx_sig, 'full')); % -> (T x N)
                case 'linear'
                    [el_weight, el_map_grd, el_map_el] = deal(varargin{:});
                    % create the advanced impulse response function with
                    % which to convolve the output
                    y = gather(... [(N x J'') x [[(J'' x J') x (J' x T')] x (T' x T | J'')]]' -> (T x N)
                        (el_map_el * (el_weight(:) .* convd(el_map_grd' * x.ux, rx_sig, 2, 'full'))).' ... % 
                        ); % [(N x J'') x (J'' x T)]' -> T x N

                otherwise, warning('Unrecognized mapping option - mapping to grid pixels by default.');
                   y = gather(convn(x.ux.', rx_sig, 'full'));
            end
        end
    
        function chd = kspaceRunSim(kspaceFirstOrderND_, ...
                kgrid, kmedium, kmedium_iso, ksource, ksensor, kwave_args_, ...
                rx_sig, elemmethod, rx_args, t0, fs_, W ...
                )

            Np = numel(ksource);
            parfor (puls = 1:Np, W)
                % TODO: make this part of some 'info' logger or something
                fprintf('\nComputing pulse %i of %i\n', puls, Np);
                tt_pulse = tic;

                % simulate
                % TODO: try to run these two in parallel on the same GPU?
                sensor_data     = kspaceFirstOrderND_(kgrid, kmedium    , ksource(puls), ksensor, kwave_args_{puls}{:}); %#ok<PFBNS>
                sensor_data_iso = kspaceFirstOrderND_(kgrid, kmedium_iso, ksource(puls), ksensor, kwave_args_{puls}{:});

                % Process the simulation data
                out{puls} = UltrasoundSystem.kspaceFirstOrderPostProc(sensor_data    , rx_sig, elemmethod, rx_args{:}) ...
                          - UltrasoundSystem.kspaceFirstOrderPostProc(sensor_data_iso, rx_sig, elemmethod, rx_args{:}); %#ok<PFBNS> data is small

                % report timing % TODO: make this part of some 'info' logger or something
                fprintf('\nFinished pulse %i of %i\n', puls, Np);
                toc(tt_pulse)
            end

            % create ChannelData objects
            chd = ChannelData('data', cat(3, out{:}), 't0', t0, 'fs', fs_);
        end
    end
    
    % Beamforming
    methods
        function b = DAS(self, chd, c0, kwargs)
            % DAS - Delay and sum beamformer
            %
            % b = DAS(us, chd) performs delay-and-sum beamforming on 
            % the ChannelData chd. The ChannelData must conform to the 
            % delays given by the Sequence us.sequence. The output is the
            % image defined on the Scan us.scan. 
            %
            % b = DAS(self, chd, c0) uses a beamforming sound speed of c0. 
            % c0 can be a scalar or an NDarray that is broadcastable to 
            % size (I1 x I2 x I3) where  [I1, I2, I3] == self.scan.size. 
            % The default is self.sequence.c0.            
            % 
            % b = DAS(..., 'apod', apod) uses an apodization matrix
            % of apod. It must be singular in the receive dimension, the 
            % transmit dimension, or all image dimensions. The apodization
            % matrix must be broadcastable to size (I1 x I2 x I3 x N x M)
            % where [I1, I2, I3] == self.scan.size, N is the number of 
            % receive elements, and M is the number of transmits.
            % 
            % b = DAS(..., 'keep_tx', true) preserves the transmit
            % dimension in the output image b.
            %
            % b = DAS(..., 'keep_rx', true) preserves the receive
            % dimension in the output image b.
            %
            % b = DAS(..., 'fmod', fc) upmixes the data at a modulation
            % frequency fc. This undoes the effect of demodulation at the
            % same freuqency.
            %
            % b = DAS(..., 'interp', method) specifies the method for
            % interpolation. Support is provided by interp1 on the CPU or
            % is restricted to one of {'nearest', 'linear', 'cubic'} on the
            % GPU. The default is 'cubic'.
            % 
            % b = DAS(..., 'prec', type) uses the datatype type for 
            % computing the image. The type can be one of 
            % {'double', 'single'*, 'halfT'}. The precision of the channel 
            % data samples, apodization weights, and the positions / times 
            % are all affected. 
            % 
            % For half precision, the channel data samples and weights are 
            % half precision but the position and  time variables are 
            % computed in single precision. The <a href="matlab:web('https://github.com/thorstone25/halfT')">halfT project</a> must be on 
            % the path. 
            % 
            % b = DAS(..., 'device', index) uses the gpuDevice with ID
            % index. If index == -1, the current device is used. Index == 0
            % avoids using a gpuDevice. The default is -1 if a gpu is
            % available or 0 if no gpu is available.
            %
            % DAS is similar to BFDAS, but is more computationally 
            % efficient at the cost of code readability and interpolation 
            % methods because it avoids calling the ChannelData/sample 
            % method.
            %
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % 
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % chd = hilbert(zeropad(singleT(chd), 0, max(0, chd.T - 2^9))); % precondition the data
            % b = DAS(us, chd); % beamform the data
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also BFDAS BFADJOINT BFEIKONAL CHANNELDATA/SAMPLE
            
            arguments
                self (1,1) UltrasoundSystem
                chd ChannelData
                c0(:,:,:,1,1) {mustBeNumeric} = self.sequence.c0
                kwargs.fmod (1,1) {mustBeNumeric} = 0 
                kwargs.prec (1,1) string {mustBeMember(kwargs.prec, ["single", "double", "halfT"])} = "single"
                kwargs.device (1,1) {mustBeInteger} = -1 * logical(gpuDeviceCount)
                kwargs.apod {mustBeNumericOrLogical} = 1
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.keep_tx (1,1) logical = false
                kwargs.keep_rx (1,1) logical = false
            end

            % parse inputs
            [sumtx, sumrx] = deal(~kwargs.keep_tx, ~kwargs.keep_rx);

            % if we want to keep tx, we must keep rx too because of
            % available DAS functions
            if ~sumtx && sumrx
                error('Unable to keep transmit dimension but not receive dimension. Try bfDAS.'); 
            end

            % make sure t0 is a scalar
            if ~isscalar(chd.t0), chd = rectifyt0(chd); end
            
            % get positions of the imaging plane 
            [X, Y, Z, image_size] = self.scan.getImagingGrid();

            % reshape into I x N x M
            apod_args = {'apod', kwargs.apod, 'modulation', kwargs.fmod};
            
            % convert to x/y/z in 1st dimension
            P_im = permute(cat(4, X, Y, Z),[4,1,2,3]); % 3 x I1 x I2 x I3 == 3 x [I]
            
            % get positions of the aperture(s)
            P_tx = self.tx.positions(); % cast(self.tx.positions(), 'like', time(end)); % 3 x M
            P_rx = self.rx.positions(); % cast(self.rx.positions(), 'like', time(end)); % 3 x N
            
            % get the beamformer arguments
            dat_args = {chd.data, chd.t0, chd.fs, c0, 'device', kwargs.device, 'position-precision', kwargs.prec}; % data args
            if isfield(kwargs, 'interp'), interp_args = {'interp', kwargs.interp}; else,  interp_args = {}; end
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
        
        function chd = focusTx(self, chd0, seq, kwargs)
            % FOCUSTX - Synthesize transmits
            %
            % chd = FOCUSTX(self, chd0) focuses the FSA ChannelData chd0 by
            % linearly synthesizing transmits (i.e. delay and sum across 
            % transmits).
            %
            % chd = FOCUSTX(self, chd0, seq) uses the Sequence seq to focus
            % the data. The default is self.sequence.
            %
            % chd = FOCUSTX(..., Name, Value, ...) uses name-value pairs to
            % specify 
            %
            % chd = FOCUSTX(..., 'interp', method) specifies the method for
            % interpolation. Support is provided by the ChannelData/sample
            % method. The default is 'cubic'.
            %
            % chd = FOCUSTX(..., 'interp', freq, 'length', L) expands the 
            % data to length L before interpolating. For frequency domain
            % interpolation, L must be large enough to avoid aliasing and
            % ringing artefacts.
            %
            % chd = FOCUSTX(..., 'interp', freq, 'length', 'pow2') expands 
            % the data to the next power of 2. This accelerates the fft
            % computation, but uses more memory.
            %
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % 
            % % Compute the data for an FSA acquistion
            % us.sequence = Sequence('type', 'FSA', 'c0', us.sequence.c0, 'numPulse', us.xdc.numel);
            % chd = greens(us, scat); % compute the response
            %
            % % Create plane-wave data by synthesizing the transmits with
            % plane-wave delays
            % seq_pw = SequenceRadial('type', 'PW', 'c0', us.sequence.c0, ...
            %  'angles', -25:0.5:25, 'ranges', 1); % plane-wave sequence
            % chd_pw = focusTx(us, chd, seq_pw); % synthesize transmits
            %
            % % Beamform the plane-wave data
            % us.sequence = seq_pw;
            % b = DAS(us, chd_pw); % beamform the data
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also CHANNELDATA/SAMPLE

            arguments
                self (1,1) UltrasoundSystem
                chd0 ChannelData
                seq (1,1) Sequence = self.sequence;
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.length string {mustBeScalarOrEmpty} = [];
            end

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
            L = kwargs.length;
            if kwargs.interp == "freq",
                if isempty(L),               L = chd.T; % default
                elseif L == "min",           L = chd.T; % strategy
                elseif L == "pow2",          L = 2^(nextpow2(chd.T)); % strategy
                elseif ~isnan(str2double(L)),L = str2double(L); % size
                    if L < chd.T, warning('Signal length may be too short!'); end % soft error
                else
                    error("L must be a length or one of {'min' | 'pow2'}"); % unrecognized input
                end
                
                % expand to match requested FFT length (do not shrink)
                chd = zeropad(chd,0,max(0, L - chd.T)); 
            end
            
            % align dimensions
            D = 1+max(3,ndims(chd.data)); % get a free dimension for M'
            tau_focal = swapdim(tau_focal, [1,2], [chd.mdim, D]); % move data
            apod      = swapdim(apod     , [1,2], [chd.mdim, D]); % move data

            % sample and store
            z = chd.sample(chd.time - tau_focal, kwargs.interp, apod, chd.mdim); % sample (perm(T' x N x 1) x F x ... x M')
            z = swapdim(z, chd.mdim, D); % replace transmit dimension (perm(T' x N x M') x F x ...)
            chd.data = z; % store output channel data % (perm(T' x N x M') x F x ...)
        end
        
        function b = bfAdjoint(self, chd, c0, kwargs)
            % BFADJOINT - Adjoint method beamformer
            %
            % b = BFADJOINT(self, chd) beamforms the ChannelData chd using
            % an adjoint matrix method. This method computes the
            % inner-product of the normalized transmitted wave with the
            % received data in the frequency domain.
            % 
            % b = BFADJOINT(self, chd, c0) uses a beamforming sound speed 
            % of c0. c0 can be a scalar or an NDarray that is broadcastable
            % to size (I1 x I2 x I3) where  [I1, I2, I3] == self.scan.size. 
            % The default is self.sequence.c0.
            %             
            % b = BFADJOINT(..., 'apod', apod) uses an apodization matrix
            % of apod. It must be singular in the receive dimension, the 
            % transmit dimension, or all image dimensions. The apodization
            % matrix must be broadcastable to size (I1 x I2 x I3 x N x M)
            % where [I1, I2, I3] == self.scan.size, N is the number of 
            % receive elements, and M is the number of transmits.
            % 
            % b = BFADJOINT(..., 'Nfft', K) uses a K-point FFT when
            % converting the received to the frequency domain.
            %
            % b = BFADJOINT(..., 'keep_tx', true) preserves the transmit
            % dimension in the output image b.
            %
            % b = BFADJOINT(..., 'keep_rx', true) preserves the receive
            % dimension in the output image b.
            %
            % b = BFADJOINT(..., 'fmod', fc) upmixes the data at a
            % modulation frequency fc. This undoes the effect of
            % demodulation at the same frequency.
            %
            % b = BFADJOINT(..., 'bsize', B) uses an block size of B when
            % vectorizing computations. A larger block size will run
            % faster, but use more memory. The default is chosen
            % heuristically.
            %
            % b = BFADJOINT(..., 'fthresh', thresh) neglects frequencies
            % where the maximum power across the apertures is less than 
            % thresh dB-down with respect to the maximum power. This can 
            % accelerate computation but with a loss of accuracy. The
            % default is -Inf.
            %
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % 
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % chd = hilbert(zeropad(singleT(chd), 0, max(0, chd.T - 2^9))); % precondition the data
            % b = bfAdjoint(us, chd); % beamform the data
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also BFEIKONAL BFDAS DAS FOCUSTX

            % TODO: test for tall types where receive diomension is tall -
            % should work ...
            arguments
                self (1,1) UltrasoundSystem
                chd ChannelData
                c0 (:,:,:,1,1) {mustBeNumeric} = self.sequence.c0
                kwargs.fmod (1,1) {mustBeNumeric} = 0 % modulation frequency
                kwargs.fthresh (1,1) {mustBeReal} = -Inf; % threshold for including frequencies
                kwargs.apod {mustBeNumericOrLogical} = 1; % apodization matrix (I1 x I2 x I3 x N x M)
                kwargs.Nfft (1,1) {mustBeInteger, mustBePositive} = chd.T; % FFT-length
                kwargs.keep_tx (1,1) logical = false % whether to preserve transmit dimension
                kwargs.keep_rx (1,1) logical = false % whether to preserve receive dimension
                kwargs.bsize (1,1) double {mustBeInteger, mustBePositive} = max(1,floor(1*(2^30 / (4*chd.N*self.scan.nPix*8)))); % vector computation block size
                kwargs.verbose (1,1) logical = true 
                % heuristic: 1/4 Gibibyte limit on the size of the delays
            end

            % validate precision: doesn't work for halfT
            % TODO: switch to single precision compute, half precision
            % storage
            if classUnderlying(chd) == "halfT"
                warning('QUPS:bfAdjoint:InsufficientPrecision', ...
                    'Half precision data is insufficient for frequency-domain beamforming.' ...
                    );
            end

            % validate sequence/transducer: doesn't work for
            % non-linear/focused - not sure why yet - seems like a
            % difficult to trace phase error bug.
            if ~isa(self.tx, 'TransducerArray') && ~isa(self.rx, 'TransducerArray') && self.sequence.type == "VS"
                warning('QUPS:bfAdjoint:UnsupportedSequence', ...
                    'This function is unsupported for focused transmits with non-linear transducers.' ...
                    );
            end

            % parse inpcccuts
            sumtx = ~kwargs.keep_tx;
            sumrx = ~kwargs.keep_rx;
            fmod = kwargs.fmod;

            % move the data to the frequency domain, being careful to
            % preserve the time axis
            K = kwargs.Nfft; % DFT length
            f = shiftdim(chd.fs * (0 : K - 1)' / K, 1-chd.tdim); % frequency axis
            df = chd.fs / K; % frequency step size
            x = chd.data; % reference the data perm(K x N x M) x ...
            x = x .* exp( 2i*pi*fmod .* chd.time); % remodulate data
            x = fft(x,K,chd.tdim); % get the fft
            x = x .* exp(-2i*pi*f    .* chd.t0  ); % phase shift to re-align time axis

            % choose frequencies to evaluate
            xmax = max(x, [], chd.tdim); % maximum value per trace
            f_val = mod2db(x) - mod2db(xmax) >= kwargs.fthresh; % threshold
            f_val = f_val & f < chd.fs / 2; % positive frequencies only
            f_val = any(f_val, setdiff(1:ndims(x), chd.tdim)); % evaluate only freqs across aperture/frames that is above threshold

            % get the pixel positions
            D = max(4, gather(ndims(chd.data))); % >= 4
            Pi = self.scan.getImagingGrid();
            Pi = cellfun(@(x) {shiftdim(x, -D)}, Pi); % place I past max data dims 5-7
            Pi = cat(1,Pi{:});     % 3 x 1 x 1 x 1 x ... x [I]
            c0 = shiftdim(c0, -D); % 1 x 1 x 1 x 1 x ... x [I]

            % get the receive apodization, spliced if it can be applied
            apod = kwargs.apod;
            [a_n, a_m, a_mn] = deal(1);
            if all(size(apod, 1:3) == 1) % image is scalar, apply to data
                a_mn = shiftdim(apod, 3); % N x V
            elseif size(apod, 4) == 1 % receive is scalar, apply over tx
                ord = [D+(1:4), 2]; % send dims 1-5 here
                ord = [ord, setdiff(1:max(ord, 5), ord)]; % complete set of dimensions
                a_m = ipermute(apod, ord); % 1 x V x 1 x 1 x ... x [I] x 1
            elseif size(apod, 5) == 1 % transmit is scalar, apply over rx
                ord = [D+(1:3), 2]; % send dims 1-4 here
                ord = [ord, setdiff(1:max(ord, 5), ord)]; % complete set of dimensions
                a_n = ipermute(apod, ord); % 1 x N x 1 x 1 x ... x [I]
            else % none of the following is true: this request is excluded for now
                error('Unable to apply apodization due to size constraints. Apodization must be scalar in the transmit dimension, receive dimension, or all image dimensions.')
            end
            
            % get the delays for the transmit/receive green's matrix
            % kernels
            tau_tx = vecnorm(self.tx.positions() - Pi,2,1) ./ c0; % 1 x M x 1 x 1 x ... x [I]
            tau_rx = vecnorm(self.rx.positions() - Pi,2,1) ./ c0; % 1 x N x 1 x 1 x ... x [I]

            % get the transmit steering vector weights and delays
            del_tx  = self.sequence.delays(self.tx);      % M x V
            apod_tx = self.sequence.apodization(self.tx); % M x V
            del_tx = del_tx - self.sequence.t0Offset(); % offset to transducer origin

            % transform to frequency step kernels
            w_rx    = exp(-2i*pi*df.*tau_rx); %  receive greens function
            w_tx    = exp(-2i*pi*df.*tau_tx); % transmit greens function
            w_steer = exp(-2i*pi*df.*del_tx); % transmit steering delays

            % cast data type for efficency
            [w_tx, w_rx, w_steer, apod_tx] = dealfun(@(w) cast(w, 'like', real(x)), w_tx, w_rx, w_steer, apod_tx);
            b = repmat(cast(0, 'like', x), [1, size(Pi,2:ndims(Pi))]);

            % splice data 
            k = gather(find(f_val)'); % skip unimportant frequencies
            k = num2cell(reshape([k(:); nan([-mod(numel(k),-kwargs.bsize),1])], kwargs.bsize, []), 1)'; % expand and make blocks
            k{end}(isnan(k{end})) = []; % delete invalid entries
            xk = cellfun(@(k) sub(x,k,chd.tdim), k, 'UniformOutput',false);
            chd_ord = [chd.ndim, chd.mdim, chd.tdim]; % permution order

            % beamform for a block of frequencies at a time
            if kwargs.verbose, hw = waitbar(0,'Beamforming ...'); end

            % DEBUG: plot
            % figure; h = imagesc(squeeze(zeros(self.scan.size))); colorbar; colormap jet;

            % TODO: parallelize for thread-based pools only
            % parfor (ik = 1:numel(k))
            for ik = 1:numel(k)
                % get discrete frequency index 
                k_ = shiftdim(k{ik}, -2); % 1 x 1 x F

                % report progress
                lbl = "Beamforming freqs: " + min(gather(f(k_))) + " - " + max(gather(f(k_)));% + " MHz";
                if kwargs.verbose && isvalid(hw), waitbar(ik/numel(k), hw, char(lbl)); end
                % fprintf(lbl + "\n");

                % data, in freq. domain (N x V x ...)
                % TODO: adapt for tall types
                xk_ = permute(xk{ik}, [chd_ord, 4:D]);
                
                % compute the greens functions on transmit/receive
                G_tx = w_tx.^(k_-1); % 1 x M x F x 1 x [I]
                G_rx = w_rx.^(k_-1); % 1 x N x F x 1 x [I]

                % compute the inverse steering vector on transmit
                T_tx = apod_tx .* w_steer.^(k_-1); % M x V x F
                A_tx = pagemtimes(G_tx, T_tx); % 1 x V x F x 1 x [I]
                Ainv_tx = pagetranspose(A_tx); % V x 1 x F x 1 x [I] % make a column vector
                Ainv_tx = Ainv_tx ./ vecnorm(Ainv_tx, 2, 1); % normalize the power
                
                % apodize, delay, and sum the data for this frequency
                % only 1 of the a_* will contain the apodization
                if sumrx, yn = a_m .*      pagemtimes(a_n .* conj(G_rx)   , a_mn .* xk_); % 1 x V x F x 1 x [I]
                else,     yn = a_m .* (pagectranspose(a_n .*     (G_rx)) .* a_mn .* xk_); % N x V x F x 1 x [I]
                end
                if sumtx, y  = pagemtimes(yn,                 conj(Ainv_tx));  % [1|N] x 1 x F x 1 x [I]
                else,     y  =           (yn .* pagectranspose(   (Ainv_tx))); % [1|N] x V x F x 1 x [I]
                end

                % integrate over all frequencies
                b = b + sum(y,3); % [1|N] x [1|V] x 1 x 1 x [I]

                % DEBUG: update display
                % h.CData(:) = mod2db(sum(y,3)); drawnow limitrate; 
            end
            if kwargs.verbose && isvalid(hw), close(hw); end

            % move to image dimensions ([I] x ... x perm([1|N] x [1|V] x 1)
            b = swapdim(ipermute(b,[chd_ord, 4:(D+3)]), 1:3, D+(1:3));
        end    

        function b = bfEikonal(self, chd, medium, cscan, kwargs)
            % BFEIKONAL - Delay-and-sum beamformer with Eikonal delays
            %
            % b = BFEIKONAL(self, chd, medium) creates a b-mode
            % image b from the ChannelData chd and Medium medium using the
            % delays given by the solution to the eikonal equation defined
            % on the ScanCartesian cscan == self.scan.
            % 
            % The transmitter and receiver must fall within the cscan. The 
            % step size in each dimension must be identical. The eikonal 
            % equation is solved via the fast marching method.
            %
            % b = BFEIKONAL(self, chd, medium, cscan) uses the given 
            % ScanCartesian cscan instead of self.scan.
            % 
            % b = BFEIKONAL(..., Name,Value, ...) defines additional
            % parameters via Name/Value pairs
            %
            % b = BFEIKONAL(..., 'keep_rx', true) preserves the receive 
            % dimension
            %
            % b = BFEIKONAL(..., 'keep_tx', true) preserves the tranmit 
            % dimension
            %
            % b = BFEIKONAL(..., 'fmod', fc) upmixes the data at a
            % modulation frequency fc. This undoes the effect of
            % demodulation at the same frequency.
            %
            % b = BFEIKONAL(..., 'apod',A) uses an apodization defined by
            % the ND-array A. A must be broadcastable to size
            % I1 x I2 x I3 x N x M where I1 x I2 x I3 is the size of the
            % image, N is the number of receivers, and M is the number of
            % transmits.
            %
            % b = BFEIKONAL(..., 'bsize',B) computes using a B 
            % frequencies at a time in order to limit memory usage. A lower
            % B uses less memory to prevent OOM errors but a higher B 
            % yields better computer performance.
            %
            % b = BFEIKONAL(..., 'parenv', clu) or
            % b = BFEIKONAL(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to
            % parallelize computations. parallel.ThreadPools will be
            % ignored due to mex function restrictions.
            % 
            % b = BFEIKONAL(..., 'parenv', 0) avoids using a 
            % parallel.Cluster or parallel.Pool. Use 0 when operating on a 
            % GPU or if memory usage explodes on a parallel.ProcessPool.
            %
            % b = BFEIKONAL(..., 'interp',method, ...) specifies the 
            % method for interpolation. Support is provided by the 
            % ChannelData/sample method. The default is 'cubic'.
            %
            % Example:
            % 
            % % This example requires kWave
            % if ~exist('kWaveGrid', 'class')
            %     warning('kWave must be on the path to run this example.');
            % end
            % 
            % % Setup a system
            % sscan = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray('numel', 16, 'fc', 3e6, 'bw', [1.5, 4.5]*1e6, 'pitch', 1.5e3/3e6);
            % seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', 1500);
            % us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
            % 
            % % Create a Medium to simulate
            % [c0, rho0] = deal(1.5e3, 1e3); 
            % [c, rho] = deal(c0*ones(sscan.size), rho0*ones(sscan.size));
            % [Xg, ~, Zg] = sscan.getImagingGrid();
            % 
            % % Define isoimpedance layers
            % z0 = rho0 * c0; % ambient impedance
            % [c(Zg > 15e-3), rho(Zg > 15e-3)] = deal(1.4e3, z0/1.4e3); % isoimpedance
            % [c(Zg > 25e-3), rho(Zg > 25e-3)] = deal(1.6e3, z0/1.6e3); % isoimpedance 
            % [c(Zg > 35e-3), rho(Zg > 35e-3)] = deal(1.4e3, z0/1.4e3); % isoimpedance
            % [c(Zg > 45e-3), rho(Zg > 45e-3)] = deal(1.5e3, z0/1.5e3); % isoimpedance
            % 
            % % Define density scatterers
            % rho(Xg == 0 & Zg == 10e-3) = rho0*2; 
            % rho(Xg == 0 & Zg == 20e-3) = rho0*2; 
            % rho(Xg == 0 & Zg == 30e-3) = rho0*2; 
            % rho(Xg == 0 & Zg == 40e-3) = rho0*2; 
            % rho(Xg == 0 & Zg == 50e-3) = rho0*2;
            % 
            % % Construct the Medium
            % med = Medium.Sampled(sscan, c, rho);
            % 
            % % Simulate the ChannelData
            % if gpuDeviceCount, dtype = 'gpuArray-single'; 
            % else, dtype = 'single'; end
            % chd = kspaceFirstOrder(us, med, sscan, 'DataCast', dtype, 'CFL_max', 0.5);
            % 
            % % Beamform
            % b_naive = DAS(us, chd);
            % b_c0    = bfEikonal(us, chd, med, sscan);
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(real(chd));
	    % title('Channel Data')
            % 
            % % Display the images
            % bim_naive = mod2db(b_naive);
            % bim_c0    = mod2db(b_c0   );
            % 
            % figure;
            % subplot(1,2,1);
            % imagesc(us.scan, bim_naive, [-80 0] + max(bim_naive(:)));
            % colormap gray; colorbar;
            % title('Naive Delay-and-Sum');
            % 
            % subplot(1,2,2);
            % imagesc(us.scan, bim_c0   , [-80 0] + max(bim_c0(:)   ));
            % colormap gray; colorbar;
            % title('Eikonal Delay-and-Sum');
            % 
            % See also DAS BFDAS BFADJOINT

            arguments
                self (1,1) UltrasoundSystem
                chd ChannelData
                medium Medium
                cscan (1,1) ScanCartesian = self.scan
                kwargs.fmod (1,1) {mustBeNumeric} = 0
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate')
                kwargs.apod {mustBeNumericOrLogical} = 1;
                kwargs.keep_rx (1,1) logical = false;
                kwargs.keep_tx (1,1) logical = false;
                kwargs.bsize (1,1) double {mustBeInteger, mustBePositive} = max(1,floor(1*(2^30 / (chd.N*self.scan.nPix*8)))); 
                kwargs.verbose (1,1) logical = true;
                % 1 Gibibyte limit on the size of the delays
            end

            % get summation options
            sumtx = ~kwargs.keep_tx;
            sumrx = ~kwargs.keep_rx;
            fmod = kwargs.fmod;

            % get cluster
            parenv = kwargs.parenv; % compute cluster/pool/threads
            % travel times cannot use threadPool because mex call
            if isempty(parenv) || isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool'), parenv = 0; end

            % get worker transfer function
            if(parenv == 0)
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
                "The simulation scan must have equally sized steps in all non-singleton dimensions." ...
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
            gi_opts = {'cubic', 'none'}; % interpolater options
            if kwargs.verbose 
                tt = tic; fprintf('\nComputing Eikonal time delays ... \n');
            end
            parfor (n = 1:chd.N, parenv)
                % fprintf('rx %i\n', n);
                [tau_map_rx] = msfm(squeeze(cnorm.Value), double(Prc(:,n))); %#ok<PFBNS> % travel time to each point
                rx_samp{n} = griddedInterpolant(grd, tau_map_rx,gi_opts{:}); %#ok<PFBNS> % make interpolator on cscan
            end
            if self.tx == self.rx % if apertures are identical, copy
                tx_samp = rx_samp;
            else % else compute for each tx
            parfor (m = 1:chd.M, parenv)
                % fprintf('tx %i\n', m);
                [tau_map_tx] = msfm(squeeze(cnorm.Value), double(Pvc(:,m))); %#ok<PFBNS> % travel time to each point
                tx_samp{m} = griddedInterpolant(grd, tau_map_tx,gi_opts{:}); %#ok<PFBNS> % make interpolator on cscan
            end
            end
            if kwargs.verbose
                fprintf('\nEikonal time delays completed in %0.3f seconds.\n', toc(tt));
            end

            % get the imaging grid
            gi = self.scan.getImagingGrid(); % {I1 x I2 x I3} each
            gi = gi(sel); % select and trim dimensions 

            % splice args
            interp_method = kwargs.interp; 
            apod = kwargs.apod;

            % get sample times for each tx/rx
            tau_rx = cellfun(@(f) f(gi{:}), rx_samp, 'UniformOutput', false); % all receive delays
            tau_tx = cellfun(@(f) f(gi{:}), tx_samp, 'UniformOutput', false); % all transmit delays
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
            kvb = kwargs.verbose; % splice
            if kvb, hw = waitbar(0,'Beamforming ...'); else, hw = []; end % create a wait bar
            parfor (m = 1:Mp, 0) % for each transmit (no cluster because parpool may cause memory issues here)
                % make the eikonal delays and apodization align with the channel data
                tau = tau_rx + tau_tx{m}; % I1 x I2 x I3 x N x M
                if isscalar(apod), a = apod{1}; else, a = apod{m}; end % recieve apodization per transmit (I1 x I2 x I3 x N x M)
                a = ipermute(a, ord); % move time delays / apodization into matching dimensions  
                tau = ipermute(tau, ord); % both: perm(1 x N x M) [x 1 x ... ] x I1 x I2 x I3
                    
                % sample and unpack the data (perm(1 x N x M) [x F x ... ] x I1 x I2 x I3)
                ym = sample(chds(m), tau, interp_method, a, sdim, fmod); % sample and sum over rx/tx maybe
                
                % sum or accumulate over transmit blocks
                if sumtx, b = b + ym; else, bm{m} = ym; end 
                if kvb && isvalid(hw), waitbar((Mp-m+1)/Mp, hw); end % update if not closed: parfor loops go backwards
            end
            if ~sumtx, b = cat(chd.mdim,bm{:}); end % combine if not summing tx
            if kwargs.verbose && isvalid(hw), close(hw); end % close waitbar if not already closed

            % move image output into lower dimension
            if istall(b), b = gather(b); end % we have to gather tall arrays to place anything in dim 1
            b = swapdim(b, 1:3, D+(1:3)); % move image dimensions down (I1 x I2 x I3 [x F x ... ] x perm(1 x N x M))
        end
    
        function b = bfDAS(self, chd, c0, kwargs)
            % BFDAS - Delay-and-sum beamformer
            %
            % b = BFDAS(self, chd, c0) creates a b-mode image b from the 
            % ChannelData chd and sound speed c0.
            % 
            % b = BFDAS(..., Name, Value, ...) passes additional Name/Value
            % pair arguments
            % 
            % b = BFDAS(..., 'apod', apod) uses an apodization matrix
            % of apod. It must be singular in the receive dimension, the 
            % transmit dimension, or all image dimensions. The apodization
            % matrix must be broadcastable to size (I1 x I2 x I3 x N x M)
            % where [I1, I2, I3] == self.scan.size, N is the number of 
            % receive elements, and M is the number of transmits.
            % 
            % b = BFDAS(..., 'keep_tx', true) preserves the transmit
            % dimension in the output image b.
            %
            % b = BFDAS(..., 'keep_rx', true) preserves the receive
            % dimension in the output image b.
            %
            % b = BFDAS(..., 'bsize', B) uses an block size of B to
            % compute at most B transmits at a time. A larger block size 
            % will run faster, but use more memory. The default is chosen
            % heuristically.
            %   
            % b = BFDAS(..., 'fmod', fc) upmixes the data at a modulation
            % frequency fc. This undoes the effect of demodulation at the
            % same frequency.
            %
            % b = BFDAS(..., 'interp', method) specifies the method for
            % interpolation. Support is provided by the ChannelData/sample
            % method. The default is 'cubic'.
            % 
            % b = BFDAS(..., 'device', index) uses the gpuDevice with ID
            % index. If index == -1, the current device is used. Index == 0
            % avoids using a gpuDevice. The default is -1 if a gpu is
            % available.
            %
            % b = BFDAS(..., 'parenv', clu) or
            % b = BFDAS(..., 'parenv', hcp) uses a parallel.Cluster clu 
            % or a parallel.Pool hcp to parallelize computations. 
            % b = BFDAS(..., 'parenv', 0) avoids using a 
            % parallel.Cluster or parallel.Pool. Use 0 when operating on a 
            % GPU or if memory usage explodes on a parallel.ProcessPool.
            % 
            % The default is the current parallel pool returned by gcp.
            % 
            % A parallel.ThreadPool will tend to perform better than a 
            % parallel.ProcessPool because threads are allowed to share 
            % memory.
            %        
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.sequence.c0); % define a point target
            % 
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % b = bfDAS(us, chd); % beamform the data
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also DAS BFADJOINT CHANNELDATA/SAMPLE PARCLUSTER PARPOOL

            arguments
                self (1,1) UltrasoundSystem
                chd ChannelData
                c0 (:,:,:,1,1) {mustBeNumeric} = self.sequence.c0
                kwargs.fmod (1,1) {mustBeNumeric} = 0 
                kwargs.device (1,1) {mustBeInteger} = -1 * logical(gpuDeviceCount)
                kwargs.apod {mustBeNumericOrLogical} = 1
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.keep_tx (1,1) logical = false
                kwargs.keep_rx (1,1) logical = false
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Pool", "parallel.Cluster", "double"])} = gcp('nocreate');
                kwargs.bsize (1,1) double {mustBeInteger, mustBePositive} = max(1,floor(1*(2^30 / (chd.N*self.scan.nPix*8)))); 
                % 1 Gibibyte limit on the size of the delays
            end

            % get cluster
            parenv = kwargs.parenv;
            if isempty(parenv) || isa(chd.data, 'gpuArray'), parenv = 0; end % run on CPU by default

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
            fmod = kwargs.fmod;

            % move image dimensions beyond the data
            D = max([ndims(chd.data), ndims(dv), ndims(dr), ndims(cinv), 5]); % highest dimension of data
            [dv, dr, apod, cinv] = dealfun(@(x) swapdim(x, 1:3, D+(1:3)), dv, dr, apod, cinv);
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
            % for m = 1:numel(chds)
            parfor (m = 1:numel(chds), parenv) % for each transmit
                tau = cinv .* (dvm{m} + dr); % get sample times (1 x 1 x 1 x N x 1 x ... x I1 x I2 x I3)
                if isscalar(am), a = am{1}; else, a = am{m}; end % (1 x 1 x 1 x N x 1 x ... x I1 x I2 x I3)

                % move to permutation of (1 x N x M) - N/M aligned with ChannelData
                tau = swapdim(tau, [4,5], [chds(m).ndim, chds(m).mdim]); 
                a   = swapdim(a  , [4,5], [chds(m).ndim, chds(m).mdim]); 

                % sample, apodize, and sum over rx if requested
                z = sample(chds(m), tau, interp_method, a, sdim, fmod); % (perm(1 x N x M) x F x ... x I1 x I2 x I3)

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
            % RECOMPILECUDA(self, defs) compiles for the compiler 
            % definition structs defs. These structures are generated by
            % the UltrasoundSystem class. The default is all the
            % definitions returned by UltrasoundSystem.getMexFileDefs().
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
        function recompileCUDA(self, defs)
            % RECOMPILECUDA - Recompile CUDA ptx files
            %
            % RECOMPILECUDA(self) recompiles all CUDA files and stores
            % them in self.tmp_folder.
            %
            % RECOMPILECUDA(self, defs) compiles for the compiler 
            % definition structs defs. These structures are generated by
            % the UltrasoundSystem class. The default is all the
            % definitions returned by UltrasoundSystem.getCUDAFileDefs().
            % 
            % See also ULTRASOUNDSYSTEM.RECOMPILEBFCONST
            % ULTRASOUNDSYSTEM.RECOMPILE ULTRASOUNDSYSTEM.RECOMPILEMEX 
            % ULTRASOUNDSYSTEM.GETCUDAFILEDEFS

            arguments
                self (1,1) UltrasoundSystem
                defs (1,:) struct = UltrasoundSystem.getCUDAFileDefs();
            end
            
            % src file folder
            src_folder = UltrasoundSystem.getSrcFolder();
            
            % get the current gpu's compute capability support
            try  
                g = gpuDevice();
                arch = "compute_" + replace(g.ComputeCapability,'.',''); % use the current GPU's CC number
            catch
                warning("Unable to access GPU.");
                arch = string.empty; % don't use this argument
            end

            % compile each
            for d = defs
                % make full command
                com = join(cat(1,...
                    "nvcc ", ...
                    "--ptx " + fullfile(src_folder, d.Source), ...
                    "-arch=" + arch + " ", ... compile for active gpu
                    "-o " + fullfile(self.tmp_folder, strrep(d.Source, '.cu', '.ptx')), ...
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
            N = self.rx.numel; % number of receiver elements
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
                    {"QUPS_T="+chd.T}; ... number of time samples
                    ];
            end
            
            % compile each
            for d = defs
                % make full command
                com = join(cat(1,...
                    "nvcc --ptx ", ...
                    "-arch=native ", ... % for half types: TODO move to compile option
                    fullfile(src_folder, d.Source) + " ", ...
                    "-o " + fullfile(self.tmp_folder, strrep(d.Source, '.cu', '.ptx')), ...
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
        function recompileGREENSCONST(self, scat)
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
            arguments
                self (1,1) UltrasoundSystem
                scat Scatterers {mustBeScalarOrEmpty} = Scatterers.empty
            end

            % src file folder
            src_folder = UltrasoundSystem.getSrcFolder();

            % get the Waveform length
            wv = conv(self.rx.impulse, ...
                conv(self.tx.impulse, self.sequence.pulse, self.fs), ...
                self.fs); % transmit waveform, convolved at US frequency
            wv.fs = self.fs;
            T = length(wv.samples);

            % get the other sizes for greens.cu
            N = self.rx.numel; % number of receiver elements
            M = self.tx.numel; % number of transmits

            % get all source code definitions
            defs = UltrasoundSystem.genCUDAdef_greens();

            % add the defined macros
            defs.DefinedMacros = cat(1, ...
                defs.DefinedMacros, ... keep the current defs
                "QUPS_" + {... prepend 'QUPS_'
                "T="+T,... elements
                "N="+N,... elements
                "M="+M,... transmits
                }');

            % if I is provided, include it
            if ~isempty(scat)
                defs.DefinedMacros = [defs.DefinedMacros; ...
                    {"QUPS_I="+scat.numScat}; ... number of time samples
                    ];
            end

            % compile each
            for d = defs
                % make full command
                com = join(cat(1,...
                    "nvcc --ptx ", ...
                    "-arch=native ", ... % for half types: TODO move to compile option
                    fullfile(src_folder, d.Source) + " ", ...
                    "-o " + fullfile(self.tmp_folder, strrep(d.Source, '.cu', '.ptx')), ...
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
        function defs = getCUDAFileDefs()
            % GETCUDAFILEDEFS - Get the CUDA compilation definition structs
            %
            % defs = GETCUDAFILEDEFS() returns a struct array with fields 
            % and values specifying compilation arguments. These are used
            % by UltrasoundSystem.recompileCUDA to recompile code for the
            % current host machine.
            %
            % See also RECOMPILECUDA ULTRASOUNDSYSTEM.GETMEXFILEDEFS
            
            % get all source code definitions
            defs = [...
                UltrasoundSystem.genCUDAdef_beamform(),...
                UltrasoundSystem.genCUDAdef_interpd(),...
                UltrasoundSystem.genCUDAdef_convd(),...
                UltrasoundSystem.genCUDAdef_greens(),...
                UltrasoundSystem.genCUDAdef_wbilerp(),...
                ];
        end

        function defs = getMexFileDefs()
            % GETMEXFILEDEFS - Get the mex compilation definition structs
            %
            % defs = GETMEXFILEDEFS() returns a struct array with fields 
            % and values specifying compilation arguments. These are used
            % by UltrasoundSystem.recompileMex to recompile code for the
            % current host machine.
            %
            % See also RECOMPILEMEX ULTRASOUNDSYSTEM.GETCUDAFILEDEFS
            
            % get all source code definitions
            defs = [...
                UltrasoundSystem.genMexdef_msfm(), ... % both msfm files
                ];
        end

    end
    
    % source file recompilation definitions
    methods(Static,Hidden)
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

        function d = genCUDAdef_wbilerp()
            % no halp :(

            % filename
            d.Source = {...
                'wbilerp.cu', ...
                }';

            d.IncludePath = {}; % include folders
            d.Libraries = {}; % libraries

            d.CompileOptions = {...  compile options
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

function tmp = mktempdir(),
tmp = tempname(); % new folder
try
    mkdir(tmp); % try to create it
catch % cannot create new folder
    tmp = cachedir(); % reference the cache folder
end
end

function d = cachedir()
    d = fullfile(fileparts(mfilename('fullpath')), 'bin');
end
