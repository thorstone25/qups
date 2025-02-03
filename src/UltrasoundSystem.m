% ULTRASOUNDSYSTEM - Comprehensive ultrasound system definition class
%
% The ULTRASOUNDSYSTEM class is a synthesis class containing the properties
% describing a medical ultrasound system. Once a Transducer, Sequence, and
% Scan are defined, the methods are provided to use a Scatterers or Medium
% to simulate ChannelData, or beamform ChannelData into an image.
%
% Multiple simulators are supported. Most simulators support the arbitrary
% Sequences and arbitrary Transducers allowed in QUPS. All simulators
% support parallel enviroments such as ProcessPools (multiple MATLAB
% instances) and parclusters (for compute clusters) with some supporting
% ThreadPools (multi-threading). External simulators must be installed
% separately. They include:
% 
% * greens - simulate point scatterers via QUPS (GPU-enabled)
% * simus  - simulate point scatterers via MUST
% * calc_scat_all - simulate point scatterers via FieldII, then synthesize transmits
% * calc_scat_multi - simulate point scatterers via FieldII, per transmit
% * kspaceFirstOrder - simulate a medium via K-wave (GPU-enabled)
% * fullwaveSim - simulate a medium via Fullwave (currently non-public)
% 
% Multiple beamformers are provided, all of which are GPU enabled, either
% natively in MATLAB or via ptx binaries, or both. They include:
% 
% * DAS - a performant conventional delay-and-sum beamformer
% * bfDAS - a generic delay-and-sum beamformer
% * bfDASLUT - a delay-and-sum beamformer for custom delays
% * bfEikonal - a delay-and-sum beamformer using eikonal equation delays
% * bfAdjoint - a frequency domain matrix adjoint beamformer
% * bfMigration - a frequency domain stolt's f-k migration beamformer
%
% Most beamformers are pixel-based, so classical scanline-based processing
% must be emulated via apodization ND-arrays. The provided receive
% apodization generators include:
%
% * apScanline - emulate a scan-line beamformer (focal sequences)
% * apMultiline - emulate a multi-line beamformer (focal sequences)
% * apTranslatingAperture - emulate a translating transmit aperture
% * apApertureGrowth - limit the f# on rx
% * apAcceptanceAngle - limit the element-to-pixel angle on rx
% * apTxParallelogram - limit to an illumination parallelogram on tx (plane-waves)
% 
% One can also synthesize new transmit sequences from full synthetic
% aperture (FSA) data or synthesizse FSA data from focused pulses. These
% utilities are provided by:
%
% * focusTx - synthesize a generic transmit sequence from FSA data
% * refocus - synthesize FSA data from a generic transmit sequence
% 
% See also CHANNELDATA TRANSDUCER SEQUENCE SCAN SCATTERERS MEDIUM 

classdef UltrasoundSystem < matlab.mixin.Copyable & matlab.mixin.CustomDisplay
    
    % objects
    properties(Dependent)
        xdc Transducer % Transducer object (if receive and transmit are identical)
    end
    properties
        tx (1,1) Transducer = TransducerArray.P4_2v() % Transducer object (transmit)
        rx (1,1) Transducer = TransducerArray.P4_2v() % Transducer object (receive)
        seq (1,1) Sequence = Sequence() % Sequence object
        scan (1,1) Scan = ScanCartesian() % Scan object
    end
    
    % parameters
    properties
        % FS - Simulation sampling frequency 
        %
        % ULTRASOUNDSYSTEM.FS sets the sampling frequency and data
        % precision for simulation routines.
        %
        % Example:
        %
        % us = UltrasoundSystem('fs', 50e6); % 50MHz
        % 
        % See also GREENS SIMUS CALC_SCAT_MULTI KSPACEFIRSTORDER
        fs (1,1) {mustBePositive} = single(40e6)   % simulation sampling frequency
    end
    
    properties(Dependent)
        fc  {mustBeNumeric} % central operating frequency(ies)
        lambda {mustBeNumeric} % wavelength at the central frequency(ies)
    end
    
    properties(Hidden,NonCopyable,SetAccess=protected)
        tmp_folder (1,1) string = tempname() % temporary folder for compiled binaries
    end

    % aliases
    properties(Hidden, Dependent)
        sequence
    end
        
    % constructor/destructor
    methods
        % constructor
        function us = UltrasoundSystem(kwargs, opts)
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
            % us = ULTRASOUNDSYSTEM(...,'seq', seq) sets the Sequence
            % to be the Sequence seq.
            %
            % us = ULTRASOUNDSYSTEM(...,'scan', scan) sets the Scan to be
            % the Scan scan.
            %
            % us = ULTRASOUNDSYSTEM(...,'fs', fs) sets the simulation 
            % sampling frequency to fs.
            %
            % us = ULTRASOUNDSYSTEM(...,'recompile', true) attempts to
            % recompile mex and CUDA files for the UltrasoundSystem object.
            % The default is false.
            %
            % us = ULTRASOUNDSYSTEM(...,'copybin', true) attempts to find
            % and copy any compiled binaries that exist on the path. If
            % 'recompile' is also true, if any binaries are missing, all
            % binaries will be recompiled. The default is false.
            %
            % Example:
            % UltrasoundSystem('recompile', true)
            % 
            % 
            % See also TRANSDUCER SEQUENCE SCAN
            arguments
                kwargs.?UltrasoundSystem
                opts.recompile (1,1) logical = false
                opts.copybin (1,1) logical = false
            end
            
            % initailize via name-Value pairs
            f = string(fieldnames(kwargs))'; % name-value fields (1 x F)

            % initialize
            for s = f, us.(s) = kwargs.(s); end
            
            % if neither transmit nor receive provided, make them
            % equivalent for convenience.
            if ~isfield(kwargs, 'tx') && ~isfield(kwargs, 'rx'), us.tx = us.rx; end

            % if sampling frequency not provided, set to Nyquist
            if ~isfield(kwargs, 'fs')
                us.fs = 2*single(max([ ...
                    2*us.tx.fc, us.tx.bw, 2*us.rx.fc, us.rx.bw ...
                    ],[],'omitnan'));
            end
            
            % create a default Sequence if none provided
            if ~isfield(kwargs, 'seq')
                us.seq = Sequence(...
                    'type','FSA',...
                    'focus', [0;0;0], ...
                    'numPulse', us.tx.numel ...
                    );
            end

            % if the scan was not provided, create one and set the
            % resolution
            if ~isfield(kwargs, 'scan')
                if isa(us.rx, 'TransducerConvex')
                    us.scan = ScanPolar('origin', us.rx.center);
                    us.scan.r = us.scan.r + norm(us.rx.center);
                    [us.scan.dr, us.scan.dy] = deal(min(us.lambda) / 4);
                    us.scan.da = 1; % every 1 degree
                elseif isa(us.rx, 'TransducerMatrix')
                    us.scan = ScanCartesian();
                    us.scan.y = us.scan.x;
                    [us.scan.dx, us.scan.dy, us.scan.dz] = deal(min(us.lambda) / 4);
                else
                    us.scan = ScanCartesian();
                    [us.scan.dx, us.scan.dy, us.scan.dz] = deal(min(us.lambda) / 4);
                end
            end

            % shadow with the (newly created) temp folder for binaries and 
            % const-compiled code
            us.tmp_folder = tempname;
            mkdir(  us.tmp_folder);
            addpath(us.tmp_folder);

            % no recompilation on a thread pool
            if parallel.internal.pool.isPoolThreadWorker
                if opts.recompile, warning("Cannot recompile on a thread worker!"); end
                return;
            end

            % short circuit if no binary action is requested
            if ~opts.copybin && ~opts.recompile, return; end
            
            % copy code or recompile it
            if gpuDeviceCount % only do CUDA stuff if there's a MATLAB-compatible GPU
                defs = us.genCUDAdefs();
                if opts.copybin
                fls = arrayfun(@(d) string(strrep(d.Source, '.cu', '.ptx')), defs);
                e = logical(arrayfun(@(x)exist(x, 'file'), fls)); % already exists?
                s = arrayfun(@(fl) copyfile(which(fl), fullfile(us.tmp_folder, fl)), fls(e)); % move there if not?
                end
                if opts.recompile && (~opts.copybin || any(~s)), us.recompileCUDA(); end % attempt to recompile code
            end

            % copy code or recompile it
            defs = us.genMexdefs();
            if opts.copybin
            fls = arrayfun(@(d) string(strrep(d.Source, 'c', mexext())), defs);
            e = logical(arrayfun(@(x)exist(x, 'file'), fls)); % already exists?
            s = arrayfun(@(fl) copyfile(which(fl), fullfile(us.tmp_folder, fl)), fls(e)); % move there if not?
            end
            if opts.recompile && (~opts.copybin || any(~s)), us.recompileMex(); end % attempt to recompile code
        end
    end

    % destructor
    methods(Access=protected)
        function delete(us)
            % DELETE - Destroy an UltrasoundSystem ... programatically.
            %
            % On object destruction, any temporary directories are removed.
            %
            % Note: argument validation will invalidate this method as a
            % destructor.
            % 
            % See also HANDLE

            % if we made a temp folder, clean it up
            tmp = us.tmp_folder;
            if ~isempty(tmp) && exist(tmp, 'dir')
                rmpath(tmp) % remove from the path
                list = dir(tmp); % all files in the folder
                nms = string({list.name}); % get the file/folder names
                nms(nms == "." | nms == "..")  = []; % delete self references
                
                % check that it's only ptx and mex files we are deleting
                assert(all(...
                    endsWith(nms, [".ptx", string(mexext())]) ...
                    ), "QUPS:UltrasoundSystem:dirtyDirectory", ...
                    "Deletion of directory " + tmp + ...
                    " aborted due to presence of unexpected folders/files:"+newline...
                    + join(nms, ", ") ...
                    );

                % rmdir(tmp, 's'); % recursive deletion - dangerous

                % safe: remove specific files first, then (attempt) the folder
                delete(fullfile(tmp, "*" + ".ptx")); % remove any ptx files
                delete(fullfile(tmp, "*" + string(mexext()))); % remove any mex files
                rmdir(tmp); % remove main folder
                % disp("[DEBUG]: Deleted " + tmp);
            end
        end
    end

    methods
        % convert to a structure to remove class info
        function s = obj2struct(us)
            % OBJ2STRUCT - Convert a QUPS object into a native MATLAB struct
            %
            % us = OBJ2STRUCT(us) converts the UltrasoundSystem us and all
            % of it's properties into native MATLAB structs.
            %
            % Example:
            %
            % % Create an UltrasoundSystem
            % us = UltrasoundSystem() 
            %
            % % convert to a MATLAB struct
            % us = obj2struct(us)
            %
            arguments
                us UltrasoundSystem
            end

            % squash warnings
            wmsg = ["MATLAB:structOnObject", "QUPS:UltrasoundSystem:syntaxDeprecated"];
            S = warning(); % warning state
            for w = wmsg, warning('off', w); end
            
            s = struct(us); % convert us to a struct
            if isfield(s,'xdc'), s = rmfield(s,{'tx','rx'}); end % remove duplicate information
            s = rmfield(s, {'fc', 'sequence', 'tmp_folder'}); % inferred / duplicate / private

            % convert all sub-classes
            for f = intersect(string(fieldnames(s))', ["tx","rx","xdc","seq","scan"])
                if isobject(s.(f))
                    s.(f) = obj2struct(s.(f)); 
                end
            end
            s.class = class(us); % append class info

            % restore warnings
            warning(S);
        end    
    end
    
    % overloading methods
    methods(Access=protected, Hidden)
        % copy 
        function other = copyElement(self)
            arguments, self (1,1) UltrasoundSystem, end % always scalar in            
            other = copyElement@matlab.mixin.Copyable(self); % shallow copy handle
            other.tmp_folder = tempname(); % new temp dir
            if exist(self.tmp_folder, 'dir')
                copyfile(self.tmp_folder, other.tmp_folder); % copy binaries over
            else
                mkdir(other.tmp_folder); % create empty
            end
            addpath(other.tmp_folder); % add to path
            % disp("[DEBUG]: Adding path " + other.tmp_folder + " with " + (numel(dir(other.tmp_folder)) - 2) + " files.");
        end
    end

    % display
    methods
        function h = plot(us, ax, im_args)
            % PLOT - Plot the geometry of the UltrasoundSystem 
            %
            % h = PLOT(us) plots the UltrasoundSystem us by plotting 
            % the transmit and receive transducer(s), the imaging pixels,
            % and representation of the transmit sequence on the same axes.
            %
            % h = PLOT(us, ax) plots on the axes ax. The default is the
            % current axes returned by gca.
            %
            % h = PLOT(..., Name, Value, ...) passes following arguments to
            % the call to plot for each of the objects i.e. to the plot
            % function for us.scan, us.seq, and us.xdc or us.rx and us.tx
            % if using separate transmit and receive transducers.
            %
            % Example:
            % % Create a default system using a focused transmit
            % xdc = TransducerArray();
            % us = UltrasoundSystem('xdc', xdc);
            % us.seq.type = 'FC';
            % us.seq.focus = [0;0;30e-3] + ...
            % linspace(-xdc.numel/4, xdc.numel/4, xdc.numel/2+1) .* [xdc.pitch;0;0];
            %
            % % plot it
            % figure;
            % h = plot(us, 'LineWidth', 2);
            %
            % See also PLOT ULTRASOUNDSYSTEM/XDC ULTRASOUNDSYSTEM/SCAN ULTRASOUNDSYSTEM/SEQUENCE

            arguments
                us UltrasoundSystem
                ax (1,1) matlab.graphics.axis.Axes = gca;
            end
            arguments
                im_args.?matlab.graphics.primitive.Line
            end

            hstate = ax.NextPlot;
            hold(ax, 'on');
            set(ax, 'ydir', 'reverse');
            title(ax, 'Geometry');
            plargs = struct2nvpair(im_args);
            hps = plot(us.scan, '.', 'Color', [1.0 0.75 0.75], 'DisplayName', 'Grid', plargs{:}); % the imaging points
            if us.tx == us.rx
                hxdc = plot(us.xdc, ax, 'Color', [1 0 1], 'DisplayName', 'Elements', plargs{:}); % elements
            else
                htx = plot(us.tx, ax, 'b', 'DisplayName', 'Transmit Elements', plargs{:}); % tx elements
                hrx = plot(us.rx, ax, 'r', 'DisplayName', 'Receive Elements', plargs{:}); % rx elements
                hxdc = [htx, hrx];
            end
            rq = max([range([hps.XData]), range([hps.YData]), range([hps.ZData])]); % scale quivers by range of the image
            % rq = max(rq, range(xlim(ax))); % include x?
            switch us.seq.type % show the transmit sequence
                case 'PW', hseq = plot(us.seq, ax, rq/2, 'k.', 'DisplayName', 'Tx Sequence', plargs{:}); % scale the vectors for the plot
                otherwise, hseq = plot(us.seq, ax,       'k.', 'DisplayName', 'Tx Sequence', plargs{:}); % plot focal points, if they exist
            end
            h = [hxdc, hps, hseq];
            legend(ax, h);
            grid(ax, 'on'); 
            grid(ax, 'minor')
            ax.NextPlot = hstate;
        end
    end

    % object display
    methods(Access=protected, Hidden)
        function propgrp = getPropertyGroups(us)
            if ~isscalar(us)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(us);
            else
                p = string(properties(us)); % get public properties
                p(end+1) = "tmp_folder";% + hidden temp folder
                if us.rx == us.tx
                    p(ismember(p, ["rx","tx"])) = []; % only display 'xdc'
                else 
                    p(ismember(p,    "xdc"   )) = []; % display 'tx' and 'rx'
                end
                propgrp = matlab.mixin.util.PropertyGroup(p);
            end
        end
    end

    % property modification
    methods
        function us = scale(us, kwargs)
            % SCALE - Scale the units of the system
            %
            % us = SCALE(us, 'dist', dscale) scales the values of 
            % distance by dscale.
            % 
            % us = SCALE(us, 'time', tscale) scales the values of
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
                us (1,1) UltrasoundSystem
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end

            us = copy(us);
            args = struct2nvpair(kwargs);
            
            % Scan: convert distance only
            if isfield(kwargs, 'dist')
                us.scan = scale(us.scan, 'dist', kwargs.dist);
            end
            % Sampling: convert time only
            if isfield(kwargs, 'time')
                us.fs = us.fs / kwargs.time;
            end

            % Transducer: convert distance and/or time/freq
            if us.tx == us.rx
                us.xdc = scale(us.xdc, args{:}); % convert in distance and time
            else
                us.tx = scale(us.tx, args{:}); % convert in distance and time
                us.rx = scale(us.rx, args{:}); % convert in distance and time
            end

            % Sequence
            us.seq = scale(us.seq, args{:}); % convert in distance and time
        end
    end

    % Modified Green's function based direct computations
    methods
        function [chd, wv] = greens(us, scat, element_subdivisions, kwargs)
            % GREENS - Simulate ChannelData via a shifted Green's function.
            % 
            % chd = GREENS(us, scat) simulates the response of the 
            % Scatterers scat from the UltrasoundSystem us and returns 
            % the corresponding ChannelData chd. It first computes the full
            % synthetic aperture data using a simple Green's function 
            % kernel applied to all sub-elements and all point scatterers 
            % and then applies the transmit Sequence via focusTx.
            %
            % chd = GREENS(us, scat, element_subdivisions) uses the 
            % length 2 array element_subdivisions to specifiy the 
            % subdivision of each element into a grid of sub-apertures 
            % in the integration. This argument is passed to FieldII to 
            % construct the subelements. FieldII must be on the path.
            % Setting the subdvisions to [1,1] avoids this behaviour. The
            % Default is [1,1].
            %
            % chd = GREENS(..., 'R0', R0) sets the minimum distance for
            % propagation loss. Setting R0 to 0 disables propagation loss.
            % The default is us.lambda.
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
            % [...] = GREENS(..., 'fsk', fsk ...) samples and convolves the
            % transmit sequence pulse and transducer impulse response
            % functions at a sampling frequency of fsk. The default is
            % us.fs.
            %
            % Note: The precision of the computation is determined by the
            % precision of the transmit pulse. Setting 
            %       `us.fs = single(us.fs);` 
            % before calling greens can dramatically accelerate computation
            % time on consumer GPUs which may have a performance ratio of
            % 32:1 or 64:1 for single:double precision.
            % 
            % About: This function provides a fast, low-fidelity simulation
            % routine based on assuming isotropic radiation from a point
            % source and simply applying a Green's function approximation.
            % 
            % 
            % 
            % Example:
            % % Simulate some data
            % us = UltrasoundSystem(); % get a default system
            % us.fs = single(us.fs); % use single precision for speed
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % chd = greens(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            %
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_MULTI ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/KSPACEFIRSTORDER 
            
            arguments
                us (1,1) UltrasoundSystem
                scat Scatterers
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
            end
            arguments
                kwargs.device (1,1) {mustBeInteger} = -1 * (logical(gpuDeviceCount()) || (exist('oclDeviceCount','file') && logical(oclDeviceCount())))
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["double","parallel.Pool","parallel.Cluster"])} = gcp('nocreate');
                kwargs.tall (1,1) logical = false; % whether to use a tall type
                kwargs.bsize (1,1) {mustBeInteger, mustBePositive} = max([1 us.seq.numPulse],[],'omitnan'); % number of simulataneous scatterers
                kwargs.verbose (1,1) logical = false;
                kwargs.fsk (1,1) {mustBePositive} = us.fs % input kernel sampling frequency
                kwargs.R0 (1,1) double {mustBeNonnegative} = max(us.lambda) % minimum distance for divide by 0
            end
            
            % get the centers of all the sub-elements
            if all(element_subdivisions == 1) % no sub-elements
                pv = us.tx.positions();
                pn = us.rx.positions();
            else % TODO: use xdc.patches
                pv = us.tx.getBaryCenters(element_subdivisions);
                pn = us.rx.getBaryCenters(element_subdivisions);
            end

            % cast dimensions to compute in parallel
            pn = permute(pn, [6,1,2,4,3,5]); % 1 x 3 x N x 1 x En x 1
            pv = permute(pv, [6,1,4,2,5,3]); % 1 x 3 x 1 x M x  1 x Em

            % get the boundaries of the transducer
            txb = us.tx.bounds();
            rxb = us.rx.bounds();

            % select each corner
            j = 1+logical([bitand(uint16([1; 2; 4]), uint16(0:7))]);
            txb = swapdim(sel(txb, j, 2),2,3); % 3 x 1 x 8
            rxb = swapdim(sel(rxb, j, 2),2,3); % 3 x 1 x 8
            
            % get maximum necessary time sample (use manhattan distance and
            % sound speed to give an upper bound)
            maxdist = @(p) max(vecnorm(p,2,1), [], 'all');
            mindist = @(p) min(vecnorm(p,2,1), [], 'all');
            taumax = arrayfun(@(scat)(maxdist(scat.pos - txb) + maxdist(scat.pos - rxb) + vecnorm(range(txb,3)) + vecnorm(range(rxb,3))) ./ min(scat.c0), shiftdim(scat(:),-3));
            taumin = arrayfun(@(scat)(mindist(scat.pos - txb) + mindist(scat.pos - rxb) - vecnorm(range(txb,3)) - vecnorm(range(rxb,3))) ./ max(scat.c0), shiftdim(scat(:),-3));
            
            % Directly convolve the Waveform objects to get the final
            % convolved kernel
            wv = conv(us.rx.impulse, ...
                conv(us.tx.impulse, us.seq.pulse, gather(kwargs.fsk)), ...
                gather(kwargs.fsk)); % transmit waveform, convolved at US frequency
            wv.fs = kwargs.fsk;
            kern = wv.samples;

            % choose device
            use_dev = kwargs.device && ismember(kwargs.interp, ["nearest", "linear", "cubic", "lanczos3"]);
            use_gdev = use_dev && exist('greens.ptx', 'file') && gpuDeviceCount(); % use the GPU kernel
            use_odev = use_dev && exist('oclDeviceCount','file') && oclDeviceCount(); % use the OpenCL kernel
            if use_odev 
                if kwargs.device > 0, oclDevice(kwargs.device); end % select device
                dev = oclDevice(); % reference
                use_odev = isscalar(dev) ...
                    && (isUnderlyingType(kern, "single") ...
                    || (isUnderlyingType(kern, "double") && dev.SupportsDouble) ...
                    || (isUnderlyingType(kern, "half"  ) && dev.SupportsHalf));
            end
                
            F = numel(scat);
            % if kwargs.verbose, hw = waitbar(0); end

            for f = F:-1:1 % for each Scatterers
            % get minimum/maximum sample times
            tmin = taumin(f) + wv.t0 - wv.duration; % duration added as a quick hack
            tmax = taumax(f) + wv.tend;

            % create time vector (T x 1)
            % this formulation is guaranteed to pass through t == 0
            n0 = floor(tmin * us.fs);
            ne = ceil(tmax * us.fs);
            t = (n0 : ne)';

            % pre-allocate output
            [T, N, M, E] = deal(numel(t), us.rx.numel, us.tx.numel, prod(element_subdivisions));
            x   = complex(zeros([1 T N M 0], 'like', kern)); % set size/type
            x(:,:,:,:,1) = 0; % pre-allocate (once, not twice)

            % splice
            c0  = scat(f).c0;
            ps = scat(f).pos; % 3 x S
            as = scat(f).amp; % 1 x S
            fso = us.fs; % output channel data sampling frequency
            if use_gdev, ps = gpuArray(ps); end

            % compute the maximum distance for each scatterer
            % sort points by geometric mean of minimum/maximum distance
            [rminrx, rmaxrx, rmintx, rmaxtx] = deal(+inf, -inf, +inf, -inf);
            for p = us.tx.positions, rminrx = min(rminrx, vecnorm(ps - p, 2, 1)); end % minimum scat time
            for p = us.tx.positions, rmaxrx = max(rmaxrx, vecnorm(ps - p, 2, 1)); end % maximum scat time
            for p = us.rx.positions, rmintx = min(rmintx, vecnorm(ps - p, 2, 1)); end % minimum scat time
            for p = us.rx.positions, rmaxtx = max(rmaxtx, vecnorm(ps - p, 2, 1)); end % maximum scat time
            [rmin, rmax] = deal(rmintx + rminrx, rmaxtx + rmaxrx);
            [~, i] = sort((rmin .* rmax)); % get sorting

            % sort points by geometric mean of minimum/maximum distance
            [ps, as, rmin, rmax] = dealfun(@(x)sub(x,i,2), ps,as,rmin,rmax); % apply to all position variables

            if use_gdev || use_odev
                % function to determine type
                isftype = @(x,T) strcmp(class(x), T) || any(arrayfun(@(c)isa(x,c),["tall", "gpuArray"])) && strcmp(classUnderlying(x), T);

                % determine the data type
                if     isftype(kern, 'double'), typ = "double"; prc = 64; suffix = "" ; cfun = @double;
                elseif isftype(kern, 'single'), typ = "single"; prc = 32; suffix = "f"; cfun = @single;
                elseif isftype(kern, 'halfT' ), typ = "halfT" ; prc = 16; suffix = "h"; cfun = @(x) alias(halfT(x));
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
                [   x, ps, as, pn, pv, kn, t0k, t0x, fso, cinv_] = dealfun(cfun, ...
                    x, ps, as, pn, pv, kern, t(1)/fso, wv.t0, fso, 1/c0 ...
                    );

                % re-map sizing
                [QI, QS, QT, QN, QM] = deal(scat(f).numScat, length(t), length(kern), N, M);

                % get the index bounds for the output time axis
                sb = ([rmin; rmax] ./ c0 + t0x - t0k) * fso + [0; QT];

                % grab the kernel reference
                if use_gdev
                    k = parallel.gpu.CUDAKernel('greens.ptx', 'greens.cu', 'greens' + suffix);
                    k.setConstantMemory( 'QUPS_S', uint64(QS) ); % always set S
                    try k.setConstantMemory('QUPS_T', uint64(QT), ...
                            'QUPS_N', uint64(QN), 'QUPS_M', uint64(QM) ...
                            ); %#ok<TRYNC> % already set by const compiler
                    end 
                    try k.setConstantMemory('QUPS_I', uint64(QI)); end %#ok<TRYNC> % might be already set by const compiler
                elseif use_odev
                    k = oclKernel('greens.cl');
                    k.macros = [k.macros, "QUPS_"+["S","T","N","M","I"]+"="+uint64([QS,QT,QN,QM,QI])]; % set constants
                    k.macros(end+1) = "QUPS_PRECISION="+prc;
                    k.defineTypes(repmat(typ,[1,3]), ["V","T","U"]); % time / data / position

                    % enforce complexity
                    [x, as, kn] = dealfun(@complex, x, as, kn);

                    % expand positions 3D -> 4D
                    [ps(4,:), pn(1,4,:), pv(1,4,:)] = deal(0);

                end                
                
                % set kernel sizing
                k.ThreadBlockSize(1) = min([k.MaxThreadsPerBlock, 32]); % smaller is better
                k.GridSize = [ceil(QS ./ k.ThreadBlockSize(1)), N, M];

                % get the computational bounds
                sblk = (0 : k.GridSize(1) - 1)' * k.ThreadBlockSize(1); % block starting time index
                blocks = sblk <= [-inf, sb(2,:), inf]; % point at which we are in bounds
                blocks = blocks(:,2:end) - blocks(:,1:end-1); % detect change
                sbk = cellfun(@(x) find(x,1,'first'), num2cell(blocks,2)); % get transition point
                blocks = (sblk + k.ThreadBlockSize(1)) < [-inf, sb(1,:), inf]; % point at which we are in bounds
                blocks = blocks(:,2:end) - blocks(:,1:end-1); % detect change 
                ebk = cellfun(@(x) find(x,1,'last'), num2cell(blocks,2)); % get transition point
                
                % call the kernel
                if kwargs.verbose, tt = tic; fprintf("Computing for " + gather(sum((ebk - sbk) + 1)) + " blocks."); end
                x = k.feval(x, ps, as, pn, pv, kn, sb, gather([sbk,ebk]'-1), [t0k, t0x, fso, wv.fs/fso, cinv_, kwargs.R0], [E,E], flagnum);
                
            else % operate in native MATLAB
                % parallel environment
                penv = kwargs.parenv;
                if isempty(penv), penv = 0; end % no parallel

                % make time in dim 2, scats in dim 1
                [ps, as] = deal(ps', as'); % transpose S x 3, S x 1
                tvec = reshape(t,1,[]); % 1 x T full time vector
                kern_ = reshape(kern,1,[]); % 1 x T' signal kernel
                K = length(kern_); % kernel length in indices

                % timing info
                t0_ = gather(wv.t0); % 1 x 1
                fs_ = us.fs;
                fsr = wv.fs / fs_; % sampling frequency ratio

                % TODO: set data types | cast to GPU/tall? | reset GPU?
                if use_gdev, gpuDevice(kwargs.device); end
                if use_gdev && ~kwargs.tall
                    [ps, pn, pv, tvec, kern_] = dealfun(@gpuArray, ...
                     ps, pn, pv, tvec, kern_ ...
                        );
                elseif kwargs.tall
                    [pn, pv] = dealfun(@tall, pn, pv);
                end

                % parallel worker size
                if isempty(p) % no pool -> 1
                    P = 1; 
                elseif isnumeric(penv) 
                    if penv, P = penv; else, P = 1; end % 0 -> 1
                else % pool
                    P = max(1,penv.NumWorkers); % ???
                end 
                
                % get block sizes: if less than M, go 1 scat at a time
                [B, S] = deal(floor(kwargs.bsize), scat(f).numScat);
                Bm = max(1,min(floor(B/P),floor(M/P))); % simulataneous transmits
                Bs = max( 1, floor(B / Bm / P) ); % simulataneous scats
                svec = num2cell( (1:Bs)'+(0:Bs:S-1), 1);
                mvec = num2cell( (1:Bm)'+(0:Bm:M-1), 1);
                svec{end} = svec{end}(svec{end} <= S);
                mvec{end} = mvec{end}(mvec{end} <= M);

                % splice
                [R0, terp, kvb] = deal(kwargs.R0, kwargs.interp, kwargs.verbose);
                if kvb
                    disp("Partitioning into " ...
                        +numel(svec)+ " blocks of " ...
                        +Bs+" scatterer(s) for " ...
                        +numel(mvec)+" sets of " ...
                        +Bm+" elements across "...
                        +P+" worker(s)." ...
                        );
                    tt = tic;
                end

                x = cell(size(mvec));% zeros(xsz,'like',kern_); % init x = 0;
                parfor(mv = 1:numel(mvec), penv)  % each transmit
                    m = mvec{mv}; % vector of tx indices
                    M_ = numel(m); % number of transmit this block

                    % pre-allocate this transmit
                    x_ = zeros([1,T,N,M_],'like',kern_);
                    x{mv} = 0;
                    
                    for sv = 1:numel(svec) % each set of scatterers
                    s = svec{sv}; % vector of scat indices
                    
                    for em = 1:E
                    for en = 1:E
                        % toc; % DEBUG
                        % tic; % DEBUG

                        % compute time delays
                        % TODO: do this via ray-path propagation through a medium
                        % S x 1 x N x /M x 1 x 1
                        r_rx = vecnorm(sub(ps,s,1) - sub(pn,     en,      5 ),2,2);
                        r_tx = vecnorm(sub(ps,s,1) - sub(pv, {m, em}, [4, 6]),2,2);
                        tau_rx = (r_rx ./ c0); % S x 1 x N x 1 x 1 x 1
                        tau_tx = (r_tx ./ c0); % S x 1 x 1 x M x 1 x 1

                        % propagation loss
                        if R0 % min distance - no loss if 0
                            r_rx = max(r_rx, R0); % avoid div by 0
                            r_tx = max(r_tx, R0); % avoid div by 0
                            att  = 1; % R0^-2; % scale by max energy?
                        else
                            [r_rx, r_tx, att] = deal(1); % no propagation loss
                        end

                        % compute the attenuation (S x 1 x [1|N] x [1|M] x 1 x 1)
                        att = att .* sub(as,s,1) ./ (r_rx .* r_tx); % propagation attenuation

                        % convert time delays to indices
                        % S x T x N x M x 1 x 1
                        tau_tx = tau_tx * fs_;
                        tau_rx = tau_rx * fs_;
                        t0     = t0_    * fs_;

                        % Get min/max samplipng times for this set
                        % compute only for this range and sum as we go
                        % tau = (tau_tx + tau_rx + t0); % S x 1 x N x M x 1 x 1
                        tmax = t0 + max(tau_tx + tau_rx,[],'all');
                        tmin = t0 + min(tau_tx + tau_rx,[],'all');
                        it = (tmin - K-1) <= tvec & tvec <= (tmax + K+1); %(1 x T')

                        % S x T x N x M x 1 x 1
                        if any(cellfun(@istall, {tau_tx, tau_rx, tvec, kern_}))
                            % compute as a tall array  % S | T x N x M x 1 x 1
                            it = gather(it);
                            t1 = fsr * (tvec(it) - tau_tx - t0); % create tall ND-array
                            t2 = - fsr * tau_rx;
                            s_ = matlab.tall.reduce( ...
                                @(t1, t2, att) wsinterpd2( ...
                                kern_, t1, t2, 2, att ./ fsr, 1, terp, 0 ...
                                ), ... map function
                                @(x)sum(x,1,'omitnan','native'), ... reduce function over scatterers
                                t1, t2, att...
                                );
                        else
                            % compute natively
                            % (1 x T' x N X M)
                            s_ = nan2zero(wsinterpd2(kern_, fsr * (tvec(it) - tau_tx - t0), - fsr * tau_rx, 2, att ./ fsr, 1, terp, 0));
                        end

                        % add contribution (1 x T x N x M')
                        x_(1,it,:,:) = x_(1,it,:,:) + s_; % (1 x T' x N x M')
                    end
                    end
                    end
                    x{mv} = x{mv} + x_; % accumulate
                    if kvb, fprintf("."); end % update
                end
                
                % unpack
                if isscalar(x), x = x{1}; else, x = cat(4, x{:}); end

                % compute for tall arrays
                if istall(x), x = gather(x); end

                % move back to GPU if requested
                if use_gdev, x = gpuArray(gather(x)); end
            end
            if kwargs.verbose, disp("Done!"); toc(tt); end
            
            % make a channel data object (T x N x M)
            x = reshape(x, size(x,[2:ndims(x) 1])); % same as shiftdim(x,1), but without memory copies on GPU
            chd(f) = ChannelData('t0', sub(t,1,1) ./ us.fs, 'fs', us.fs, 'data', x);

            % truncate the data if possible
            iszero = all(chd(f).data == 0, 2:ndims(chd(f).data)); % true if 0 for all tx/rx/targs
            n0 = gather(find(cumsum(~iszero, 'forward'), 1, 'first'));
            T_ = gather(find(cumsum(~iszero, 'reverse'), 1, 'last' ));
            chd(f) = subD(chd(f), n0:T_, chd(f).tdim);

            % synthesize linearly
            [chd(f)] = us.focusTx(chd(f), us.seq, 'interp', kwargs.interp, 'bsize', kwargs.bsize, 'verbose', kwargs.verbose);
            end

            % combine all frames
            chd = join(chd, 4);
        end
    end

    % USTB interop
    methods
        function [uchannel_data, uscan] = QUPS2USTB(us, chd, fmod)
            % QUPS2USTB - Create a USTB channel data object
            %
            % uchannel_data = QUPS2USTB(us, chd) creates a USTB
            % compatible uff.channel_data object from the UltrasoundSystem
            % us and ChannelData chd. USTB must be on the path.
            %
            % uchannel_data = QUPS2USTB(..., fmod) sets the modulation
            % frequency to fmod. The default is 0.
            %
            % [uchannel_data, uscan] = QUPS2USTB(...) additionally returns
            % a USTB compatible uff.scan.
            %
            % Example:
            % % Create an UltrasoundSystem and ChannelData
            % xdc = TransducerArray.L12_3v();
            % seq = SequenceRadial('type','PW','angles',-10 : 10);
            % us = UltrasoundSystem('seq', seq, 'xdc', xdc);
            % chd = ChannelData('data',zeros([128,xdc.numel,seq.numPulse]));
            % 
            % % Convert to USTB
            % uchannel_data = QUPS2USTB(us, chd);
            % 
            % See also ULTRASOUNDSYSTEM.UFF CHANNELDATA.QUPS2USTB
            arguments
                us (1,1) UltrasoundSystem
                chd (1,1) ChannelData
                fmod (1,1) {mustBeReal, mustBeFinite} = 0
            end
            uchannel_data = QUPS2USTB(chd, us.seq, us.xdc, fmod);
            uscan         = QUPS2USTB(us.scan);
        end
    end

    methods(Static)
        function [us, chd] = UFF(uchannel_data, uscan)
            % UFF - Create an UltrasoundSystem and ChannelData from a uff.channel_data object
            %
            % [us, chd] = UFF(uchannel_data) creates an
            % UltrasoundSystem us and ChannelData chd from the
            % uff.channel_data object uchannel_data.
            %
            % [us, chd] = UFF(uchannel_data, uscan) additionally sets the
            % Scan us.scan from the uff.scan uscan.
            %
            % See also ChannelData.UFF Sequence.UFF Transducer.UFF Scan.UFF

            arguments
                uchannel_data (1,1) uff.channel_data
                uscan uff.scan {mustBeScalarOrEmpty} = uff.scan.empty
            end
            fs = uchannel_data.sampling_frequency; % sampling frequency
            seq = Sequence.UFF(uchannel_data.sequence, uchannel_data.sound_speed); % Sequence
            xdc = Transducer.UFF(uchannel_data.probe); % Transducer
            us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'fs', fs);
            if ~isempty(uscan), us.scan = Scan.UFF(uscan); end
            if nargout >= 2, chd = ChannelData.UFF(uchannel_data, us.seq, us.xdc); end % data
        end
    end

    % Verasonics interop
    methods(Static)
        function [us, chd] = Verasonics(Trans, TX, TW, kwargs)
            % Verasonics - Import a Verasonics Vantage system and data
            % 
            % us = Verasonics(Trans, TX, TW) creates an UltrasoundSystem us
            % with a Sequence and Transducer defined from the Verasonics 
            % structs `Trans`, `TX`, and `TW`.
            % 
            % us = Verasonics(Trans, TX) or 
            % us = Verasonics(Trans, TX, struct.empty) uses a
            % Waveform.Delta as the excitation pulse function.
            % 
            % us = Verasonics(..., 'c0', c0) additionally specifies the
            % sound speed in m/s. This should match the Verasonics property
            % `Resource.Parameters.speedOfSound`. The default is 1540.
            % 
            % us = Verasonics(..., 'PData', PData) additionally defines a
            % Scan from the Verasonics struct `PData`.
            % 
            % [us, chd] = Verasonics(..., 'Receive', Receive, 'RcvData', RcvData)
            % additionally returns a ChannelData chd defined by the
            % Verasonics `Receive` and `RcvData` structs.
            % 
            % If `RcvData` is omitted or empty, chd will contain data of
            % size (time x acquisitions x receivers x frames x 0).
            % 
            % If `Receive` is omitted or empty, chd will contain data of
            % size (0 x transmits x receivers x 0 x 0) and assumes a
            % supported sampling frequency of approximately 4 * us.xdc.fc.
            % 
            % Example:
            % 
            % load('my_VSX_data.mat')
            % c0 = Resource.Parameters.speedOfSound;
            % us = UltrasoundSystem.Verasonics(Trans, TX, TW, 'c0', c0);
            %
            % d = load('my_other_data.mat'); % load into a struct
            % [us, chd] = UltrasoundSystem.Verasonics(...
            %     d.Trans, d.TX, d.TW ... required
            %     , 'c0', d.Resource.Parameters.speedOfSound ... Matching sound speed
            %     , 'Receive', d.Receive, 'RcvData', d.RcvData ... Import ChannelData
            %     , 'PData', d.PData ... Import Scan
            % );
            % chd = singleT(chd); % int16 -> single
            % 
            % See also CHANNELDATA.VERASONICS TRANSDUCER.VERASONICS
            % SEQUENCE.VERASONICS SCAN.VERASONICS
            arguments
                Trans (1,1) struct
                TX struct
                TW (1,:) struct = struct.empty
                kwargs.PData struct {mustBeScalarOrEmpty} = struct.empty
                kwargs.RcvData cell = {}
                kwargs.Receive struct = struct.empty
                kwargs.c0 (1,1) double = 1540
            end
            
            % wavelength
            lbda = kwargs.c0 / Trans.frequency / 1e6; 
            
            % create US
            xdc = Transducer.Verasonics(Trans, kwargs.c0);
            if numel(TW) > 1, TW = TW(unique([TX.waveform])); end
            [seq, t0q] = Sequence.Verasonics(TX, Trans, TW, 'c0', kwargs.c0, 'xdc', xdc);
            if isempty(kwargs.PData), scan_args = {};
            else, scan_args = {'scan', Scan.Verasonics(kwargs.PData, lbda)};
            end
            us = UltrasoundSystem('seq', seq, 'xdc', xdc, scan_args{:}, 'recompile',false);
            
            
            % construct ChannelData
            if ~isempty(kwargs.Receive)
                chd = ChannelData.Verasonics(kwargs.RcvData, kwargs.Receive, Trans);
                us.fs = chd.fs;
            else
                fs = 250./(4:100); % supported frequencies
                fs = fs(argmin(abs(4*Trans.frequency - fs))); % assume 200BW, closest frequency
                sz = [0, numel(TX), Trans.numelements, 0, 0]; % output data sizing
                chd = ChannelData('data', zeros(sz,'int16'), 'order','TMN', 'fs', fs);
            end
            
            % refine t0
            if isfield(Trans, 'lensCorrection')
                t0x = 2*Trans.lensCorrection; else, t0x = 0;
            end
            if isfield(TW,'peak') && ~isempty(TW.peak)
                t0p = TW.peak; else, t0p = 0;
            end
            chd.t0 = chd.t0 + swapdim(t0q,2,chd.mdim) - (t0p + t0x) / (Trans.frequency * 1e6);
        end
    end
    
    % Fullwave calls (Hidden)
    methods(Hidden)
        function conf = fullwaveConf(us, medium, cgrd, kwargs)
            % FULLWAVECONF - Generate a Fullwave simulation configuration
            %
            % conf = FULLWAVECONF(us, medium) creates a simulation
            % configuration struct to be used with fullwaveJob to simulate 
            % the response from the Medium medium using the 
            % UltrasoundSystem us.
            % 
            % conf = FULLWAVECONF(us, medium, cgrd) uses the  
            % ScanCartesian cgrd as the simulation region. The default is
            % us.scan.
            %
            % conf = FULLWAVECONF(..., 'f0', f0) uses a reference frequency
            % of f0 to configure the simulation. The default is us.tx.fc.
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
            % Note: This functionality is not publicly supported.
            % 
            % Example:
            % grd = ScanCartesian('x', 1e-3*[-35, 35], 'z', 1e-3*[-12, 27]);
            % xdc = TransducerConvex.C5_2v();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', grd, 'xdc', xdc, 'seq', seq);
            % [grd.dx, grd.dz] = deal(us.lambda / 4);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(grd.size), 1000*ones(grd.size));
            % pg = grd.positions();
            % rho(argmin(vecnorm(pg - [0 0 20e-3]',2,1))) = 1000*2; % double the density
            % med = Medium.Sampled(grd, c, rho);
            % 
            % % Create a configuration struct
            % us.scan = grd; % set the simulation grid
            % conf = fullwaveConf(us, med), % configure the sim
            % 
            % See also ULTRASOUNDSYSTEM/FULLWAVEJOB

            arguments
                us (1,1) UltrasoundSystem
                medium Medium
                cgrd ScanCartesian = us.scan
                kwargs.f0 (1,1) {mustBeNumeric} = us.tx.fc; % center frequency of the transmit / simulation
                kwargs.CFL_max (1,1) {mustBeReal, mustBePositive} = 0.3 % maximum CFL
                kwargs.txdel (1,1) string {mustBeMember(kwargs.txdel, ["discrete", "continuous", "interpolate"])} = 'interpolate';
            end            

            %% Configuration variables

            % basic vars
            c0       = medium.c0;       % speed of sound (m/s)
            omega0   = 2*pi*kwargs.f0;  % center radian frequency of transmitted wave
            dur      = diff(cgrd.zb)*2.3/c0; % duration of simulation (s) TODO: make this more general, or an input?

            % determine other grid vars
            dX   = min(abs([cgrd.dx, cgrd.dy, cgrd.dz]), [], 'omitnan');  % limit of spatial step size - will be the same in both dimensions
            fs_  = us.fs;         % data sampling frequency
            cfl0 = (c0*(1/fs_)/dX); % cfl at requested frequency
            modT = ceil(cfl0 / kwargs.CFL_max); % scaling for desired cfl
            dT   = (1/fs_)/modT;    % simulation sampling interval
            cfl  = c0 * dT / dX;    % Courant-Friedrichs-Levi condition
            ppw  = c0 / kwargs.f0 / dX; % points per wavelength - determines grid spacing
            dur  = ceil(dur / dT) * dT; % snap duration to temporal grid
            

            % DEBUG 1: this must be consistent with itself!!!
            % cfl,c0,ppw,omega0,dX,dY all go together!
            % at the end: conf.sim  = {c0,omega0,dur,ppw,cfl,maps2,xdcfw,nTic,modT};

            %% Define the Transducer(s)
            xdcfw = us.tx.getFullwaveTransducer(cgrd);
            if us.tx == us.rx
                rxfw = xdcfw; else, rxfw = us.rx.getFullwaveTransducer(cgrd);
            end

            %% Define the Transmit Delays

            % get transmit delays
            tau_tx = -us.seq.delays(us.tx); % M x V

            % forward to all subelement delays, synchronized
            tau_tx_pix = zeros([xdcfw.nInPx, size(tau_tx,2)]);
            i = logical(xdcfw.incoords(:,4)); % find non-zero indices - they map to an element
            tau_tx_pix(i,:) = tau_tx(xdcfw.incoords(i,4),:); % apply synchronized delays

            
            %% Define the Transmit Apodization

            % get apodization per element
            apod = us.seq.apodization(us.tx);

            % map to sub-elements
            tx_apod = zeros(size(tau_tx_pix)); % pre-allocate
            i = logical(xdcfw.incoords(:,4)); % find non-zero indices - they map to an element
            tx_apod(i,:) = apod(xdcfw.incoords(i,4),:); % apply synchronized delays

            %% Define the Transmit Pulse
            t0_xdc = min(min(tau_tx_pix(i,:),[],1)); % reference true start time (1 x M)
            tau_tx_pix = (tau_tx_pix - t0_xdc); % 0-base the delays for the sim (I x M)
            tau_tx_max = ceil(max(tau_tx_pix, [],'all') / dT) .* dT; % maxmium additional delay
            wv_tx = conv(us.tx.impulse, us.seq.pulse, 10/dT); % get the waveform transmitted into the medium
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
                    wv_tx.fs = 10*fs_;
                    t_up = wv_tx.time;
                    icmat = interp1(t_up, wv_tx.sample(t_up), t(:) - tau_tx_pix, 'spline', 0); 
            end

            % apply transmit apodization
            icmat = icmat .* shiftdim(tx_apod,-1); % (T' x nInPx x nTx)

            %% Define the Medium
            % get maps in some order
            maps = medium.getFullwaveMap(cgrd);
            assert(isscalar(cgrd.y), 'Unable to simulate in 3D (y-axis is not scalar).');

            % switch to X x Z x Y(=1)
            ord = arrayfun(@(c) find(cgrd.order == c), 'XZY');
            for f = string(fieldnames(maps))', maps.(f) = permute(maps.(f), ord); end

            %% Store the configuration variables
            conf.sim = {c0,omega0,dur,ppw,cfl,maps,rxfw,nTic,modT}; % args for full write function
            conf.tx  = real(icmat);  % transmit data
            conf.t0  = shiftdim(t0_xdc, -1); % time delay offset (1 x 1 x M)
            conf.f0  = kwargs.f0;     % simulation frequency
            conf.fs  = fs_;     % sampling frequency
            conf.tstart = t(1); % signal start time (0 is the peak)
            conf.outmap = rxfw.outcoords(:,4); % 4th column maps pixels to elements
        end
        
        function [chd, conf] = fullwaveSim(us, medium, grid, kwargs)
            % FULLWAVESIM - Simulate channel data via Fullwave
            %
            % chd = FULLWAVESIM(us, medium, grid) simulates the Medium 
            % medium on the simulation grid grid and returns a ChannelData
            % object chd. The simulation grid should be large and fine 
            % enough such that all elements of the Transducer can be placed
            % onto the grid.
            %
            % chd = FULLWAVESIM(us, medium) uses us.scan as the
            % simulation grid.
            %
            % chd = FULLWAESIM(..., 'simdir', dir) uses the directory dir
            % to store temporary simulation files. The default is a folder
            % in the working directory.
            %
            % chd = FULLWAESIM(..., 'parenv', clu) uses the
            % parallel.Cluster clu to compute each pulse. The simulation
            % directory must be accesible by the parallel.Cluster clu.
            % 
            % [chd, conf] = FULLWAVESIM(...) also returns the configuration
            % structure used to launch the simulation.
            %
            % Note: This functionality is not publicly supported.
            % 
            % Example:
            % 
            % % Setup a system
            % grid = ScanCartesian('x', 1e-3*[-35, 35], 'z', 1e-3*[-12, 27]);
            % xdc = TransducerConvex.C5_2v();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', grid, 'xdc', xdc, 'seq', seq);
            % [grid.dx, grid.dz] = deal(us.lambda / 4);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(grid.size), 1000*ones(grid.size));
            % pg = grid.positions();
            % rho(argmin(vecnorm(pg - [0 0 20e-3]',2,1))) = 1000*2; % double the density
            % med = Medium.Sampled(grid, c, rho);
            % 
            % % Simulate the ChannelData
            % chd = fullwaveSim(us, med);
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(hilbert(chd));
            % dbr echo 60;
            % 
            % See also ULTRASOUNDSYSTEM.KSPACEFIRSTORDER
            
            arguments
                us (1,1) UltrasoundSystem
                medium Medium
                grid ScanCartesian = us.scan
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["double","parallel.Pool","parallel.Cluster"])} = gcp('nocreate');
                kwargs.simdir (1,1) string = tempname(); % simulation directory
                kwargs.f0 (1,1) {mustBeNumeric} = us.tx.fc; % center frequency of the transmit / simulation
                kwargs.CFL_max (1,1) {mustBeReal, mustBePositive} = 0.5 % maximum CFL
                kwargs.txdel (1,1) string {mustBeMember(kwargs.txdel, ["discrete", "continuous", "interpolate"])} = 'interpolate';
            end            

            % create the configuration
            conf_args = rmfield(kwargs, setdiff(fieldnames(kwargs), {'f0', 'CFL_max', 'txdel'}));
            conf_args_ = struct2nvpair(conf_args);
            conf = fullwaveConf(us, medium, grid, conf_args_{:});

            % dispatch
            if isa(kwargs.parenv, "parallel.Cluster")
                % create a job to process it
                job = UltrasoundSystem.fullwaveJob(conf, kwargs.parenv, 'simdir', kwargs.simdir);
                submit(job); % submit the job
                wait(job); % wait for it to finish
            else % run local
                simdir = kwargs.simdir;
                if ~exist(simdir, "dir"), mkdir(simdir); end

                penv = kwargs.parenv;
                if isa(penv, "parallel.ThreadPool"), penv = 0; end % no thread pools for system call
                prj_rt = fullfile(fileparts(mfilename('fullpath')), '..'); % QUPS root
               
                % Write Simulation Files
                write_fullwave_sim(simdir, conf.sim{:}); % all the .dat files
                copyfile(fullfile(prj_rt, 'bin', 'fullwave2_executable'), simdir); % the executable must be here

                % run sim per tx
                wvtx = conf.tx; % splice
                parfor(n = 1:size(wvtx,3), penv)
                    runFullwaveTx(wvtx(:,:,n), simdir, fullfile(simdir, string(n))); 
                end
            end

            % read in the data
            chd = UltrasoundSystem.readFullwaveSim(kwargs.simdir, conf);
        end
    end
    methods(Static,Hidden)
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
            % Note: This functionality is not publicly supported.
            % 
            % Example:
            % 
            % % Setup a system
            % grid = ScanCartesian(...
            % 'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            % 'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', grid, 'xdc', xdc, 'seq', seq);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(grid.size), 1000*ones(grid.size));
            % [Xg, ~, Zg] = grid.getImagingGrid();
            % rho(Xg == 0 & Zg == 30e-3) = 1000*2; % double the density
            % med = Medium.Sampled(grid, c, rho);
            % 
            % % Simulate the ChannelData
            % time_stamp = string(datetime('now'), 'yyyy-MM-dd_HH-mm-ss');
            % simdir = fullfile(pwd, "fwsim", time_stamp); % simulation directory
            % conf = fullwaveConf(us, med, grid, 'CFL_max', 0.5); % configure the sim
            % job = UltrasoundSystem.fullwaveJob(conf, 'simdir', simdir), % create a job
            % submit(job); % submit the job
            % ... do other things
            % wait(job); % wait for the job to finish processing
            % chd = UltrasoundSystem.readFullwaveSim(simdir); % extract the ChannelData
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(hilbert(chd));
            % dbr echo 60;
            % 
            % See also READFULLWAVESIM FULLWAVECONF PARCLUSTER PARALLEL.CLUSTER

            arguments
                conf (1,1) struct
                clu (1,1) parallel.Cluster = parcluster()
                kwargs.simdir (1,1) string = tempname()
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
            % Note: This functionality is not publicly supported.
            % 
            % Example:
            % 
            % % Setup a system
            % grid = ScanCartesian('x', 1e-3*[-35, 35], 'z', 1e-3*[-12, 27]);
            % xdc = TransducerConvex.C5_2v();
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('scan', grid, 'xdc', xdc, 'seq', seq);
            % [grid.dx, grid.dz] = deal(us.lambda / 4);
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(grid.size), 1000*ones(grid.size));
            % pg = grid.positions();
            % rho(argmin(vecnorm(pg - [0 0 20e-3]',2,1))) = 1000*2; % double the density
            % med = Medium.Sampled(grid, c, rho);
            % 
            % % Simulate the ChannelData
            % time_stamp = string(datetime('now'), 'yyyy-MM-dd_HH-mm-ss');
            % simdir = fullfile(pwd, "fwsim", time_stamp); % simulation directory
            % conf = fullwaveConf(us, med, grid, 'CFL_max', 0.5); % configure the sim
            % job = UltrasoundSystem.fullwaveJob(conf, 'simdir', simdir), % create a job
            % submit(job); % submit the job
            % ... do other things
            % wait(job); % wait for the job to finish processing
            % chd = UltrasoundSystem.readFullwaveSim(simdir), % extract the ChannelData
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(hilbert(chd));
            % dbr echo 60;
            % 
            % See also RUNFULLWAVETX ULTRASOUNDSYSTEM/FULLWAVEJOB
            arguments
                simdir (1,1) string
                conf struct {mustBeScalarOrEmpty} = struct.empty
            end

            % see if we can find a conf file in the directory if not given to us
            if isempty(conf)
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
            for m = 1:M % each transmit
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
        end
    end

    % SIMUS calls
    methods
        function chd = simus(us, scat, kwargs, simus_kwargs)
            % SIMUS - Simulate channel data via MUST
            %
            % chd = SIMUS(us, scat) simulates the Scatterers scat and
            % returns a ChannelData object chd.
            %
            % When calling this function the transmit sequence pulse is 
            % ignored. SIMUS only supports tone bursts at the central 
            % frequency of the transducer.
            %
            % The transmit and receive transducer must be identical i.e. 
            % us.rx == us.tx must be true.
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
            % chd = SIMUS(...,Name, Value, ...) passes selected additional
            % Name/Value pairs to simus.m
            % 
            % References:
            % [1] Damien Garcia
            % SIMUS: An open-source simulator for medical ultrasound imaging. 
            % Part I: Theory & examples,
            % Computer Methods and Programs in Biomedicine,
            % Volume 218, 2022, 106726, ISSN 0169-2607,
            % doi: <a href ="matlab:web('https://doi.org/10.1016/j.cmpb.2022.106726')">10.1016/j.cmpb.2022.106726</a>.
            % 
            % [2] Cigier A, Varray F, Garcia D. 
            % SIMUS: An open-source simulator for medical ultrasound imaging. 
            % Part II: Comparison with four simulators. 
            % Comput Methods Programs Biomed. 2022 Jun;220:106774. 
            % doi: <a href="matlab:web('10.1016/j.cmpb.2022.106774')">10.1016/j.cmpb.2022.106774</a>.
            % Epub 2022 Mar 25. PMID: 35398580.
            % 
            % [3] D. Garcia, 
            % "Make the most of MUST, an open-source Matlab UltraSound Toolbox," 
            % 2021 IEEE International Ultrasonics Symposium (IUS), 
            % Xi'an, China, 2021, pp. 1-4, 
            % doi: <a href="matlab:web('10.1109/IUS52206.2021.9593605')">10.1109/IUS52206.2021.9593605</a>.
            % 
            % Example:
            % 
            % % Simulate some data
            % us = UltrasoundSystem(); % get a default system
            % us.fs = single(8 * us.xdc.fc); % must be multiple of 4 for MUST
            % scat = Scatterers('pos', 1e-3*[0 0 30]', 'c0', us.seq.c0); % define a point target
            % chd = simus(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(chd);
            % colorbar;
            % 
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_ALL FOCUSTX GREENS
            
            arguments
                us (1,1) UltrasoundSystem
                scat Scatterers
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate')
                kwargs.periods (1,1) {mustBePositive} = 1
                kwargs.dims {mustBeScalarOrEmpty, mustBeMember(kwargs.dims, [2,3])} = []
                simus_kwargs.FullFrequencyDirectivity (1,1) logical = false % use central freq as reference
                simus_kwargs.ElementSplitting (1,2) {mustBeInteger, mustBePositive} = 1 % element subdivisions
                simus_kwargs.dBThresh (1,1) {mustBeReal} % = -100 % threshold for computing each frequency
                simus_kwargs.FrequencyStep (1,1) {mustBeReal} % = df frequency domain resolution                
                simus_kwargs.WaitBar (1,1) logical = false % add wait bar
            end

            % set options
            simus_kwargs.ParPool = false; % disable parpool within pfield.m

            % load options
            if isempty(kwargs.parenv), kwargs.parenv = 0; end % select 0 workers if empty

            % data typing
            proto = us.fs; % data prototype

            % TODO: check the transmit/receive/sequence impulse: they 
            % cannot be satisfied if not a Delta or empty
            %if ~ismember("periods", varargin(cellfun(@ischar,varargin) | cellfun(@isstring, varargin))) % was this an input?
                %warning("QUPS:UltrasoundSystem:simus:unsatisfiable", "Transmit sequence determined by 'periods', property.");
            %end
            
            % get the points and the dimensions of the simulation(s)
            [X, Y, Z, A] = arrayfun(@(scat) ...
                deal(sub(scat.pos,1,1), sub(scat.pos,2,1), sub(scat.pos,3,1), scat.amp), ...
                scat, 'UniformOutput',false);
            [X,Y,Z,A] = dealfun(@(x) cellfun(@(x){cast(x, 'like', proto)},x), X,Y,Z,A);
            if isempty(kwargs.dims) 
                if all(cellfun(@(Y)all(Y == 0,'all'),Y), 'all') ...
                        && ~isa(us.xdc, 'TransducerMatrix')
                    kwargs.dims = 2; 
                    [Y{:}] = deal([]); % don't simulate in Y if in 1D and it is all zeros 
                else
                    kwargs.dims = 3; 
                end
            end
            if kwargs.dims == 2 && any(cellfun(@(Y)any(Y ~= 0, 'all'),Y))
                warning("QUPS:UltrasoundSystem:simus:casting", "Projecting all points onto Y == 0 for a 2D simulation.");
            end

            % get all other param struct values 
            % (implicitly force same transducer by calling .xdc)
            p = {getSIMUSParam(us.xdc)};
            p = cellfun(@struct2nvpair, p, 'UniformOutput', false);
            p = cat(2, p{:});
            p = struct(p{:});

            % set transmit sequence ... the only way we can
            % TODO: forward arguments to transmit parameters
            p.fs    = us.fs;
            p.TXnow = kwargs.periods; % number of wavelengths
            p.TXapodization = zeros([us.xdc.numel,1], 'like', proto); % set tx apodization
            p.RXdelay = zeros([us.xdc.numel,1], 'like', proto); % receive delays (none)
            
            % set options per Scatterers
            p = repmat(p, [1,1,1,numel(scat)]); % (1 x 1 x 1 x F)
            pxdc = arrayfun(@(scat) {getSIMUSParam(scat)}, scat); % properties per Scatterers
            for f = 1:numel(scat) 
                for fn = string(fieldnames(pxdc{f}))' 
                    p(f).(fn) = pxdc{f}.(fn); % join the structs manually
                end 
            end

            % get transducer offset to offset positions
            off = -us.xdc.offset;
            if isa(us.xdc, 'TransducerConvex')
                off(3) = off(3) + range(sub(us.xdc.positions,3,1)); % offset to the chord
            end

            % choose simus/simus3 handle
            if isa(us.xdc, 'TransducerMatrix'),   mustfun = @simus3;
            else,                                   mustfun = @simus;
            end
 
            % force ElementSplitting scalar for 1D arrays
            if ~isa(us.xdc, 'TransducerMatrix')
                if simus_kwargs.ElementSplitting(2) ~= simus_kwargs.ElementSplitting(1)
                warning( ...
                    "The 2nd value of 'ElementSplitting' is ignored for " ...
                    + class(us.xdc) + " types." ...
                    );
                end
                simus_kwargs.ElementSplitting = simus_kwargs.ElementSplitting(1);
            end

            % select the computing environment
            parenv = kwargs.parenv;
            if isempty(parenv), parenv = 0; end % empty pool -> 0
            isloc = ~isa(parenv, 'parallel.Pool') || ~isa(parenv, 'parallel.Cluster'); % local or parpool
            if isloc, [pclu, parenv] = deal(parenv, 0); else, pclu = 0; end % cluster or local

            % call the sim: FSA approach
            [M, F] = deal(us.xdc.numel, numel(scat)); % splice
            for f = F:-1:1 % per scat
                argf = {X{f}+off(1),Y{f}+off(2),Z{f}+off(3),A{f},zeros([M,1]),p(f),simus_kwargs}; % args per scat
                parfor (m = 1:M, pclu) % use parallel rules, but execute on main thread
                    args = argf; % copy settings for this frame
                    args{6}.TXapodization(m) = 1; % transmit only on element m
                    if isloc, rf{m,f} = mustfun(args{:}); %#ok<PFBNS> % local compute
                    else, out(m,f) = parfeval(parenv, mustfun, 1, args{:}); % add cluster job
                    end
                end
            end
            
            % gather outputs when finished
            if ~isloc, rf = fetchOutputs(out, "UniformOutput",false); end

            % create the output QUPS ChannelData object
            chd = cellfun(@(x) ChannelData('data', x), rf); % per transmit/frame (M x F) object array
            chd = arrayfun(@(f) join(chd(:,f), 3), 1:F); % join over transmits (1 x F) object array
            chd = join(chd, 4); % join over frames (1 x 1) object array
            chd.fs = us.fs; % set the smapling frequency
            chd.t0 = 0; % already offset within call to simus

            % synthesize transmit pulses
            chd = us.focusTx(chd, us.seq, 'interp', kwargs.interp);
        end
    end

    % Field II calls
    methods
        function chd = calc_scat_all(us, scat, element_subdivisions, kwargs)
            % CALC_SCAT_ALL - Simulate channel data via FieldII
            %
            % chd = CALC_SCAT_ALL(us, scat) simulates the Scatterers
            % scat and returns a ChannelData object chd.
            % 
            % chd = CALC_SCAT_ALL(us, scat, element_subdivisions)
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
            % us = UltrasoundSystem('fs', 50e6); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % chd = calc_scat_all(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            % 
            % See also ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/CALC_SCAT_MULTI FOCUSTX

            arguments
                us (1,1) UltrasoundSystem
                scat Scatterers
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])}
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
            end
            
            % set default parenv
            if ~isfield(kwargs, 'parenv')
                hcp = gcp('nocreate');
                if isa(hcp, 'parallel.ThreadPool')
                    kwargs.parenv = 0; % don't default to a thread pool
                else
                    kwargs.parenv = hcp;
                end
            end
            
            % helper function
            vec = @(x) x(:); % column-vector helper function

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = copy(us.tx.impulse); % transmitter impulse
            wv_rx = copy(us.rx.impulse); % receiver impulse
            wv_pl = copy(us.seq.pulse);
            
            % get the time axis (which passes through t == 0)
            [wv_tx.fs, wv_rx.fs, wv_pl.fs] = deal(us.fs);
            t_tx = wv_tx.time;
            t_rx = wv_rx.time;
            t_pl = wv_pl.time;

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.samples)')));
            rx_imp = gather(double(real(vec(wv_rx.samples)')));
            tx_pls = gather(double(real(vec(wv_pl.samples)')));

            % choose the cluster to operate on: avoid running on ThreadPools
            parenv = kwargs.parenv;
            if isempty(parenv), parenv = 0; end
            if isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool') 
                warning('QUPS:InvalidParenv','calc_scat_all cannot be run on a thread-based pool'); % Mex-files ...
                parenv = 0; 
            end

            % splice
            [~, F] = deal(us.seq.numPulse, numel(scat)); % number of transmits/frames
            [fs_, tx_, rx_] = deal(gather(us.fs), us.tx, us.rx); % splice
            [c0, pos, amp] = arrayfun(@(t)deal(t.c0, {t.pos}, {t.amp}), scat); % splice

            % Make position/amplitude and transducers constants across the workers
            if isa(parenv, 'parallel.Pool') 
                cfun = @parallel.pool.Constant;
            else 
                cfun = @(x)struct('Value', x);
                [pos, amp] = deal({pos},{amp}); % for struct to work on cell arrays
            end
            [pos_, amp_, tx_, rx_] = dealfun(cfun, pos, amp, tx_, rx_);

            % for each scat (frame)
            parfor (f = 1:F, parenv) % each transmit pulse
            % reinitialize field II
            field_init(-1);
            
            % set sound speed/sampling rate
            set_field('fs', double(fs_));
            set_field('c',c0(f));

            % get Tx/Rx apertures
            p_focal = [0;0;0];
            Tx = tx_.Value.getFieldIIAperture(element_subdivisions, p_focal.'); %#ok<PFBNS> % constant over workers
            Rx = rx_.Value.getFieldIIAperture(element_subdivisions, p_focal.'); %#ok<PFBNS> % constant over workers

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
            chd = us.focusTx(chd, us.seq, 'interp', kwargs.interp);
        end        

        function [chd, rfun] = calc_scat_multi(us, scat, element_subdivisions, kwargs)
            % CALC_SCAT_MULTI - Simulate channel data via FieldII
            %
            % chd = CALC_SCAT_MULTI(us, scat) simulates the Scatterers 
            % scat and returns a ChannelData object chd.
            %
            % chd = CALC_SCAT_MULTI(us, scat, element_subdivisions)
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
            % job = CALC_SCAT_MULTI(..., 'parenv', clu, 'job', true)
            % returns a job on the parcluster clu for each scatterer. The
            % job can then be run with `submit(job)`.
            % 
            % [job, rfun] = CALC_SCAT_MULTI(...) also returns a function
            % rfun to read the output of the job into a ChannelData object
            % once the job has been completed.
            % 
            % job = CALC_SCAT_MULTI(..., 'job', true) uses the default
            % cluster returned by parcluster().
            % 
            % job = CALC_SCAT_MULTI(..., 'job', true, 'type', 'Independent')
            % creates an independent job via `createJob`.
            % 
            % job = CALC_SCAT_MULTI(..., 'job', true, 'type', 'Communicating')
            % creates a communicating job via `createCommunicatingJob`.
            % This is the default.
            % 
            % An independent job transfers and stores a separate set of
            % data for each task, allowing each transmit to be simulated on
            % a different worker which may exist on different nodes, but it
            % requires a full copy of the inputs for each worker which may
            % incur a heavy data transfer and storage cost.
            % 
            % A communicating job shares all resources defined by the
            % parcluster across all workers, which reduces data transfer
            % and storage costs, but may require the entire job to fit on a
            % single node.
            % 
            % The parcluster should be configured according to the job. For
            % a communicating job, the NumWorkers and any memory options
            % should be configured for all transmits. For an independent
            % job, the NumWorkers and any memory options should be
            % allocated for an individual transmit.
            % 
            % References:
            % [1] Jensen, Jørgen Arendt.
            % <a href="matlab.web('https://orbit.dtu.dk/en/publications/field-a-program-for-simulating-ultrasound-systems')">"Field: A program for simulating ultrasound systems."</a>
            % Medical & Biological Engineering & Computing 
            % 34.sup. 1 (1997): 351-353.
            % 
            % Example:
            % % Simulate some data
            % us = UltrasoundSystem('fs', 50e6); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % chd = calc_scat_multi(us, scat); % simulate the ChannelData
            % 
            % % Display the data
            % figure;
            % imagesc(real(chd));
            % colorbar;
            % 
            % % Run on a parcluster instead
            % clu = parcluster();
            % [job, rfun] = calc_scat_multi(us, scat, 'job', true);
            % 
            % ... modify the job ...
            % 
            % % simulate data
            % submit(job); % launch
            % wait(job); % wait for the job to finish
            % chd2 = rfun(job); % read data
            % 
            % isalmostn(chd.data, chd2.data)
            % 
            % See also ULTRASOUNDSYSTEM/SIMUS ULTRASOUNDSYSTEM/FOCUSTX

            arguments
                us (1,1) UltrasoundSystem
                scat Scatterers
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
                kwargs.fc double {mustBeNonnegative} = mean([us.rx.fc, us.tx.fc]); % attenuation central frequency
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])}
                kwargs.job (1,1) logical = false
                kwargs.type (1,1) string {mustBeMember(kwargs.type, ["Communicating", "Independent"])} = "Communicating";
                kwargs.verbose (1,1) logical = false
            end
            % alias
            [kvb, fcs] = deal(kwargs.verbose, kwargs.fc);

            % validate inputs
            if ~any(numel(kwargs.fc) == [1 numel(scat)])
                error("QUPS:calc_scat_multi:vectorCentralFrequency", ...
                    "The number of central frequencies ("+numel(kwargs.fc)...
                    +") must be scalar or equivalent to the number of Scatterers (" ...
                    +numel(scat)+")");
            elseif isscalar(kwargs.fc)
                kwargs.fc = repmat(kwargs.fc, size(scat)); % vectorize
            end

            % set default parenv
            if ~isfield(kwargs, 'parenv')
                hcp = gcp('nocreate');
                if kwargs.job
                    kwargs.parenv = parcluster(); % default cluster
                elseif isa(hcp, 'parallel.ThreadPool')
                    kwargs.parenv = 0; % don't default to a thread pool
                else
                    kwargs.parenv = hcp; % use the current pool
                end
            end

            % make a (communicating) job(s) instead
            if kwargs.job
                switch kwargs.type
                    case "Communicating"
                        job = createCommunicatingJob(kwargs.parenv, 'AutoAddClientPath', true, 'AutoAttachFiles',true, 'Type', 'Pool');
                        job.createTask(@(us,scat) ...
                            us.calc_scat_multi(scat, element_subdivisions ...
                            , 'parenv', Inf, 'fc', fcs, 'verbose', kvb ...
                            ), 1, {us, scat}, 'CaptureDiary',true);

                        % anonymous function to read in data
                        rfun = @(job) job.Tasks(1).OutputArguments{1};

                    case "Independent"
                        for s = 1:numel(scat) % for each Scatterers
                            seqs = splice(us.seq,1); % split up into individual sequences
                            job(s) = createJob(kwargs.parenv, 'AutoAddClientPath', true, 'AutoAttachFiles',true); %#ok<AGROW>
                            us_ = copy(us); % copy semantics
                            for m = 1:us.seq.numPulse % each tx
                                us_.seq = seqs(m); % switch to next tx
                                job(s).createTask(@(us,scat,fc) ...
                                    us.calc_scat_multi(scat, element_subdivisions ...
                                    , 'parenv', 0, 'fc', fc, 'verbose', kvb ...
                                    ), 1, {us, scat(s), fcs(s)}, 'CaptureDiary',true);
                            end

                            % anonymous function to read in data
                            rfun = @(jobs) ...
                                join( arrayfun(@(job) ... for each job
                                join( arrayfun(@(tsk) ... for each tsk
                                tsk.OutputArguments{1}, ... extract data
                                job.Tasks), 3), ... and join tasks (txs) in dim 3
                                jobs), 4); % and join jobs (scats) in dim 4
                        end

                end
                chd = job; % alias
                return;
            end

            % helper function
            vec = @(x) x(:); % column-vector helper function

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = copy(us.tx.impulse); % transmitter impulse
            wv_rx = copy(us.rx.impulse); % receiver impulse
            wv_pl = copy(us.seq.pulse);

            % get the time axis (which passes through t == 0)
            [wv_tx.fs, wv_rx.fs, wv_pl.fs] = deal(gather(us.fs));
            t_tx = wv_tx.time;
            t_rx = wv_rx.time;
            t_pl = wv_pl.time;

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.samples)')));
            rx_imp = gather(double(real(vec(wv_rx.samples)')));
            tx_pls = gather(double(real(vec(wv_pl.samples)')));

            % get the apodization and time delays across the aperture
            apd_tx = gather( us.seq.apodization(us.tx)); % N x M
            tau_tx = gather(-us.seq.delays(     us.tx)); % N x M
            tau_offset = min(tau_tx, [], 1); % (1 x M)
            tau_tx = tau_tx - tau_offset; % 0-base the delays for FieldII

            % choose the parallel environment to operate on: avoid running on ThreadPools
            parenv = kwargs.parenv;
            if isempty(parenv), parenv = 0; end
            if isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool')
                warning("QUPS:calc_scat:invalidParalelEnvironment", ...
                    "calc_scat_multi cannot be run on a "+class(parenv)+"."); % Mex-files ...
                parenv = 0; 
            end

            % splice
            [M, F] = deal(size(tau_tx,2), numel(scat)); % number of transmits/frames
            [fs_, tx_, rx_] = deal(gather(us.fs), us.tx, us.rx); % splice
            [c0, att, pos, amp] = arrayfun(@(t)deal(t.c0, t.alpha0, {t.pos}, {t.amp}), scat); % splice
            if ~isscalar(kwargs.fc), fc_ = kwargs.fc; else, fc_ = repmat(kwargs.fc, size(scat)); end

            % Make position/amplitude and transducers constants across the workers
            if isa(parenv, 'parallel.Pool')
                cfun = @parallel.pool.Constant;
                if kwargs.verbose, tb = ticBytes(parenv); end
            else
                cfun = @(x)struct('Value', x);
                [pos, amp] = deal({pos},{amp}); % for struct to work on cell arrays
            end
                 
            kvb = kwargs.verbose; % splice
            [pos_, amp_, tx_, rx_] = dealfun(cfun, pos, amp, tx_, rx_);
            [voltages, ts] = deal(cell(M,F)); % pre-allocate
            parfor (m = 1:M, parenv) % each transmit pulse
                % for (m = 1:M) % each transmit pulse %%% DEBUG %%%
                % (re)initialize field II
                field_init(-1);

                % get Tx/Rx apertures
                p_focal = [0;0;0];
                Tx = tx_.Value.getFieldIIAperture(element_subdivisions, p_focal.'); %#ok<PFBNS> % constant over workers
                Rx = rx_.Value.getFieldIIAperture(element_subdivisions, p_focal.'); %#ok<PFBNS> % constant over workers

                % set the impulse response function / excitation function
                xdc_impulse   (Tx, tx_imp);
                xdc_impulse   (Rx, rx_imp);
                xdc_excitation(Tx, tx_pls);

                % nullify the response on the receive aperture
                xdc_times_focus(Rx, 0,  zeros([1, rx_.Value.numel])); % set the delays manually

                % for each transmit, set the transmit time delays and
                % apodization
                xdc_times_focus(Tx, 0, double(tau_tx(:,m)')); % set the delays
                xdc_apodization(Tx, 0, double(apd_tx(:,m)')); % set the apodization

                for (f = 1:F) %#ok<NO4LP> % each scat frame
                    % set sound speed/sampling rate
                    set_field('fs', double(fs_));
                    set_field('c', c0(f)); %#ok<PFBNS> % this array is small

                    % set the attenuation
                    use_att = isfinite(att(f)) && isfinite(fc_(f)); %#ok<PFBNS> this array is small too
                    if use_att % TODO: no need to branch?
                        set_field('att',fc_(f)*att(f)); % attenuation, frequency independent
                        set_field('freq_att',  att(f)); % attenuation, frequency   dependent
                        set_field('att_f0'  ,  fc_(f)); % attenuation, central frequency
                        set_field('use_att' , double(use_att)); % turn on attenuation modelling                        
                    end
                    if kvb, disp("Simulating pulse "+m+" of frame "+f+" "+sub(["without","with"],1+use_att) + " attenuation."); end


                    % call the sim
                    [voltages{m,f}, ts{m,f}] = calc_scat_multi(Tx, Rx, pos_.Value{f}.', amp_.Value{f}.'); %#ok<PFBNS> % constant over workers
                end
            end

            if isa(parenv, 'parallel.Pool') && kwargs.verbose, tocBytes(parenv, tb); end

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
            chd.fs = us.fs;
            chd.t0 = chd.t0 + t0;
        end

        function [p, t] = calc_hp(us, scat, element_subdivisions, kwargs)
            % CALC_HP - Simulate the pressure field/sensitivity via FieldII
            %
            % p = CALC_HP(us) simulates the transmit pressure field p
            % emitted by the UltrasoundSystem us. p is returned as a 
            % (T x I x M) ND-array where T is time, I is the pixels of
            % us.scan, and M is each transmit.
            % 
            % p = CALC_HP(us, scat) simulates the transmit pressure field p
            % emitted by the UltrasoundSystem us into the homogeneous
            % Medium or Scatterers scat. Only the ambient sound speed
            % and attenuation affect the pressure field calculation. The
            % default is Medium('c0', us.seq.c0).
            % 
            % p = CALC_HP(us, med, element_subdivisions) specifies the
            % number of subdivisions in width and height for each element.
            % The default is [1, 1].
            %
            % p = CALC_HP(..., "ap", aperture) specifies whether to
            % simulate the transmit aperture, receive aperture, or the
            % both (2-way) apertures. The latter calls `calc_hhp()`. The
            % default is "tx".
            %
            % p = CALC_HP(..., "rxseq", rxseq) specifies a receive Sequence
            % rxseq for the receive beam focusing (or analagously delays
            % and apodization). The default is us.seq (identical to the
            % transmit sequence).
            %
            % [p, t] = CALC_HP(...) additionally returns the time axes t.
            % 
            % [...] = CALC_HP(...,'c0', c0) sets the sound speed c0.
            % The default is us.seq.c0.
            %
            % [...] = CALC_HP(..., 'parenv', clu) or 
            % [...] = CALC_HP(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to
            % parallelize computations. parallel.ThreadPools are invalid
            % due to mex function restrictions.
            % 
            % [...] = CALC_HP(..., 'parenv', 0) avoids using a 
            % parallel.Cluster or parallel.Pool. Use 0 when operating on a 
            % GPU or if memory usage explodes on a parallel.ProcessPool.
            %
            % The default is the current pool returned by gcp.
            % 
            % job = CALC_HP(..., 'parenv', clu, 'job', true)
            % returns a job on the parcluster clu for each scatterer. The
            % job can then be run with `submit(job)`.
            % 
            % job = CALC_HP(..., 'job', true) uses the default
            % cluster returned by parcluster().
            % 
            % job = CALC_HP(..., 'job', true, 'type', 'Independent')
            % creates an independent job via `createJob`.
            % 
            % job = CALC_HP(..., 'job', true, 'type', 'Communicating')
            % creates a communicating job via `createCommunicatingJob`.
            % This is the default.
            % 
            % An independent job transfers and stores a separate set of
            % data for each task, allowing each transmit to be simulated on
            % a different worker which may exist on different nodes, but it
            % requires a full copy of the inputs for each worker which may
            % incur a heavy data transfer and storage cost.
            % 
            % A communicating job shares all resources defined by the
            % parcluster across all workers, which reduces data transfer
            % and storage costs, but may require the entire job to fit on a
            % single node.
            % 
            % The parcluster should be configured according to the job. For
            % a communicating job, the NumWorkers and any memory options
            % should be configured for all transmits. For an independent
            % job, the NumWorkers and any memory options should be
            % allocated for an individual transmit.
            % 
            % References:
            % [1] Jensen, Jørgen Arendt.
            % <a href="matlab.web('https://orbit.dtu.dk/en/publications/field-a-program-for-simulating-ultrasound-systems')">"Field: A program for simulating ultrasound systems."</a>
            % Medical & Biological Engineering & Computing 
            % 34.sup. 1 (1997): 351-353.
            % 
            % Example:
            % % Setup a Focused Transmit
            % xdc = TransducerArray.L12_3v(); % Verasonics L12-3v
            % seq = Sequence('type', 'FC', 'c0', 1540); % Sequence @ 1540 m/s
            % us = UltrasoundSystem('xdc', xdc, 'seq', seq); % get a system
            % 
            % % Configure
            % [us.scan.dx, us.scan.dz] = deal(us.lambda / 2); % set field resolution
            % us.seq.focus = [0 0 xdc.el_focus]'; % set the transmit focus
            % 
            % % Simulate
            % [p, t] = calc_hp(us); % simulate the transmit pressure field
            % p = hilbert(p); % analytic
            %
            % % Display the data
            % figure;
            % hi = imagesc(us.scan, reshape(p(end,:), us.scan.size));
            % cbd = mod2db(prctile(abs(p(:)), [0.01 99.99])); % color axes
            % dbr echo; clim(gather(cbd)); % set colors
            %
            % % Animate over time
            % pi = reshape(p, [size(p,1), us.scan.size]); % as an image
            % ttls = "Transmit Field Energy" + newline ...
            %      + "Time: t = "+(1e6*t)+" us"; % titles
            % animate(pi, hi, "title", ttls, "loop", false, "fs", Inf);
            % 
            % % Example 2:
            % % Run the sim on a parcluster instead
            % clu = parcluster();
            % job = calc_hp(us, 'job', true, 'parenv', clu);
            % 
            % ... modify the job ...
            % 
            % % simulate data
            % submit(job); % launch
            % wait(job); % wait for the job to finish
            % out = fetchOutputs(job); % read data
            % [p2, t2] = deal(out{:}); % parse
            % 
            % isalmostn(p, p2) % verify
            % 
            % See also ULTRASOUNDSYSTEM/CALC_SCAT_ALL ULTRASOUNDSYSTEM/CALC_SCAT_MULTI

            arguments
                us (1,1) UltrasoundSystem
                scat {mustBeA(scat, ["Medium","Scatterers"])} = Medium('c0', us.seq.c0);
                element_subdivisions (1,2) double {mustBeInteger, mustBePositive} = [1,1]
                kwargs.ap (1,1) string {mustBeMember(kwargs.ap, ["tx","rx","2way"])} = "tx"
                kwargs.rxseq (1,1) Sequence = us.seq % foci for the receive beams
                kwargs.c0 (1,1) double = us.seq.c0 % sound speed
                kwargs.fc double {mustBeNonnegative} % attenuation central frequency
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])}
                kwargs.job (1,1) logical = false
                kwargs.type (1,1) string {mustBeMember(kwargs.type, ["Communicating", "Independent"])} = "Communicating";
                kwargs.verbose (1,1) logical = false
            end
            
            % validate inputs
            if ~isfield(kwargs, "fc") % initialize if not already
                switch kwargs.ap
                    case "tx", kwargs.fc =       us.tx.fc            ;
                    case "rx", kwargs.fc =                 us.rx.fc  ;
                    otherwise, kwargs.fc = mean([us.tx.fc, us.rx.fc]);
                end
            end

            % one freq per scat
            if ~any(numel(kwargs.fc) == [1 numel(scat)])
                error("QUPS:calc_hp:vectorCentralFrequency", ...
                    "The number of central frequencies ("+numel(kwargs.fc)...
                    +") must be scalar or equivalent to the number of Scatterers (" ...
                    +numel(scat)+")");
            elseif isscalar(kwargs.fc)
                kwargs.fc = repmat(kwargs.fc, size(scat)); % vectorize
            end

            % must be the same number of rxseq as txseq
            if kwargs.ap ~= "tx" && (us.seq ~= kwargs.rxseq || us.seq.numPulse ~= kwargs.rxseq.numPulse)
                error("QUPS:calc_hp:invalidRxSequence", ...
                    "The receive focal sequence must have an identical number of foci ("+us.seq.numPulse ...
                    +") as the transmit focal sequence ("+kwargs.rxseq.numPulse+")" ...
                    );
            end

            % set default parenv
            if ~isfield(kwargs, 'parenv')
                hcp = gcp('nocreate');
                if kwargs.job
                    kwargs.parenv = parcluster(); % default cluster
                elseif isa(hcp, 'parallel.ThreadPool')
                    kwargs.parenv = 0; % don't default to a thread pool
                else
                    kwargs.parenv = hcp; % use the current pool
                end
            end

            % alias
            [kvb, fcs, rxseq] = deal(kwargs.verbose, kwargs.fc, kwargs.rxseq);

            % make a (communicating) job(s) instead
            if kwargs.job
                fkeep = ["c0", "ap", "verbose"]; % fields to keep
                opts = namedargs2cell(rmfield(kwargs, setdiff(string(fieldnames(kwargs)), fkeep))); % all other fields
                switch kwargs.type
                    case "Communicating"
                        job = createCommunicatingJob(kwargs.parenv, 'AutoAddClientPath', true, 'AutoAttachFiles',true, 'Type', 'Pool');
                        job.createTask(@(us, scat) ...
                            us.calc_hp(scat, element_subdivisions ...
                            , 'parenv', Inf, 'fc', fcs, 'rxseq', rxseq, opts{:} ...
                            ), 2, {us, scat}, 'CaptureDiary',true);

                        % anonymous function to read in data
                        rfun = @fetchOutputs;
                        % rfun = @(job) job.Tasks(1).OutputArguments{1};

                    case "Independent"
                        for s = 1:numel(scat) % for each Scatterers
                            seqs = splice(us.seq,1); % split up into individual sequences
                            rsqs = splice(rxseq ,1); % split up into individual sequences
                            job(s) = createJob(kwargs.parenv, 'AutoAddClientPath', true, 'AutoAttachFiles',true); %#ok<AGROW>
                            us_ = copy(us); % copy semantics
                            for m = 1:us.seq.numPulse % each tx
                                us_.seq = seqs(m); % switch to next tx
                                job(s).createTask(@(us,scat,fc,rsq) ...
                                    us.calc_hp(scat, element_subdivisions ...
                                    , 'parenv', 0, 'fc', fc, 'rxseq', rsq, opts{:} ...
                                    ), 2, {us, scat(s), fcs(s), rsqs(s)}, 'CaptureDiary',true);
                            end

                            % anonymous function to read in data
                            rfun = @fetchOutputs;
                            % rfun = @(jobs) ...
                            %     join( arrayfun(@(job) ... for each job
                            %     join( arrayfun(@(tsk) ... for each tsk
                            %     tsk.OutputArguments{1}, ... extract data
                            %     job.Tasks), 3), ... and join tasks (txs) in dim 3
                            %     jobs), 4); % and join jobs (scats) in dim 4
                        end

                end
                [p, t] = deal(job, rfun); % alias
                return;
            end

            % helper function
            vec = @(x) x(:); % column-vector helper function

            % get the Tx/Rx impulse response function / excitation function
            wv_tx = copy(us.tx.impulse); % transmitter impulse
            wv_rx = copy(us.rx.impulse); % receiver impulse
            wv_pl = copy(us.seq.pulse); % excitation pulse

            % get the time axis (which passes through t == 0)
            [wv_tx.fs, wv_rx.fs, wv_pl.fs] = deal(gather(us.fs));
            t_tx = wv_tx.time;
            t_rx = wv_rx.time;
            t_pl = wv_pl.time;

            % define the impulse and excitation pulse
            tx_imp = gather(double(real(vec(wv_tx.samples)')));
            rx_imp = gather(double(real(vec(wv_rx.samples)')));
            tx_pls = gather(double(real(vec(wv_pl.samples)')));

            % get the apodization and time delays across the aperture(s)
            apd_tx = gather( us.seq.apodization(us.tx)); % N x M
            tau_tx = gather(-us.seq.delays(     us.tx)); % N x M
            tau_offset_tx = min(tau_tx, [], 1); % (1 x M)
            tau_tx = tau_tx - tau_offset_tx; % 0-base the delays for FieldII

            apd_rx = gather( rxseq.apodization(us.rx)); % N x M
            tau_rx = gather(-rxseq.delays(     us.rx)); % N x M
            tau_offset_rx = min(tau_rx, [], 1); % (1 x M)
            tau_rx = tau_rx - tau_offset_rx; % 0-base the delays for FieldII

            % choose the parallel environment to operate on: avoid running on ThreadPools
            parenv = kwargs.parenv;
            if isempty(parenv), parenv = 0; end
            if isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool')
                warning("QUPS:calc_scat:invalidParalelEnvironment", ...
                    "calc_scat_multi cannot be run on a "+class(parenv)+"."); % Mex-files ...
                parenv = 0; 
            end

            % splice
            [M, F] = deal(size(tau_tx,2), numel(scat)); % number of transmits/frames
            [fs_, tx_, rx_, ap] = deal(gather(us.fs), us.tx, us.rx, kwargs.ap); % splice
            [c0, att] = arrayfun(@(t) deal(t.c0, t.alpha0), scat); % splice
            pos = {reshape(cast(us.scan.positions(), 'double'), 3, [])}; % field positions {(3 x I)}
            if ~isscalar(kwargs.fc), fc_ = kwargs.fc; else, fc_ = repmat(kwargs.fc, size(scat)); end

            % Make position/amplitude and transducers constants across the workers
            if isa(parenv, 'parallel.Pool')
                cfun = @parallel.pool.Constant;
                if kwargs.verbose, tb = ticBytes(parenv); end
            else
                cfun = @(x)struct('Value', x);
                [pos] = {pos}; % for struct to work on cell arrays
            end
                 
            [pos_, tx_, rx_] = dealfun(cfun, pos, tx_, rx_);
            [pres, ts] = deal(cell(M,F)); % pre-allocate
            parfor (m = 1:M, parenv) % each transmit pulse
            % for (m = 1:M) % each transmit pulse %%% DEBUG %%%
                % (re)initialize field II
                field_init(-1);

                % get Tx/Rx apertures
                p_focal = [0;0;0];
                Tx = tx_.Value.getFieldIIAperture(element_subdivisions, p_focal.'); %#ok<PFBNS> % constant over workers
                Rx = rx_.Value.getFieldIIAperture(element_subdivisions, p_focal.'); %#ok<PFBNS> % constant over workers

                % set the impulse response function / excitation function
                xdc_impulse   (Tx, tx_imp);
                xdc_impulse   (Rx, rx_imp);
                xdc_excitation(Tx, tx_pls);

                % nullify the response on the receive aperture
                % xdc_times_focus(Rx, 0,  zeros([1, rx_.Value.numel])); % set the delays manually
                % for each receive, set the receive delays and apodization
                xdc_times_focus(Rx, 0, double(tau_rx(:,m)')); % set the delays
                xdc_apodization(Rx, 0, double(apd_rx(:,m)')); % set the apodization

                % for each transmit, set the transmit time delays and
                % apodization
                xdc_times_focus(Tx, 0, double(tau_tx(:,m)')); % set the delays
                xdc_apodization(Tx, 0, double(apd_tx(:,m)')); % set the apodization

                for (f = F:-1:1) %#ok<NO4LP> % each scat frame
                    % set sound speed/sampling rate
                    set_field('fs', double(fs_));
                    set_field('c', c0(f)); %#ok<PFBNS> % this array is small

                    % set the attenuation
                    use_att = isfinite(att(f)) && isfinite(fc_(f)); %#ok<PFBNS> this array is small too
                    if use_att % TODO: no need to branch?
                        set_field('att',fc_(f)*att(f)); % attenuation, frequency independent
                        set_field('freq_att',  att(f)); % attenuation, frequency   dependent
                        set_field('att_f0'  ,  fc_(f)); % attenuation, central frequency
                        set_field('use_att' , double(use_att)); % turn on attenuation modelling                        
                    end
                    if kvb, disp("Simulating pulse "+m+" of frame "+f+" "+sub(["without","with"],1+use_att) + " attenuation."); end


                    % call the sim
                    % [voltages{m,f}, ts{m,f}] = calc_scat_multi(Tx, Rx, pos_.Value{f}.', amp_.Value{f}.'); %#ok<PFBNS> % constant over workers
                    switch ap
                        case "tx"  , [pres{m,f}, ts{m,f}] = calc_hp( Tx   , pos_.Value{f}.');  %#ok<PFBNS> % constant over workers
                        case "rx"  , [pres{m,f}, ts{m,f}] = calc_hp( Rx   , pos_.Value{f}.');              % constant over workers
                        case "2way", [pres{m,f}, ts{m,f}] = calc_hhp(Tx,Rx, pos_.Value{f}.');              % constant over workers
                    end
                end
            end

            if isa(parenv, 'parallel.Pool') && kwargs.verbose, tocBytes(parenv, tb); end

            % adjust start time based on signal time definitions
            switch ap
                case "tx", tau_offset = tau_offset_tx                ;
                case "rx", tau_offset =                 tau_offset_rx;
                otherwise, tau_offset = tau_offset_tx + tau_offset_rx;
            end
            t0 = ... cell2mat(ts) + ... % fieldII start time (1 x 1 x M x F)
                (t_pl(1) + t_tx(1) + t_rx(1)) ... signal delays for impulse/excitation
                + shiftdim(tau_offset,-1) ... 0-basing the delays across the aperture
                ; % 1 x 1 x M
            
            % create the output QUPS ChannelData object 
            chd = cellfun(@(x, t0) ChannelData('data', x, 't0', t0, 'fs', us.fs), pres, ts); % per transmit/frame (M x F) object array
            chd = arrayfun(@(f) join(chd(:,f), 3), 1:F); % join over transmits (1 x F) object array
            if ~isscalar(chd), chd = join(chd, 4); end % join over frames (1 x 1) object array

            % set sampling frequency and transmit times for all
            chd.fs = us.fs;
            chd.t0 = chd.t0 + t0;

            % output (T x I1 x I2 x I3 x M x F)
            [p, t] = deal(chd.data, chd.time);
        end
    end
    
    % k-Wave calls
    methods
        function [chd, readfun, args] = kspaceFirstOrder(us, med, cgrd, varargin, kwargs, karray_args, kwave_args, ksensor_args)
            % KSPACEFIRSTORDER - Simulate channel data via k-Wave
            % 
            % chd = KSPACEFIRSTORDER(us, med) simulates the Medium med 
            % and returns a ChannelData object chd via k-Wave.
            %
            % chd = KSPACEFIRSTORDER(us, med, cgrd) operates using
            % the simulation region defined by the ScanCartesian cgrd. The
            % default is us.scan.
            % 
            % If using an element mapping method provided by kWaveArray,
            % the step sizes in all dimensions must be identical.
            %
            % chd = KSPACEFIRSTORDER(..., 'T', T) runs the simulation
            % until time T. The default is computed heuristically to
            % include a two-way propagation of a signal at the slowest
            % sound speed.
            % 
            % chd = KSPACEFIRSTORDER(..., 'isosub', true) runs an 
            % additional simulation per pulse of an analagous medium with 
            % constant impdeance in order to subtract any background 
            % artefacts from the data. The default is false.
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
            % lead to instability. The default is 0.25.
            % 
            % chd = KSPACEFIRSTORDER(..., 'ElemMapMethod', method) 
            % specifies the computational method used for mapping 
            % transducer elements to the computational grid. Must be one of 
            % {'nearest*, 'linear', 'karray-direct', 'karray-depend'}. 
            % 
            % The 'nearest' method uses the nearest pixel. The 'linear' 
            % method uses linear interpolation weights. The 'karray-direct'
            % method uses kWaveArray methods but avoids recomputing 
            % intermediate results. The 'karray-depend' method always uses 
            % the kWaveArray methods, but can be slower.
            % 
            % chd = KSPACEFIRSTORDER(..., 'parenv', clu) or 
            % chd = KSPACEFIRSTORDER(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to compute
            % each pulse in parallel via parfor. 
            % 
            % When using a parallel.Cluster clu, the settings for clu must
            % be setup for an independent parallel.Job where an indepedent
            % task is created for each pulse. For example, if using a SLURM
            % cluster with multiple GPUs, clu.SubmitArguments should
            % contain ' --gpus=1' so that 1 GPU is requested per pulse.
            % 
            % [job, readfun] = KSPACEFIRSTORDER(..., 'parenv', clu) instead 
            % returns a parallel.CommunicatingJob job and a function to
            % read the ChannelData object from the completed job, readfun. 
            % 
            % The settings for clu should must be setup for a
            % parallel.CommunicatingJob as a single MPI job is created for
            % all pulses. For example, if using a SLURM cluster with
            % multiple GPUs, clu.SubmitArguments should contain the number
            % of GPUs desired for execution e.g. ' --gpus=4' and
            % clu.NumWorkers should equal the maximum number of
            % simulataneous pulses to use e.g. 8. If more CPUs than GPUs
            % are requested, GPUs will be shared between multiple CPUs.
            % Sharing GPUs between multiple CPUs can be faster, but
            % requires sufficient GPU RAM.
            % 
            % Use `submit(job)` to begin running the job. 
            % Use `wait(job)` to wait for the job to complete. 
            % When the job has completed successfully, use 
            % `chd = readfun(job)` to extract the ChannelData object.
            %
            % A MATLAB parallel.Job is saved until it is deleted, so
            % simulations can be recalled later from the job reference.
            % 
            % chd = KSPACEFIRSTORDER(..., 'parenv', 0) avoids using a
            % parallel.Cluster or parallel.Pool. 
            %
            % The default is the current pool returned by `gcp('nocreate')`.
            % 
            % [...] = KSPACEFIRSTORDER(..., 'binary', true) calls k-Wave's 
            % accelerated binary functions, which are appended with 'C' for
            % the CPU version and 'G' for the GPU version. The selection is
            % determined by whether the 'DataCast' option represent a GPU 
            % or CPU data type. The path to the binaries should be
            % specified using the 'BinaryPath' option.
            %
            % Note: If using this option with a cluster, the binaries must 
            % be accesible and executable from the worker session, not the 
            % current MATLAB session.
            %
            % Note: The binary option is (currently) incompatible with a
            % parallel.ThreadPool. If `gcp('nocreate')` returns a
            % parallel.ThreadPool, you must set the 'parenv' option
            % explicitly.
            % 
            % [...] = KSPACEFIRSTORDER(..., Name, Value, ...)
            % specifies other Name/Value pairs that are valid for 
            % kWaveArray's constructor or kWave's kspaceFirstOrderND. 
            % 
            % kWave's kspaceFirstOrderND 'PMLSize' and 'PMLInside'
            % arguments are invalid as they are overwritten by this method.
            %
            % References:
            % [1] B. E. Treeby and B. T. Cox, 
            % "k-Wave: MATLAB toolbox for the simulation and reconstruction of photoacoustic wave-fields,"
            % J. Biomed. Opt., vol. 15, no. 2, p. 021314, 2010.
            % <a href="matlab:web('https://doi.org/10.1117/1.3360308')">DOI: 10.1117/1.3360308</a>
            % 
            % [2] B. E. Treeby, J. Budisky, E. S. Wise, J. Jaros, and B. T. Cox, 
            % "Rapid calculation of acoustic fields from arbitrary continuous-wave sources," 
            % J. Acoust. Soc. Am., vol. 143, no. 1, pp. 529-537, 2018.
            % <a href="matlab:web('https://doi.org/10.1121/1.5021245')">DOI: 10.1121/1.5021245</a>
            % 
            % Example:
            % 
            % % Setup a system
            % xdc = TransducerArray();
            % cgrd = ScanCartesian( ...
            %     'x', 1e-3*(-10 : 0.1 : 10), ...
            %     'z', 1e-3*(  0 : 0.1 : 30) ...
            % );
            % seq = SequenceRadial('angles', 0, 'ranges', 1); % plane-wave
            % us = UltrasoundSystem('xdc', xdc, 'seq', seq, 'fs', single(16*xdc.fc));
            % 
            % % Create a Medium to simulate
            % [c, rho] = deal(1500*ones(cgrd.size), 1000*ones(cgrd.size));
            % [Xg, ~, Zg] = cgrd.getImagingGrid();
            % rho(Xg == 0 & Zg == 10e-3) = 1000*2; % add a density scatterer
            % med = Medium.Sampled(cgrd, c, rho, 'c0', 1500, 'rho0', 1000);
            % 
            % % Simulate the ChannelData
            % chd = kspaceFirstOrder(us, med, cgrd);
            % 
            % % Display the ChannelData
            % figure;
            % imagesc(hilbert(chd));
            % 
            % See also ULTRASOUNDSYSTEM/FULLWAVESIM PARALLEL.JOB/FETCHOUTPUTS
            arguments % required arguments
                us (1,1) UltrasoundSystem
                med (1,1) Medium = Medium("c0", us.seq.c0, "rho0", us.seq.c0 / 1.5);
                cgrd (1,1) ScanCartesian = us.scan
            end
            arguments(Repeating)
                varargin % unidentified arguments passed to k-Wave directly
            end
            arguments % keyword arguments for this function
                kwargs.T double {mustBeScalarOrEmpty} = [], % simulation time (s)
                kwargs.PML (1,:) double = [20 56], % (one-sided) PML size range
                kwargs.CFL_max (1,1) double {mustBePositive} = 0.25, % maximum cfl number (for stability)
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate'), % parallel environment for running simulations
                kwargs.ElemMapMethod (1,1) string {mustBeMember(kwargs.ElemMapMethod, ["nearest","linear","karray-direct", "karray-depend"])} = 'nearest', % one of {'nearest'*,'linear','karray-direct', 'karray-depend'}
                kwargs.el_sub_div (1,2) double = max(us.tx.getLambdaSubDiv(med.c0), us.rx.getLambdaSubDiv(med.c0)), % element subdivisions (width x height)
                kwargs.bsize (1,1) double {mustBeInteger, mustBePositive} = us.seq.numPulse % block size
                kwargs.binary (1,1) logical = false; % whether to use the binary form of k-Wave
                kwargs.isosub (1,1) logical = false; % whether to subtract the background using an isoimpedance sim
                kwargs.gpu (1,1) logical = logical(gpuDeviceCount()) % whether to employ gpu during pre-processing
            end
            arguments % kWaveArray arguments - these are passed to kWaveArray
                karray_args.UpsamplingRate (1,1) double =  10
                karray_args.BLITolerance (1,1) double = 0.05
                karray_args.BLIType (1,1) string {mustBeMember(karray_args.BLIType, ["sinc", "exact"])} = 'sinc' % stencil - exact or sinc
            end
            arguments % kWave 1.1 arguments - these are passed to kWave
                kwave_args.BinaryPath (1,1) string = getkWavePath('binaries')
                kwave_args.CartInterp (1,1) string {mustBeMember(kwave_args.CartInterp, ["linear", "nearest"])}
                kwave_args.CreateLog (1,1) logical
                kwave_args.DataCast (1,1) string = [repmat('gpuArray-',[1,logical(gpuDeviceCount())]), 'single']
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
            arguments
                ksensor_args.directivity_pattern (1,1) string {mustBeMember(ksensor_args.directivity_pattern, ["pressure", "gradient", "total", "none"])} = "none"
            end

            % parse
            if isinf(cgrd.dx), kwargs.el_sub_div(1) = 1; end % don't use sub-elements laterally for 1D sims
            if isinf(cgrd.dy), kwargs.el_sub_div(2) = 1; end % don't use sub-elements in elevation for 2D sims

            % start measuring total execution time
            tt_kwave = tic;

            % get the kWaveGrid
            % TODO: check that the step sizes are all equal - this is
            % required by kWaveArray
            [kgrid, Npml] = getkWaveGrid(cgrd, 'PML', kwargs.PML);

            % get the kWave medium struct
            kmedium = getMediumKWave(med, cgrd);

            % make an iso-impedance medium maybe
            if kwargs.isosub
                kmedium_iso = kmedium;
                kmedium_iso.density = (med.c0 * med.rho0) ./ kmedium.sound_speed;
            else
                kmedium_iso = repmat(kmedium, [1, 0]);
            end

            % get the minimum necessary time step from the CFL
            dt_cfl_max = kwargs.CFL_max * min([cgrd.dz, cgrd.dx, cgrd.dy]) / max(kmedium.sound_speed,[],'all');
            
            % use a time step that aligns with the sampling frequency
            dt_us = inv(us.fs);
            cfl_ratio = dt_cfl_max / dt_us;
        
            if cfl_ratio < 1
                warning("QUPS:kspaceFirstOrder:upsampling",'Upsampling the kWave time step for an acceptable CFL.');
                time_step_ratio = ceil(inv(cfl_ratio));            
            else
                % TODO: do (or don't) downsample? Could be set for temporal reasons
                time_step_ratio = 1;
            end
            kgrid.dt = gather(dt_us / time_step_ratio);
            fs_ = us.fs * time_step_ratio;
            
            % kgrid to cgrd coordinate translation mapping
            kgrid_origin = cellfun(@median, {cgrd.z, cgrd.x, cgrd.y});

            % get the source signal
            txsig = conv(us.tx.impulse, us.seq.pulse, 4*fs_); % transmit waveform, 4x intermediate convolution sampling
            apod =  us.seq.apodization(us.tx); % transmit apodization (M x V)
            del  = -us.seq.delays(us.tx); % transmit delays (M x V)
            txn0 = floor((txsig.t0   + min(del(:))) * fs_); % minimum time sample - must pass through 0
            txne = ceil ((txsig.tend + max(del(:))) * fs_); % maximum time sample - must pass through 0
            t_tx = shiftdim((txn0 : txne)' / fs_, -2); % transmit signal time indices (1x1xT')
            txsamp = permute(apod .* txsig.sample(t_tx - del), [3,1,2]); % transmit waveform (T' x M x V)
            
            % define the source and on the grid 
            % get the direction weights
            [~,~,wnorm] = us.tx.orientations;
            wnorm(1:3,:) = wnorm([3 1 2],:); % k-Wave coordinate mapping
            wnorm = swapdim(wnorm,1,5); % 1 x M x 1 x 1 x 3

            % generate {psig, mask, elem_weights}
            if us.tx == us.rx, aps = "xdc"; else, aps = ["tx", "rx"]; end
            for ap = aps % always do rx last: rx variables used in post
            switch kwargs.ElemMapMethod
                case 'nearest'
                    pg = cgrd.positions(); % -> 3 x Z x X x Y
                    pn = us.(ap).positions; % element positions
                    if kwargs.gpu, [pg, pn] = dealfun(@gpuArray, pg, pn); end
                    [Nx, Ny, Nz] = dealfun(@(n) n + (n==0), kgrid.Nx, kgrid.Ny, kgrid.Nz); % kwave sizing ( 1 if sliced )
                    mask = false(Nx, Ny, Nz); % grid size
                    assert(all(size(mask,1:3) == cgrd.size), 'kWave mask and Scan size do not correspond.');
                    clear ind; % init
                    for n = us.(ap).numel:-1:1 % get nearest pixel for each element
                        ind(n) = argmin(vecnorm(pn(:,n) - pg,2,1),[],'all', 'linear');
                    end
                    assert(numel(unique(ind)) == us.(ap).numel, 'Elements are not unique; the Scan spacing or sizing may be too small.');
                    mask(ind) = true;
                    psig = pagetranspose(txsamp .* wnorm); % -> (J' x T' x V x 1 x 3) with M == J'
                    elem_weights = eye(us.(ap).numel) ; % J' x M with ( J' == M )
                    clear pg;

                case 'linear'
                    % get ambient sound speed
                    c0map = med.c0;

                    % get a mapping of delays and weights to all
                    % (sub-)elements (J' x M)
                    [mask, el_weight, el_dist, el_ind] = us.(ap).elem2grid(cgrd, kwargs.el_sub_div);% perm(X x Y x Z), (J' x M)
                    el_map_grd = sparse((1:nnz(mask))' == el_ind(:)'); % matrix mapping (J' x J'')

                    % apply to transmit signal: for each element
                    if (ap == "tx" || ap == "xdc")
                    [del, apod, t_tx] = dealfun(@(x)shiftdim(x, -1), del, apod, t_tx); % (1 x M x V), (1 x 1 x 1 x T')
                    V = size(del,3); % number of transmits
                    B = min(V,kwargs.bsize); % block size for processing
                    for b = flip(0:B:V-1) % for each block (offset)
                    parfor(v = 1:min(B,V-b)) % B-threaded computation - could be moved to process / remote server % 64-threads / 8 workers
                        [delv, apodv] = deal(sub(del,v+b,3), sub(apod,v+b,3));
                        psigv = swapdim(t_tx,3,4) - swapdim(delv,3,4) - el_dist/c0map; % J''' x M x T' x V (time delay)
                        psigv = sample(txsig, psigv); % sample the time delays (J''' x M x T' x V)
                        psigv = wnorm .* el_weight .* swapdim(apodv,3,4) .* psigv; % per sub-element transmit waveform (J''' x M x T' x V x 3)
                        psigv = reshape(psigv, [prod(size(psigv,1:2)), size(psigv,3:4), 1, size(psigv,5)]); % per element transmit waveform (J'' x T' x V x 1 x 3)
                        psigv = double(psigv); % in-place ?
                        psigv = reshape(el_map_grd * psigv(:,:), [size(el_map_grd,1), size(psigv,2:5)]); % per grid-point transmit waveform (J' x T' x V x 1 x 3)
                        psigv = cast(psigv, 'like', psigv([])); % in-place ?
                        psig(:,:,v+b,:,:) = gather(psigv); % store
                    end
                    end
                    psig = reshape(psig, [size(psig,1:3), 1 3]);
                    end

                case {'karray-direct', 'karray-depend'}
                    karray_args.BLIType = char(karray_args.BLIType);
                    karray_opts = struct2nvpair(karray_args);
                    karray = kWaveArray(us.(ap), kgrid.dim, kgrid_origin, karray_opts{:});
                    mask = karray.getArrayBinaryMask(kgrid);

                    % assign source for each transmission (J' x T' x V)
                    % TODO: abstract this logic to come from transducer directly so
                    % it can be used in fullwave or some other FDTD method
                    switch kwargs.ElemMapMethod
                        case 'karray-direct'
                            % define locally
                            vec = @(x) x(:);

                            % get the offgrid source weights
                            elem_weights = arrayfun(@(i){sparse(vec(karray.getElementGridWeights(kgrid, i)))}, 1:us.(ap).numel);  % (J x {M})
                            elem_weights = cat(2, elem_weights{:}); % (J x M)
                            elem_weights = elem_weights(mask(:),:); % (J' x M)
                            if (ap == "tx" || ap == "xdc")
                                psig = pagemtimes(full(elem_weights), 'none', txsamp, 'transpose'); % (J' x M) x (T' x M x V) -> (J' x T' x V)
                            end

                            % get the offgrid source sizes
                            elem_meas = arrayfun(@(i)karray.elements{i}.measure, 1:us.(ap).numel);
                            elem_dim  = arrayfun(@(i)karray.elements{i}.dim    , 1:us.(ap).numel);
                            elem_norm = elem_meas ./ (kgrid.dx) .^ elem_dim; % normalization factor
                            elem_weights = elem_weights ./ elem_norm; 

                        case 'karray-depend' % compute one at a time and apply casting rules
                            if (ap == "tx" || ap == "xdc")

                            psig = cellfun(@(x) ...
                                {cast(karray.getDistributedSourceSignal(kgrid, x.'), 'like', x)}, ...
                                num2cell(real(txsamp), [1,2]) ...
                                );
                            psig = cat(3, psig{:}); % (J' x T' x V)
                            end
                    end
                    
            end

            if (ap == "rx" || ap == "xdc")
                % define the sensor
                ksensor.mask = mask; % pixels to record

                % set the directivity
                switch ksensor_args.directivity_pattern
                    case {'none'}
                    otherwise
                        switch kwargs.ElemMapMethod
                            case {'linear', 'nearest'}
                                elem_dir = zeros(size(mask)); % receive element directivity
                                theta = deg2rad(us.rx.orientations());
                                if kwargs.ElemMapMethod == "nearest"
                                    elem_dir(ind) = theta; % set element directions
                                else
                                    for n = 1:us.rx.numel, elem_dir(el_ind(:,n)) = theta(n); end % set element directions
                                end

                                ksensor.directivity_angle = elem_dir; % 0 is most sensitive in x
                                ksensor.directivity_size  = us.rx.width; % size of the element for producing the directivity
                                ksensor.directivity_pattern = char(ksensor_args.directivity_pattern); % type of pattern
                            case {'karray-direct', 'karray-depend'}
                                warning('Unable to set directivity for karray element mapping methods.');
                        end
                end
            end

            if (ap == "tx" || ap == "xdc")
                % define the source if it's not empty i.e. all zeros
                switch kwargs.ElemMapMethod
                    case {'nearest', 'linear'} % use vector velocity source
                        if any(sub(psig,1,5),'all'), for v = us.seq.numPulse:-1:1, ksource(v).ux = real(sub(psig,{v,1},[3,5])); end, end % set transmit pulses (J' x T' x V x 1 x 1)
                        if any(sub(psig,2,5),'all'), for v = us.seq.numPulse:-1:1, ksource(v).uy = real(sub(psig,{v,2},[3,5])); end, end % set transmit pulses (J' x T' x V x 1 x 1)
                        if any(sub(psig,3,5),'all'), for v = us.seq.numPulse:-1:1, ksource(v).uz = real(sub(psig,{v,3},[3,5])); end, end % set transmit pulses (J' x T' x V x 1 x 1)
                        [ksource.u_mask] = deal(mask); % set transmit aperture mask
                    case {'karray-direct', 'karray-depend'} % use a pressure source
                        if any(sub(psig,1,5),'all'), for v = us.seq.numPulse:-1:1, ksource(v).p = real(sub(psig,v,3)); end, end % set transmit pulses (J' x T' x V x 1 x 1)
                        [ksource.p_mask] = deal(mask); % set transmit aperture mask
                end

            end

            end

            % set the total simulation time: default to a single round trip at ambient speed
            if isempty(kwargs.T)
                kwargs.T = 2 * (vecnorm(range([cgrd.xb; cgrd.yb; cgrd.zb], 2),2,1) ./ med.c0);
            end
            Nt = 1 + gather(floor((kwargs.T + max(range(t_tx))) / kgrid.dt)); % number of steps in time
            kgrid.setTime(Nt, kgrid.dt);

            % get the receive impulse response function
            rx_imp = copy(us.rx.impulse);
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
            
            % dispatch to the appropriate function
            if kwargs.binary
                if contains(kwave_args.DataCast, "gpuArray")
                    switch kgrid.dim
                        case 1, kspaceFirstOrderND_ = @kspaceFirstOrder1D;
                        case 2, kspaceFirstOrderND_ = @kspaceFirstOrder2DG;
                        case 3, kspaceFirstOrderND_ = @kspaceFirstOrder3DG;
                        otherwise, error("Unsupported dimension size (" + kgrid.dim  +").");
                    end
                else
                    switch kgrid.dim
                        case 1, kspaceFirstOrderND_ = @kspaceFirstOrder1D;
                        case 2, kspaceFirstOrderND_ = @kspaceFirstOrder2DC;
                        case 3, kspaceFirstOrderND_ = @kspaceFirstOrder3DC;
                        otherwise, error("Unsupported dimension size (" + kgrid.dim  +").");
                    end
                end
            else
                switch kgrid.dim
                    case 1, kspaceFirstOrderND_ = @kspaceFirstOrder1D;
                    case 2, kspaceFirstOrderND_ = @kspaceFirstOrder2D;
                    case 3, kspaceFirstOrderND_ = @kspaceFirstOrder3D;
                    otherwise, error("Unsupported dimension size (" + kgrid.dim  +").");
                end
                kwave_args = rmfield(kwave_args, 'BinaryPath');
            end

            if us.seq.numPulse > 1
                % make a unique movie name for each pulse
                [fld, nm, ext] = fileparts(kwave_args.MovieName);        
                mv_nm = cellstr(fullfile(fld, {nm} + "_pulse" + (1:us.seq.numPulse) + ext));
            else
                mv_nm = cellstr(kwave_args.MovieName);
            end
            
            kwave_args = repmat(kwave_args, [us.seq.numPulse, 1]);
            [kwave_args.MovieName] = deal(mv_nm{:});
   
            % get all arguments
            kwave_args_ = arrayfun(...
                @(pulse) struct2nvpair(kwave_args(pulse)), ...
                1:us.seq.numPulse, ... 
                'UniformOutput', false ...
                );

            % add the unspecified arguments
            if ~isempty(varargin), for v = 1:numel(kwave_args_), kwave_args_{v} = cat(1, kwave_args_{v}(:), varargin(:)); end, end

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
                    N = us.rx.numel;
                    assert((rx_imp.tend - rx_imp.t0) > 1/fs_, "Cannot use 'linear' element mapping method when the receiver impulse response function is a Delta function.");
                    rx_sig = gather(real(sample(rx_imp, t_rx(:)' + el_dist(:)/c0map))); % J'' x T''
                    el_map_el = sparse((1:N) == vec(ones(size(el_ind)) .* (1:N)))'; % map from convolved samples to elements
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
                readfun = [];

                % TODO: make reports optional
                fprintf(string(us.seq.type) + " k-Wave simulation completed in %0.3f seconds.\n", toc(tt_kwave));
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
            
            if isstruct(x)
            f = string(fieldnames(x)); % recorded fields
            switch method
                case {} % {'linear', 'nearest'} % particle velocity
                    f = f(ismember(f, {'ux', 'uy', 'uz'})); % velocity fields
                    u = arrayfun(@(f) {x.(f)}, f); % extract each velocity field
                    u = cat(3, u{:}); % J' x T x D, D = 2 or 3
                case {'nearest', 'linear','karray-depend', 'karray-direct'}
                    u = x.p; % pressure (J' x T)
            end
            else % the data is given directly, not in a struct
                u = x; % (J' x T)
            end

            switch method
                case 'karray-depend'
                    [kgrid, karray] = deal(varargin{:});
                    for d = size(u,3):-1:1
                        y(:,:,d) = gather(convn( karray.combineSensorData(kgrid, u(:,:,d)).', rx_sig, 'full'));
                    end
                case {'karray-direct'}
                    elem_weights = deal(varargin{:});
                    for d = size(u,3):-1:1
                        y(:,:,d) = gather(convn(transpose(elem_weights' * double(u(:,:,d))), rx_sig, 'full')); % (T x T' | N) X [(N x J') x (J' x T' | D)]' -> T x N | D
                    end
                case {'nearest'}
                    y = gather(convn(pagetranspose(u), rx_sig, 'full')); % -> (T x N | D)
                case 'linear'
                    [el_weight, el_map_grd, el_map_el] = deal(varargin{:});
                    
                    % create the advanced impulse response function with
                    % which to convolve the output
                    D = size(u,3);
                    el_map_el  = sparse(double(el_map_el ));
                    el_map_grd = sparse(double(el_map_grd));
                    
                    % gather here for OOM issues
                    % TODO: only gather if actually running out of memory
                    usegpu = isa(u, 'gpuArray');
                    u = gather(u);
                    for d = D:-1:1 
                        % [(N x J'') x [[(J'' x J') x (J' x T' | D)] x (T' x T | J'')]]' -> (T x N | D)
                        y(:,:,d) = gather(cast(transpose( ...
                            el_map_el * double(el_weight(:) .* convd( ...
                            gather(cast(el_map_grd' * double(u(:,:,d)), 'like', u)), ...
                            rx_sig, 2, 'full', 'gpu', usegpu ...
                            )) ...
                            ), 'like', u(1)*rx_sig(1)));
                    end

                otherwise, warning('Unrecognized mapping option - mapping to grid pixels by default.');
                   y = gather(convn(u, rx_sig, 'full'));
            end

            % multiply by a normal vector in the direction of each element
            % to combine the velocities
            switch method
                case {}% {'linear', 'nearest'}
                    % y = sum(y .* wnorm, 3); % (T x N x D) * (1 x N x D) -> (T x N)
            end
        end
    
        function chd = kspaceRunSim(kspaceFirstOrderND_, ...
                kgrid, kmedium, kmedium_iso, ksource, ksensor, kwave_args_, ...
                rx_sig, elemmethod, rx_args, t0, fs_, W ...
                )

            % splice
            setgpu = isa(W, 'parallel.Cluster'); % on implicit pools, gpu must be set explicitly
            runiso = ~isempty(kmedium_iso); % if empty, no need to simulate
            Np = numel(ksource); % number of sims
            parfor (puls = 1:Np, W)
                % TODO: make this part of some 'info' logger or something
                fprintf('\nComputing pulse %i of %i\n', puls, Np);
                tt_pulse = tic;

                % configure gpu selection
                if gpuDeviceCount() % if gpus exist
                    % get reference to the desired gpuDevice
                    if setgpu % set the current gpuDevice on implicit pools
                        % TODO: allow input to specify gpu dispatch
                        tsk = getCurrentTask();
                        g = gpuDevice(1+mod(tsk.ID-1, gpuDeviceCount())); % even split across tasks 
                    else % get current device on expliciti pools or currrent session
                        g = gpuDevice();
                    end

                    % specify the device input explicitly if using a GPU binary
                    if endsWith(func2str(kspaceFirstOrderND_), 'G')
                        % explicitly set for kWave binary
                        kwave_args_{puls} = [kwave_args_{puls}, {'DeviceNum'; (g.Index - 1)}]; %#ok<PFOUS>
                    end
                end

                % simulate
                sensor_data     = kspaceFirstOrderND_(kgrid, kmedium    , ksource(puls), ksensor, kwave_args_{puls}{:}); % data is small

                % Post-process the simulation data
                out{puls} = UltrasoundSystem.kspaceFirstOrderPostProc(sensor_data, rx_sig, elemmethod, rx_args{:}); %#ok<PFBNS> data is small

                if runiso
                % simulate
                sensor_data = kspaceFirstOrderND_(kgrid, kmedium_iso, ksource(puls), ksensor, kwave_args_{puls}{:});

                % Post-process the simulation data, and remove background
                out{puls} = out{puls} ...
                          - UltrasoundSystem.kspaceFirstOrderPostProc(sensor_data, rx_sig, elemmethod, rx_args{:}); 
                end

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
        function [b, k, PRE_ARGS, POST_ARGS] = DAS(us, chd, varargin, kwargs)
            % DAS - Delay-and-sum beamformer (compute-optimized)
            %
            % b = DAS(us, chd) performs delay-and-sum beamforming on 
            % the ChannelData chd. The ChannelData must conform to the 
            % delays given by the Sequence us.seq. The output is the
            % image defined on the Scan us.scan. 
            %
            % b = DAS(us, chd, A) uses an apodization defined by the
            % ND-array A.
            %
            % b = DAS(..., 'apod', A) is supported for backwards
            % compatability.
            %
            % b = DAS(us, chd, A1, A2, ..., An) uses multiple separable
            % apodization matrices A1, A2, ... An defined by point-wise
            % multiplication. Each matrix Ai must be broadcastable to size
            % I1 x I2 x I3 x N x M where I1 x I2 x I3 is the size of the
            % image, N is the number of receivers, and M is the number of
            % transmits. In other words, for each Ai,
            % ```
            % assert(all( ...
            %     size(Ai,1:5) == 1 | ...
            %     size(Ai,1:5) == [us.scan.size, us.xdc.numel, us.seq.numPulse] ...
            % ));
            %  ```
            % 
            % NOTE: Multiple apodization matrices are not supported for
            % OpenCL kernels.
            % 
            % b = DAS(us, chd, ..., 'c0', c0) uses a beamforming sound speed
            % of c0. c0 can be a scalar or an ND-array that is broadcastable
            % to size (I1 x I2 x I3) where  [I1, I2, I3] == us.scan.size.
            % The default is us.seq.c0.            
            % 
            % b = DAS(..., 'keep_tx', true) preserves the transmit
            % dimension in the output image b.
            %
            % b = DAS(..., 'keep_rx', true) preserves the receive
            % dimension in the output image b.
            %
            % b = DAS(..., 'fmod', fc) upmixes the data at a modulation
            % frequency fc. This undoes the effect of
            % demodulation/downmixing at the same freuqency.
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
            % [b, k, PRE_ARGS, POST_ARGS] = DAS(...) when the CUDA ptx is
            % used returns the parallel.gpu.CUDAKernel k as well as the
            % arguments for calling the data PRE_ARGS and POST_ARGS. The
            % kernel can then be called per frame f of ChannelData chd as
            %
            %     b{f} = k.feval(PRE_ARGS{:}, chd.data(:,:,:,f), POST_ARGS{:});
            %
            % The ChannelData chd must have the order == 'TNM'. If chd.data
            % is a gpuArray, it must have the same type as was used to
            % create the parallel.gpu.CUDAKernel k. This is useful for
            % processing many identical frames with minimal overhead.
            %
            % NOTE: if the input data is smaller than was used to create
            % the parallel.gpu.CUDAKernel k, an illegal address error may
            % occur, requiring MATLAB to be restarted!
            % 
            % DAS is similar to BFDAS, but is more computationally 
            % efficient at the cost of code readability and available
            % interpolation method because it avoids calling the
            % ChannelData.sample() method.
            %
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % 
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % if isreal(chd) chd = hilbert(chd); end % complex analog
            % chd = zeropad(singleT(chd), 0, max(0, chd.T - 2^9)); % precondition the data
            % b = DAS(us, chd); % beamform the data
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also BFDAS BFADJOINT BFEIKONAL CHANNELDATA.SAMPLE
            
            arguments
                us (1,1) UltrasoundSystem
                chd ChannelData
            end
            arguments (Repeating)
                varargin (:,:,:,:,:) {mustBeNumericOrLogical}
            end
            arguments
                kwargs.c0 (:,:,:,1,1) {mustBeNumeric} = us.seq.c0
                kwargs.fmod (1,1) {mustBeNumeric} = 0 
                kwargs.prec (1,1) string {mustBeMember(kwargs.prec, ["single", "double", "halfT"])} = "single"
                kwargs.device (1,1) {mustBeInteger} = -1 * (logical(gpuDeviceCount()) || (exist('oclDeviceCount','file') && logical(oclDeviceCount())))
                kwargs.apod {mustBeNumericOrLogical} = 1
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.keep_tx (1,1) logical = false
                kwargs.keep_rx (1,1) logical = false
            end

            % parse inputs
            [sumtx, sumrx, c0] = deal(~kwargs.keep_tx, ~kwargs.keep_rx, kwargs.c0);

            % get concatenation dimension(s) for multiple ChannelData
            chds = chd; % ND-array of ChannelData
            chddim = find(size(chds) > 1);
            assert(isempty(chddim) || isscalar(chddim), 'ChannelData array can only contain up to one non-scalar dimension.');
            D = max(max(arrayfun(@(chd)ndims(chd.data), chds)), ndims(chds)); % max dimension of data

            % apodization and modulation
            if ~isequal(1, kwargs.apod), varargin = [varargin, {kwargs.apod}]; end % append if non-default
            apod_args = [[repmat({'apod'},size(varargin)); varargin], {'modulation'; kwargs.fmod}];
            apod_args = apod_args(:)';

            % get positions of the imaging plane
            P_im = us.scan.positions(); % 3 x I1 x I2 x I3 == 3 x [I]

            % get positions of the receive aperture
            P_rx = us.rx.positions(); % 3 x N

            % choose beamforming flag
            if      sumtx &&  sumrx, fun = 'DAS';
            elseif  sumtx && ~sumrx, fun = 'SYN';
            elseif ~sumtx &&  sumrx, fun = 'MUL';
            elseif ~sumtx && ~sumrx, fun = 'BF';
            end

            % process for each ChannelData
            for i = numel(chds):-1:1
                % choose ChannelData
                chd = chds(i);

                % make sure t0 is a scalar in all dims except transmit
                if ~all(size(chd.t0, setdiff(1:ndims(chd.t0), chd.mdim)) == 1), warning("Resampling data for a scalar t0."); chd = rectifyt0(chd); end

                % data must be ordered T x perm(N x M) x ...
                if chd.tdim ~= 1, chd = rectifyDims(chd); end % reorder if necessary

                % get the beamformer arguments
                dat_args = {chd.data, gather(chd.t0), gather(chd.fs), c0, 'device', kwargs.device, 'input-precision', kwargs.prec, 'transpose', (chd.ndim>chd.mdim)}; % data args
                if isfield(kwargs, 'interp'), interp_args = {'interp', kwargs.interp}; else,  interp_args = {}; end
                ext_args = [interp_args, apod_args]; % extra args
    
                switch us.seq.type
                    case 'FSA'
                        [~,~,nf] = us.tx.orientations(); % normal vectors
                        pos_args = {P_im, P_rx, us.tx.positions(), nf};
                        ext_args{end+1} = 'diverging-waves'; %#ok<AGROW>
                    case 'PW'
                        pos_args = {P_im, P_rx, [0;0;0], us.seq.focus}; % TODO: use origin property in tx sequence
                        ext_args{end+1} = 'plane-waves'; %#ok<AGROW> 
                    case {'VS', 'FC', 'DV'}
                        nf = us.seq.focus - us.tx.offset; % normal vector
                        pos_args = {P_im, P_rx, us.seq.focus, nf ./ norm(nf)};
                        if us.seq.type == "DV", ext_args{end+1} = 'diverging-waves'; end %#ok<AGROW>
                end

                % request the CUDA kernel?
                if nargout > 1, ext_ret = cell(1,3); else, ext_ret = {}; end

                % beamform the data (I1 x I2 x I3 x perm(N x M) x F x ...)
                [b, ext_ret{:}] = das_spec(fun, pos_args{:}, dat_args{:}, ext_args{:});

                % move data dimension, back down raise aperture dimensions (I1 x I2 x I3 x F x ... x perm(N x M))
                b = permute(b, [1:3,6:(D+2),4:5]);

                % store
                bi{i} = b;
            end

            % unpack
            if isempty(chddim), b = bi{1}; else, b = cat(chddim, bi{:}); end

            % map outputs
            if nargout > 1, [k, PRE_ARGS, POST_ARGS] = deal(ext_ret{:}); end
        end
        
        function chd = focusTx(us, chd, seq, kwargs)
            % FOCUSTX - Synthesize transmits
            %
            % chd = FOCUSTX(us, chd0) focuses the FSA ChannelData chd0 by
            % linearly synthesizing transmits (i.e. delay and sum across 
            % transmits).
            %
            % chd = FOCUSTX(us, chd0, seq) uses the Sequence seq to focus
            % the data. The default is us.seq.
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
            % chd = FOCUSTX(..., 'bsize', B) uses a maximum block size of B
            % output transmits at a time when vectorizing computations. A 
            % larger block size will run faster, but use more memory. The
            % default is seq.numPulse.
            %
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % 
            % % Compute the data for an FSA acquistion
            % us.seq = Sequence('type', 'FSA', 'c0', us.seq.c0, 'numPulse', us.xdc.numel);
            % chd = greens(us, scat); % compute the response
            %
            % % Create plane-wave data by synthesizing the transmits with plane-wave delays
            % 
            % seq_pw = SequenceRadial('type', 'PW', 'c0', us.seq.c0, ...
            %  'angles', -25:5:25, 'ranges', 1); % plane-wave sequence
            % chd_pw = focusTx(us, chd, seq_pw); % synthesize transmits
            %
            % % Beamform the plane-wave data
            % us.seq = seq_pw;
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
                us (1,1) UltrasoundSystem
                chd (1,1) ChannelData
                seq (1,1) Sequence = us.seq;
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.length string {mustBeScalarOrEmpty} = string.empty;
                kwargs.buffer (1,1) {mustBeNumeric, mustBeInteger} = 0
                kwargs.bsize (1,1) {mustBePositive, mustBeInteger} = max(1,seq.numPulse*us.tx.numel)
                kwargs.verbose (1,1) logical = false;
            end

            % Copy semantics
            chd = copy(chd);

            % dist/time to receiver
            tau = - seq.delays(     us.tx); %    M  x    M'
            apd =   seq.apodization(us.tx); % [1|M] x [1|M']

            % nothing to do for (true) FSA acquisitions: all 0 delays
            % identity matrix apodization
            switch seq.type, case 'FSA' 
                if ~nnz(tau) && isequal(apd, eye(us.tx.numel)), return; end 
            end

            % resample only within the window where we currently have data.
            i = logical(apd) | false(size(tau)); % non-zero apodization indices (broadcasted)
            nmin = floor(min(tau(i),[],'all','omitnan') .* chd.fs); % minimum sample time
            nmax =  ceil(max(tau(i),[],'all','omitnan') .* chd.fs); % maximum sample time
            chd.t0 = chd.t0 + nmin / chd.fs; % shift time axes forwards to meet minimum sample time
            tau = tau - nmin / chd.fs; % shift delays backwards to meet time axes
            chd = zeropad(chd, 0, (nmax - nmin) + kwargs.buffer); % expand time axes to capture all data
            
            % legacy: pick new signal length
            L = kwargs.length;
            if kwargs.interp == "freq"
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
            D = 1+max([3,ndims(chd.data), ndims(chd.t0)]); % get a free dimension for M'
            assert(D > chd.mdim, "Transmit must be in the first 3 dimensions (" + chd.mdim + ").");
            tau = swapdim(tau,[1 2],[chd.mdim D]); % move data
            apd = swapdim(apd,[1 2],[chd.mdim D]); % move data
            
            % splicing
            B = max(1, floor(kwargs.bsize / us.tx.numel)); % simultaneous blocks of TxN, as best we can manage
            id  = num2cell((1:B)' + (0:B:seq.numPulse-1),1); % output sequence indices
            id{end}(id{end} > seq.numPulse) = []; % delete OOB indices

            if kwargs.verbose
                disp("Focusing with "+numel(id)+" blocks of "+B+" transmit(s) each.")
            end

            % sample and store
            for i = numel(id):-1:1
                z{i} = chd.sample2sep(chd.time, - sub(tau,id{i},D), kwargs.interp, sub(apd,id{i},D), chd.mdim); % sample ({perm(T' x N x 1) x F x ...} x M')
                z{i} = swapdim(z{i}, chd.mdim, D); % replace transmit dimension (perm({T' x N} x M') x {F x ...})
            end
            z = cat(chd.mdim, z{:}); % unpack transmit dimension (perm(T' x N x M') x F x ...)
            chd.data = z; % store output channel data % (perm(T' x N x M') x F x ...)
        end
        
        function [chd, Hi] = refocus(us, chd, seq, kwargs)
            % REFOCUS - Recreate full-synthetic aperture data
            %
            % chd = refocus(us, chd) refocuses the ChannelData chd
            % captured with the UltrasoundSystem us into FSA data.
            %
            % [chd, Hi] = refocus(...) additionally returns the decoding
            % matrix Hi. Hi is of size (V x M x T) where V is the number
            % of transmit pulses, M is the number of transmit elements, and
            % T is the number of time samples.
            % 
            % chd = refocus(us, chd, seq) refocuses the ChannelData chd
            % assuming that the ChannelData chd was captured using Sequence
            % seq instead of the Sequence us.seq.
            %
            % [...] = refocus(..., 'method', 'tikhonov') uses a tikhonov
            % inversion scheme to compute the decoding matrix. This is
            % best employed for focused transmit seuqneces
            %
            % [...] = refocus(..., 'method', 'tikhonov', 'gamma', gamma)
            % uses a regularization parameter of gamma for the tikhonov
            % inversion. If the method is not 'tikhonov', the argument is
            % ignored. The regularization  parameter is scaled by the
            % maximum singular value for each frequency. The default is
            % (chd.N / 10) ^ 2.
            %
            % References:
            % [1] Ali, R.; Herickhoff C.D.; Hyun D.; Dahl, J.J.; Bottenus, N. 
            % "Extending Retrospective Encoding For Robust Recovery of the Multistatic Dataset". 
            % IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
            % vol. 67, no. 5, pp. 943-956, Dec. 2019.
            % <a href="matlab:web('https://doi.org/10.1109/TUFFC.2019.2961875')">DOI 10.1109/TUFFC.2019.2961875</a>
            % 
            % [2] Bottenus, N. "Recovery of the complete data set from focused transmit beams".
            % IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
            % vol. 65, no. 1, pp. 30-38, Jan. 2018.
            % <a href="matlab:web('https://doi.org/10.1109/TUFFC.2017.2773495')">DOI 10.1109/TUFFC.2017.2773495</a>
            % 
            % [3] Hyun, D; Dahl, JJ; Bottenus, N. 
            % "Real-Time Universal Synthetic Transmit Aperture Beamforming with 
            % Retrospective Encoding for Conventional Ultrasound Sequences (REFoCUS)".
            % 2021 IEEE International Ultrasonics Symposium (IUS), pp. 1-4.
            % <a href="matlab:web('https://doi.org/10.1109/IUS52206.2021.9593648')">DOI: 10.1109/IUS52206.2021.9593648</a>
            % 
            % 
            % Example:
            % % Define the setup - make plane waves
            % c0 = 1500;
            % us = UltrasoundSystem();
            % seq = Sequence(        'type', 'FSA', 'numPulse', us.xdc.numel, 'c0', c0);
            % seqpw = SequenceRadial('type', 'PW', 'angles', -45:1.5:45, 'c0', c0);
            % seqfc = SequenceRadial('type', 'FC', 'focus',[0 0 20]' + [1 0 0]'.*sub(us.xdc.positions,1,1), 'c0', c0);
            % scat = Scatterers('pos', [5 0 30]'*1e-3, 'c0', seqpw.c0); % define a point target
            %
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % b = DAS(us, chd); % beamform the data
            %
            % % Focus into plane waves and beamform
            % chdpw = focusTx(us, chd, seqpw);
            % uspw = copy(us);
            % uspw.seq = seqpw;
            % bpw = DAS(uspw, chdpw);
            %
            % % Refocus back to FSA, and beamform
            % chdfsa1 = refocus(us, chdpw, seqpw);
            % bfsa1 = DAS(us, chdfsa1);
            %
            % % Focus at focal points and beamform
            % chdfc = focusTx(us, chd, seqfc);
            % usfc = copy(us);
            % usfc.seq = seqfc;
            % bfc = DAS(usfc, chdfc);
            %
            % % Refocus back to FSA, and beamform
            % chdfsa2 = refocus(us, chdpw, seqpw);
            % bfsa2 = DAS(us, chdfsa2);
            %
            % % Display the channel data
            % figure('Name', 'Channel Data');
            % tiledlayout(3,2);
            % chds = [chd, chd, chdpw, chdfsa1, chdfc, chdfsa2];
            % tnms = ["Original", "Original", "PW", "PW-REF", "FC", "FC-REF"];
            % for i = 1:numel(chds)
            %     himc = imagesc(chds(i),ceil(chds(i).M/2),nexttile()); 
            %     title(tnms(i));
            %     colorbar; colormap default; caxis(max(caxis) + [-50 0]);
            % end
            % linkaxes([himc.Parent]);
            %
            % % Display the images
            % figure('Name', 'B-mode');
            % tiledlayout(3,2);
            % bs = {b, b, bpw, bfsa1, bfc, bfsa2};
            % bpow = gather(mod2db(cellfun(@(b)max(b,[],'all'), bs)));
            % for i = 1:numel(bs)
            %     himb(i) = imagesc(us.scan, bs{i}, nexttile(), bpow(i) + [-80 0]);
            %     colormap gray; colorbar; 
            %     title(tnms(i));
            % end
            % linkaxes([himb.Parent]);
            %
            % See also FOCUSTX
            arguments
                us (1,1) UltrasoundSystem
                chd (1,1) ChannelData
                seq (1,1) Sequence = us.seq
                kwargs.gamma (1,1) {mustBeNumeric} = (chd.N / 10)^2 % heuristically chosen
                kwargs.method (1,1) string {mustBeMember(kwargs.method, "tikhonov")} = "tikhonov"
            end

            % dispatch pagenorm/pagemrdivide function based on MATLAB version
            if     exist('pagenorm', 'file')
                pagenorm2 = @(x) pagenorm(x,2);
            elseif exist('pagesvd', 'builtin')
                pagenorm2 = @(x) sub(pagesvd(x),{1,1},1:2);
            else
                pagenorm2 = @(x) cellfun(@(x) svds(x,1), num2cell(double(x),1:2));
            end
            if exist('pagemrdivide','builtin')
                pagemrdivide_ = @pagemrdivide;
            else
                pagemrdivide_ = @(x,y) cell2mat(cellfun(@mrdivide,num2cell(x,1:2),num2cell(y,1:2),'UniformOutput',false));
            end

            % get the apodization / delays from the sequence
            tau = -seq.delays(     us.tx); % (M x V)
            a   =  seq.apodization(us.tx); % (M x V)

            % get the frequency vectors
            f = gather(chd.fftaxis); % perm(... x T x ...)

            % construct the encoding matrix (M x V x T)
            H = a .* exp(+2j*pi*swapdim(f,chd.tdim,3).*tau);

            % compute the pagewise tikhonov-inversion inverse (V x M x T)
            % TODO: there are other options according to the paper - they
            % should be options here
            switch kwargs.method
                case "tikhonov"
                    % TODO: option to use pinv, as it is (slightly)
                    % different than matrix division
                    A = real(pagemtimes(H, 'ctranspose', H, 'none')) + (kwargs.gamma * pagenorm2(gather(H)).^2 .* eye(chd.M)); % A = (H'*H + gamma * I)
                    Hi = pagetranspose(pagemrdivide_(gather(H), gather(A))); % Hi = (A^-1 * H)' <=> (H / A)'
                    Hi = cast(Hi, 'like', H);
            end

            % move inverse matrix to matching data dimensions
            D = max(ndims(chd.data), ndims(chd.t0));
            ord = [chd.mdim, D+1, chd.tdim];
            ord = [ord, setdiff(1:D, ord)]; % all dimensions
            Hi = ipermute(Hi, ord);

            % move data to the frequency domain
            x = chd.data;
            x = fft(x,chd.T,chd.tdim); % get the fft in time
            omega0 = exp(-2i*pi*f.*chd.t0); % time-alignment phase
            x = x .* omega0; % phase shift to re-align time axis
            
            % apply to the data - this is really a tensor-times-matrix
            % operation, but it's not natively supported yet. 
            for v = chd.N:-1:1
                y{v} = sum(sub(Hi,v,D+1) .* x, chd.mdim);
            end
            y = cat(chd.mdim, y{:});

            % move back to the time domain
            t0 = min(chd.t0,[],'all');
            y = y .* exp(+2i*pi*f.*t0); % re-align time axes
            y = ifft(y, chd.T, chd.tdim);
            
            % copy semantics
            chd = copy(chd);
            chd.data = y;
            chd.t0 = t0;
        end

        function [b, tau_rx, tau_tx, tau_foc] = bfAdjoint(us, chd, varargin, kwargs)
            % BFADJOINT - Adjoint method beamformer
            %
            % b = BFADJOINT(us, chd) beamforms the ChannelData chd using
            % an adjoint matrix method. This method computes the
            % inner-product of the normalized transmitted wave with the
            % received data in the frequency domain.
            % 
            % b = BFADJOINT(us, chd, A) uses an apodization 
            % defined by the ND-array A.
            %
            % b = BFADJOINT(..., 'apod', A) is supported for backwards
            % compatability.
            %
            % b = BFADJOINT(us, chd, A1, A2, ..., An) uses multiple
            % separable apodization matrices A1, A2, ... An defined by
            % point-wise multiplication. Each matrix Ai must be
            % broadcastable to size I1 x I2 x I3 x N x M where I1 x I2 x I3 
            % is the size of the image, N is the number of receivers, and M 
            % is the number of transmits. In other words, for each Ai,
            % ```
            % assert(all( ...
            %     size(Ai,1:5) == 1 | ...
            %     size(Ai,1:5) == [us.scan.size, us.xdc.numel, us.seq.numPulse] ...
            % ));
            %  ```
            % 
            % b = BFADJOINT(..., 'c0', c0) uses a beamforming sound speed 
            % of c0. c0 can be a scalar or an ND-array that is broadcastable
            % to size (I1 x I2 x I3) where [I1, I2, I3] == us.scan.size. 
            % The default is us.seq.c0.
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
            % demodulation/downmixing at the same frequency.
            %
            % b = BFADJOINT(..., 'bsize', B) uses a maximum block size of B
            % when vectorizing computations. A larger block size will run 
            % faster, but use more memory. The default is chosen
            % heuristically.
            %
            % b = BFADJOINT(..., 'fthresh', thresh) neglects frequencies
            % where the maximum power across the apertures is less than 
            % thresh dB-down with respect to the maximum power. This can 
            % accelerate computation but with a loss of accuracy. The
            % default is -Inf.
            %
            % b = BFADJOINT(..., 'parenv', clu) or
            % b = BFADJOINT(..., 'parenv', pool) uses the
            % parallel.Cluster clu or the parallel.Pool pool to
            % parallelize computations. parallel.ThreadPools will be
            % ignored due to mex function restrictions.
            % 
            % [b, tau_rx] = BFADJOINT(...) additionally returns the receive
            % delays used in the receive propagation step.
            % 
            % [b, tau_rx, tau_tx, tau_foc] = BFADJOINT(...) additionally 
            % returns the transmitter delays and focal delays used in the 
            % transmit propagation step. 
            % 
            % The ND-array tau_tx represents the propagation delay from 
            % element to pixel. The matrix tau_foc represents the transmit
            % element delays per transmit of the pulse sequence. Note that
            % these have a different interpretation than in bfDAS due to
            % how the adjoint method is defined. 
            % 
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % 
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % chd = zeropad(singleT(chd), 0, max(0, chd.T - 512)); % precondition the data
            % if isreal(chd), chd = hilbert(chd); end % complex analog
            % b = bfAdjoint(us, chd); % beamform the data
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also BFEIKONAL BFDAS DAS FOCUSTX

            % TODO: test for tall types where receive dimension is tall -
            % should work ...
            arguments
                us (1,1) UltrasoundSystem
                chd (1,1) ChannelData
            end
            arguments(Repeating)
                varargin (:,:,:,:,:) {mustBeNumericOrLogical} % apodization matrices (I1 x I2 x I3 x N x M)
            end
            arguments
                kwargs.c0 (:,:,:,1,1) {mustBeNumeric} = us.seq.c0
                kwargs.fmod (1,1) {mustBeNumeric, mustBeFinite} = 0 % modulation frequency
                kwargs.fthresh (1,1) {mustBeReal, mustBeNegative} = -Inf; % threshold for including frequencies
                kwargs.apod {mustBeNumericOrLogical} = 1; % apodization matrix (I1 x I2 x I3 x N x M)
                kwargs.Nfft (1,1) {mustBeInteger, mustBePositive} = chd.T; % FFT-length
                kwargs.keep_tx (1,1) logical = false % whether to preserve transmit dimension
                kwargs.keep_rx (1,1) logical = false % whether to preserve receive dimension
                kwargs.bsize (1,1) {mustBeFloat, mustBeInteger, mustBePositive} = max(1,floor(1*(2^30 / (4*chd.N*us.scan.nPix*8)))); % vector computation block size
                kwargs.verbose (1,1) logical = true 
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])}
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
            if ~isa(us.tx, 'TransducerArray') && ~isa(us.rx, 'TransducerArray') && any(us.seq.type == ["DV","FC","VS"])
                warning('QUPS:bfAdjoint:UnsupportedSequence', ...
                    'This function is unsupported for focused transmits with non-linear transducers.' ...
                    );
            end

            % set the default parallel environment
            if ~isfield(kwargs, 'parenv')
                % use a thread pool if activated and no gpu variables in use
                if isa(gcp('nocreate'), 'parallel.ThreadPool') && ~isa(chd.data, 'gpuArray')
                    kwargs.parenv = gcp('nocreate');
                else 
                    kwargs.parenv = 0;
                end
            end

            % parse inputs
            sumtx = ~kwargs.keep_tx;
            sumrx = ~kwargs.keep_rx;
            fmod = kwargs.fmod;
            kvb = kwargs.verbose;
            varargin = [{kwargs.apod}, varargin]; % backwards compatibility

            % move the data to the frequency domain, being careful to
            % preserve the time axis
            K = kwargs.Nfft; % DFT length
            f = chd.fftaxis; % frequency axis
            df = chd.fs / double(K); % frequency step size
            x = chd.data; % reference the data perm(K x N x M) x ...
            x = x .* exp(+2i*pi*fmod .* chd.time); % remodulate data
            x = fft(x,K,chd.tdim); % get the fft
            x = x .* exp(-2i*pi*f    .* chd.t0  ); % phase shift to re-align time axis
            x = x .* exp(+2i*pi*f .* shiftdim(us.seq.t0Offset(),2-chd.mdim)); % offset to transducer origin

            % choose frequencies to evaluate
            xmax = max(x, [], chd.tdim); % maximum value per trace
            f_val = mod2db(x) - mod2db(xmax) >= kwargs.fthresh; % threshold
            f_val = f_val & f < chd.fs / 2; % positive frequencies only
            f_val = any(f_val, setdiff(1:ndims(x), chd.tdim)); % evaluate only freqs across aperture/frames that is above threshold

            % get the pixel positions
            D = gather(max([4, ndims(chd.data), ndims(chd.t0)])); % >= 4
            Pi = us.scan.positions();
            Pi = swapdim(Pi, 1:4, [1, D+(1:3)]); % place I after data dims (3 x 1 x 1 x 1 x ... x [I])
            c0 = shiftdim(kwargs.c0, -D); % 1 x 1 x 1 x 1 x ... x [I]

            % get the receive apodization, spliced if it can be applied
            [a_n, a_m, a_mn] = deal(1);
            for s = 1:numel(varargin)
            apod = varargin{s};
            if all(size(apod, 1:3) == 1) % image is scalar, apply to data
                a_mn = a_mn .* shiftdim(apod, 3); % N x V
            elseif size(apod, 4) == 1 % receive is scalar, apply over tx
                ord = [D+(1:4), 2]; % send dims 1-5 here
                ord = [ord, setdiff(1:max(ord, 5), ord)]; %#ok<AGROW> % complete set of dimensions
                a_m = a_m .* ipermute(apod, ord); % 1 x V x 1 x 1 x ... x [I] x 1
            elseif size(apod, 5) == 1 % transmit is scalar, apply over rx
                ord = [D+(1:3), 2]; % send dims 1-4 here
                ord = [ord, setdiff(1:max(ord, 5), ord)]; %#ok<AGROW> % complete set of dimensions
                a_n = a_n .* ipermute(apod, ord); % 1 x N x 1 x 1 x ... x [I]
            else % none of the following is true: this request is excluded for now
                error("Unable to apply apodization ("+s+") due to size constraints. Apodization must be scalar in the transmit dimension, receive dimension, or all image dimensions.")
            end
            end
            
            % get the delays for the transmit/receive green's matrix
            % kernels
            tau_tx = vecnorm(us.tx.positions() - Pi,2,1) ./ c0; % 1 x M x 1 x 1 x ... x [I]
            tau_rx = vecnorm(us.rx.positions() - Pi,2,1) ./ c0; % 1 x N x 1 x 1 x ... x [I]

            % get the transmit steering vector weights and delays
            del_tx  = us.seq.delays(us.tx);      % M x V
            apod_tx = us.seq.apodization(us.tx); % M x V
            del_tx = del_tx + us.seq.t0Offset(); % offset to transducer origin

            % transform to frequency step kernels
            w_rx    = exp(-2i*pi*df.*tau_rx); %  receive greens function
            w_tx    = exp(-2i*pi*df.*tau_tx); % transmit greens function
            w_steer = exp(-2i*pi*df.*del_tx); % transmit steering delays

            % cast data type for efficency
            [w_tx, w_rx, w_steer, apod_tx] = dealfun(@(w) cast(w, 'like', real(x)), w_tx, w_rx, w_steer, apod_tx);
            b = zeros(1,'like', x);

            % splice data 
            k = gather(find(f_val)'); % skip unimportant frequencies
            k = num2cell(reshape([k(:); nan([-mod(numel(k),-kwargs.bsize),1])], kwargs.bsize, []), 1)'; % expand and make blocks
            k{end}(isnan(k{end})) = []; % delete invalid entries
            xk = cellfun(@(k) sub(x,k,chd.tdim), k, 'UniformOutput',false);
            chd_ord = [chd.ndim, chd.mdim, chd.tdim]; % permution order
            ichd_ord(chd_ord) = 1:3; % inverse permutation

            if kwargs.verbose, tt = tic; fprintf("Beamforming for " + numel(k) + " frequency bands "); end % . Completed: ["); end
            % DEBUG: plot
            % figure; h = imagesc(squeeze(zeros(us.scan.size))); colorbar; colormap jet;

            % TODO: parallelize for thread-based pools only
            parfor (ik = 1:numel(k), kwargs.parenv)
            % for ik = 1:numel(k)
                % get discrete frequency index - must be float due to pow
                k_ = shiftdim(k{ik}, -2); % 1 x 1 x F

                % report progress
                % lbl = "Beamforming freqs: " + min(gather(f(k_))) + " - " + max(gather(f(k_)));% + " MHz";
                % fprintf(lbl + "\n");

                % data, in freq. domain (N x V x ...)
                % TODO: adapt for tall types
                xk_ = permute(xk{ik}, [chd_ord, 4:D]);
                
                % compute the greens functions on transmit/receive
                % G_tx = w_tx.^(k_-1); % 1 x M x F x 1 x [I]
                % G_rx = w_rx.^(k_-1); % 1 x N x F x 1 x [I]

                % compute the inverse steering vector on transmit
                A_tx = apod_tx .* w_steer.^(k_-1); % M x V x F
                A_tx = pagemtimes(w_tx.^(k_-1), A_tx); % 1 x V x F x 1 x [I]
                Ainv_tx = pagetranspose(A_tx); % V x 1 x F x 1 x [I] % make a column vector
                Ainv_tx = Ainv_tx ./ vecnorm(Ainv_tx, 2, 1); % normalize the power
                
                % apodize, delay, and sum the data for this frequency
                % only 1 of the a_* will contain the apodization
                if sumrx, yn = a_m .*      pagemtimes(a_n .* conj(w_rx.^(k_-1))   , a_mn .* xk_); % 1 x V x F x 1 x [I]
                else,     yn = a_m .* (pagectranspose(a_n .*     (w_rx.^(k_-1))) .* a_mn .* xk_); % N x V x F x 1 x [I]
                end
                if sumtx, y  = pagemtimes(yn,                 conj(Ainv_tx));  % [1|N] x 1 x F x 1 x [I]
                else,     y  =           (yn .* pagectranspose(   (Ainv_tx))); % [1|N] x V x F x 1 x [I]
                end

                % integrate over all frequencies
                b = b + sum(y,3); % [1|N] x [1|V] x 1 x 1 x [I]

                % DEBUG: update display
                % h.CData(:) = mod2db(sum(y,3)); drawnow limitrate; 
            
                % if kvb, fprintf(string(ik) +","); end
                if kvb, fprintf("."); end
            end
            % if kvb, fprintf("\b]\nDone! "); toc(tt), end
            if kvb, fprintf("\nDone! "); toc(tt), end

            % move to image dimensions ([I] x ... x perm([1|N] x [1|V] x 1)
            b = swapdim(b, [D+(1:3), 3+(1:D-3), ichd_ord], 1:D+3);

            % set delays to proper output sizing
            if nargout >= 2
                tau_rx = permute(reshape(tau_rx, [us.rx.numel, us.scan.size]), [2:4,1,5]);
                tau_tx = permute(reshape(tau_tx, [us.tx.numel, us.scan.size]), [2:4,5,1]);
                tau_foc = del_tx; 
            end
        end    

        function [b, tau_rx, tau_tx] = bfEikonal(us, chd, med, cgrd, varargin, kwargs)
            % BFEIKONAL - Delay-and-sum beamformer with Eikonal delays
            %
            % b = BFEIKONAL(us, chd, med, cgrd) creates a b-mode
            % image b from the ChannelData chd and Medium med using the
            % delays given by the solution to the eikonal equation defined
            % on the ScanCartesian cgrd.
            % 
            % The transmitter and receiver must fall within the cgrd. The 
            % step size in each dimension must be identical. The eikonal 
            % equation is solved via the fast marching method.
            %
            % b = BFEIKONAL(us, chd, med) uses the ScanCartesian us.scan as
            % sound speed grid for computing time-delays.
            % 
            % The grid spacing for each dimension must be (almost)
            % identical e.g. assuming cgrd.y = 0, 
            % abs(cgrd.dx - cgrd.dz) < eps must hold.
            % 
            % A good heuristic is a grid spacing of < lambda / 10 to 
            % avoid accumulated phase errors.
            % 
            % b = BFEIKONAL(us, chd, med, cgrd, A) uses an apodization 
            % defined by the ND-array A.
            %
            % b = BFEIKONAL(..., 'apod', A) is supported for backwards
            % compatability.
            %
            % b = BFEIKONAL(us, chd, med, cgrd, A1, A2, ..., An) uses 
            % multiple separable apodization matrices A1, A2, ... An
            % defined by point-wise multiplication. Each matrix Ai must be
            % broadcastable to size I1 x I2 x I3 x N x M where I1 x I2 x I3
            % is the size of the image, N is the number of receivers, and M 
            % is the number of transmits. In other words, for each Ai,
            % ```
            % assert(all( ...
            %     size(Ai,1:5) == 1 | ...
            %     size(Ai,1:5) == [us.scan.size, us.xdc.numel, us.seq.numPulse] ...
            % ));
            %  ```
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
            % demodulation/downmixing at the same frequency.
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
            % b = BFEIKONAL(..., 'bsize', B) uses a maximum block size of B
            % when vectorizing computations. A larger block size will run 
            % faster, but use more memory. The default is chosen by bfDASLUT.
            %
            % [b, tau_rx, tau_tx] = BFEIKONAL(...) additionally returns the
            % receive and transmit time delays tau_rx and tau_tx are size
            % (I1 x I2 x I3 x N x 1) and (I1 x I2 x I3 x 1 x M), where 
            % I1 x I2 x I3 is the size of the image, N is the number of
            % receivers, and M is the number of transmits.
            % 
            % [...] = BFEIKONAL(..., 'delay_only',true) computes delays but
            % avoids computing the image.
            % 
            % References: 
            % [1] Hassouna MS, Farag AA. 
            % Multi-stencils fast marching methods: a highly accurate solution to the eikonal equation on cartesian domains.
            % IEEE Trans Pattern Anal Mach Intell. 2007 Sep;29(9):1563-74. 
            % <a href="matlab:web('https://doi.org/10.1109/tpami.2007.1154')">DOI: 10.1109/TPAMI.2007.1154</a>. PMID: 17627044.
            % 
            % Example:
            % % This example requires kWave
            % if ~exist('kWaveGrid', 'class')
            %     warning('kWave must be on the path to run this example.');
            % end
            % 
            % % Setup a system
            % cgrd = ScanCartesian(...
            %   'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
            %   'z', 1e-3*linspace(-02, 58, 1+60*2^3) ...
            % );
            % xdc = TransducerArray('numel', 16, 'fc', 3e6, 'bw', [1.5, 4.5]*1e6, 'pitch', 1.5e3/3e6);
            % seq = Sequence('type', 'FSA', 'numPulse', xdc.numel, 'c0', 1500);
            % us = UltrasoundSystem('scan', cgrd, 'xdc', xdc, 'seq', seq, 'fs', 10*xdc.fc);
            % 
            % % Create a Medium to simulate
            % [c0, rho0] = deal(1.5e3, 1e3); 
            % [c, rho] = deal(c0*ones(cgrd.size), rho0*ones(cgrd.size));
            % [Xg, ~, Zg] = cgrd.getImagingGrid();
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
            % sct = Scatterers('pos', [0 0 1e-3]'*(10 : 10 : 50));
            % 
            % % Construct the Medium
            % med = Medium.Sampled(cgrd, c, rho);
            % 
            % % Simulate the ChannelData
            % us.fs = single(us.fs); % accelerate
            % chd = kspaceFirstOrder(us, med, cgrd, 'CFL_max', 0.50);
            % chd = hilbert(chd); % get analytic signal
            % 
            % % Beamform
            % b_naive = DAS(us, chd);
            % b_c0    = bfEikonal(us, chd, med, cgrd);
            % 
            % % Display the ChannelData
            % figure; imagesc(chd); dbr echo; title('Channel Data');
            % 
            % % Display the images
            % figure;
            % h    = nexttile(); imagesc(us.scan, b_naive); title(  'Naive Delay-and-Sum');
            % dbr b-mode; hold on; plot(us.xdc, 'c'); plot(sct, 'r.'); 
            % h(2) = nexttile(); imagesc(us.scan, b_c0   ); title('Eikonal Delay-and-Sum');
            % dbr b-mode; hold on; plot(us.xdc, 'c'); plot(sct, 'r.'); 
            %
            % % formatting
            % linkprop(h, 'CLim');
            % linkaxes(h);
            % arrayfun(@legend, h)
            % 
            % See also DAS BFDASLUT BFDAS BFADJOINT

            arguments
                us (1,1) UltrasoundSystem
                chd ChannelData = ChannelData.empty
                med (1,1) Medium = Medium('c0', us.seq.c0)
                cgrd (1,1) ScanCartesian = us.scan
            end
            arguments(Repeating)
                varargin (:,:,:,:,:) {mustBeNumericOrLogical}
            end
            arguments
                kwargs.fmod (1,1) {mustBeNumeric} = 0
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                kwargs.parenv {mustBeScalarOrEmpty, mustBeA(kwargs.parenv, ["parallel.Cluster", "parallel.Pool", "double"])} = gcp('nocreate')
                kwargs.apod {mustBeNumericOrLogical} = 1
                kwargs.keep_rx (1,1) logical = false;
                kwargs.keep_tx (1,1) logical = false;
                kwargs.verbose (1,1) logical = true;
                kwargs.bsize (1,1) {mustBePositive, mustBeInteger}
                kwargs.delay_only (1,1) logical = isempty(chd); % compute only delays
            end

            % oops?
            if kwargs.delay_only && (nargout <= 1)
                warning( ...
                    "QUPS:bfEikonal:unexpectedOutput", ...
                    "Computing delays only, but only " + nargout + " outputs were requested. Was this intentional?" ...
                    )
            elseif ~kwargs.delay_only && isempty(chd)
                warning( ...
                    "QUPS:bfEikonal:emptyChannelData", ...
                    "An image was requested , but the ChannelData is empty. Was this intentional?" ...
                    )
            end

            % get dimensions of ChannelData
            assert(sum(size(chd) > 1) <= 1, "ChannelData array may only have 1 non-singleton dimension (" + join(string([size(chd)])," x ") + ").");

            % check the data is FSA and matches transmit / receive
            % transducers
            if ~kwargs.delay_only
                assert(all(us.tx.numel == [chd.M]), 'Number of transmits must match number of transmitter elements.')
                assert(all(us.rx.numel == [chd.N]), 'Number of receives must match number of receiver elements.')
            end

            % get cluster
            parenv = kwargs.parenv; % compute cluster/pool/threads
            % travel times cannot use threadPool with mex call
            if isempty(parenv) || (endsWith(which("msfm"+nnz(cgrd.size>1)+"d"),mexext()) && (isa(parenv, 'parallel.ThreadPool') || isa(parenv, 'parallel.BackgroundPool'))), parenv = 0; end

            % get worker transfer function - mark constant to save memory
            if isa(parenv, 'parallel.ProcessPool') && isscalar(gcp('nocreate')) && (parenv == gcp('nocreate')) 
                constfun = @parallel.pool.Constant; % only use if parenv is the current process pool
            else 
                constfun = @(x) struct('Value', x);
            end

            % get the grid definition
            dp = [cgrd.dx, cgrd.dy, cgrd.dz]; % step size each dim
            grd = {cgrd.x, cgrd.y, cgrd.z}; % grid definition
            nsdims = (~isinf(dp)); % non-singleton dimensions
            dp = uniquetol(dp(nsdims)); % spatial step size
            ord = argn(2, @ismember, cgrd.order, 'XYZ'); % order of grid
            sel = ord(ismember(ord, find(nsdims))); % selection of grid indices
            grd = grd(sel); % re-order and trim the grid
            
            % check the grid spacing is identical in non-singleton 
            % dimensions
            assert(numel(dp) == 1, ...
                "The simulation scan must have equally sized steps in all non-singleton dimensions." ...
                );

            % get the transmit, receive
            Pv = us.tx.positions();
            Pr = us.rx.positions();
            
            % get the sound speed in the Medium
            c = props(med, cgrd, 'c');

            % convert positions to sound speed grid coordinates (1-based)
            og = [cgrd.x(1); cgrd.y(1); cgrd.z(1)];
            [Pvc, Prc] = dealfun(@(x) sub((x - og) ./ dp, sel, 1) + 1, Pv, Pr); % ([2|3] x [N|M])

            % enforce type and send data to workers maybe
            cnorm = constfun(double(gather(c ./ dp)));

            % get one-way delays within the field then generate samplers, using
            % reduced dimensions ([M|N] x {Cx x Cy x Cz})
            gi_opts = {'cubic', 'none'}; % interpolater options
            if kwargs.verbose 
                tt = tic; fprintf('\nComputing Eikonal time delays ... \n');
            end
            parfor (n = 1:us.rx.numel, parenv)
                % fprintf('rx %i\n', n);
                [tau_map_rx] = msfm(squeeze(cnorm.Value), double(Prc(:,n))); %#ok<PFBNS> % travel time to each point
                rx_samp{n} = griddedInterpolant(grd, tau_map_rx,gi_opts{:}); %#ok<PFBNS> % make interpolator on cgrd
            end
            if us.tx == us.rx % if apertures are identical, copy
                tx_samp = rx_samp;
            else % else compute for each tx
            parfor (m = 1:us.tx.numel, parenv)
                % fprintf('tx %i\n', m);
                [tau_map_tx] = msfm(squeeze(cnorm.Value), double(Pvc(:,m))); %#ok<PFBNS> % travel time to each point
                tx_samp{m} = griddedInterpolant(grd, tau_map_tx,gi_opts{:}); %#ok<PFBNS> % make interpolator on cgrd
            end
            end
            if kwargs.verbose
                fprintf('\nEikonal time delays completed in %0.3f seconds.\n', toc(tt));
            end

            % get the imaging grid
            gi = us.scan.getImagingGrid(); % {I1 x I2 x I3} each
            gi = gi(sel); % select and trim dimensions 

            % get sample times for each tx/rx (I1 x I2 x I3 x N x M)
            tau_rx = cellfun(@(f) f(gi{:}), rx_samp, 'UniformOutput', false); % all receive delays
            tau_tx = cellfun(@(f) f(gi{:}), tx_samp, 'UniformOutput', false); % all transmit delays
            tau_rx = cat(4, tau_rx{:}); % use all at a time
            tau_tx = cat(5, tau_tx{:}); % reconstruct matrix

            % short-circuit - empty matrix with output image sizing
            if kwargs.delay_only, b = zeros([us.scan.size, 1+kwargs.keep_rx*(us.rx.numel-1), 1+kwargs.keep_tx*(us.tx.numel-1), 0]); return; end

            % extract relevant arguments
            lut_args = ["apod", "fmod", "interp", "keep_tx", "keep_rx", "bsize"];
            args = namedargs2cell(rmfield(kwargs, setdiff(fieldnames(kwargs), lut_args)));

            % beamform
            b = bfDASLUT(us, chd, tau_rx, tau_tx, varargin{:}, args{:});
        end
    
        function [b, tau_rx, tau_tx] = bfDAS(us, chd, varargin, kwargs)
            % BFDAS - Delay-and-sum beamformer
            %
            % b = BFDAS(us, chd) creates a b-mode image b from the 
            % ChannelData chd.
            % 
            % b = BFDAS(..., Name, Value, ...) passes additional Name/Value
            % pair arguments
            % 
            % b = BFDAS(us, chd, A) uses an apodization 
            % defined by the ND-array A.
            %
            % b = BFDAS(..., 'apod', A) is supported for backwards
            % compatability.
            %
            % b = BFDAS(us, chd, A1, A2, ..., An) uses 
            % multiple separable apodization matrices A1, A2, ... An
            % defined by point-wise multiplication. Each matrix Ai must be
            % broadcastable to size I1 x I2 x I3 x N x M where I1 x I2 x I3
            % is the size of the image, N is the number of receivers, and M 
            % is the number of transmits. In other words, for each Ai,
            % ```
            % assert(all( ...
            %     size(Ai,1:5) == 1 | ...
            %     size(Ai,1:5) == [us.scan.size, us.xdc.numel, us.seq.numPulse] ...
            % ));
            %  ```
            % 
            % b = BFDAS(..., 'c0', c0) uses a beamforming sound speed 
            % of c0. c0 can be a scalar or an ND-array that is broadcastable
            % to size (I1 x I2 x I3) where [I1, I2, I3] == us.scan.size. 
            % The default is us.seq.c0.
            %             
            % b = BFDAS(..., 'keep_tx', true) preserves the transmit
            % dimension in the output image b.
            %
            % b = BFDAS(..., 'keep_rx', true) preserves the receive
            % dimension in the output image b.
            %
            % b = BFDAS(..., 'fmod', fc) upmixes the data at a modulation
            % frequency fc. This undoes the effect of
            % demodulation/downmixing at the same frequency.
            %
            % b = BFDAS(..., 'interp', method) specifies the method for
            % interpolation. Support is provided by the
            % ChannelData/sample2sep method. The default is 'cubic'.
            % 
            % b = BFDAS(..., 'bsize', B) uses a maximum block size of B
            % when vectorizing computations. A larger block size will run 
            % faster, but use more memory. The default is chosen by bfDASLUT.            
            %
            % [b, tau_rx, tau_tx] = BFDAS(...) additionally returns the
            % receive and transmit time delays tau_rx and tau_tx are size
            % (I1 x I2 x I3 x N x 1) and (I1 x I2 x I3 x 1 x M), where 
            % I1 x I2 x I3 is the size of the image, N is the number of
            % receivers, and M is the number of transmits.
            % 
            % [...] = BFDAS(..., 'delay_only',true) computes delays but
            % avoids computing the image.
            %        
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
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
            % See also DAS BFDASLUT BFADJOINT CHANNELDATA/SAMPLE PARCLUSTER PARPOOL

            arguments
                us (1,1) UltrasoundSystem
                chd ChannelData = ChannelData.empty
            end
            arguments(Repeating)
                varargin (:,:,:,:,:) {mustBeNumericOrLogical} % apodization matrices (I1 x I2 x I3 x N x M)
            end
            arguments
                kwargs.c0 (:,:,:,1,1) {mustBeNumeric} = us.seq.c0
                kwargs.apod {mustBeNumericOrLogical} = 1
                kwargs.fmod (1,1) {mustBeNumeric} = 0 
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "lanczos3"])} = 'cubic'
                kwargs.keep_tx (1,1) logical = false
                kwargs.keep_rx (1,1) logical = false
                kwargs.delay_only (1,1) logical = isempty(chd); % compute only delays
                kwargs.bsize (1,1) {mustBePositive, mustBeInteger}
            end

            % get image pixels, outside of range of data
            Pi = us.scan.positions; % 3 x I1 x I2 x I3

            % get the transmit, receive
            Pr = us.rx.positions(); % receiver positions

            % get virtual source or plane wave geometries
            switch us.seq.type
                case 'FSA',             [Pv, Nv] = deal(us.tx.positions(), argn(3, @()us.tx.orientations));
                case {'VS','FC','DV'},  [Pv, Nv] = deal(us.seq.focus, normalize(us.seq.focus - us.tx.offset,1,"norm"));
                case 'PW',              [Pv, Nv] = deal([0;0;0], us.seq.focus); % TODO: use origin property in tx sequence
            end
            Pr = swapdim(Pr,2,5); % move N to dim 5
            [Pv, Nv] = deal(swapdim(Pv,2,6), swapdim(Nv,2,6)); % move M to dim 6

            % get receive delays
            dr = vecnorm(Pi - Pr,2,1); % 1 x I1 x I2 x I3 x N x 1
                
            % transmit sensing vector
            dv = Pi - Pv; % 3 x I1 x I2 x I3 x 1 x M
            switch us.seq.type
                case {'DV','FSA'},  dv = vecnorm(dv, 2, 1)                         ;
                case {'VS','FC'},   dv = vecnorm(dv, 2, 1) .* sign(sum(dv .* Nv,1));
                case 'PW',          dv =                           sum(dv .* Nv,1) ;
            end % 1 x I1 x I2 x I3 x 1 x M

            % bring to I1 x I2 x I3 x 1 x M
            [dv, dr] = deal(reshape(dv,size(dv,2:6)), reshape(dr,size(dr,2:6)));

            % convert to time and alias
            dv = dv ./ kwargs.c0; 
            dr = dr ./ kwargs.c0;
            tau_rx = dr;
            tau_tx = dv;

            % short-circuit - empty matrix with output image sizing
            if kwargs.delay_only, b = zeros([us.scan.size, 1+kwargs.keep_rx*(us.rx.numel-1), 1+kwargs.keep_tx*(us.tx.numel-1), 0]); return; end

            % extract relevant arguments
            lut_args = ["apod", "fmod", "interp", "keep_tx", "keep_rx", "bsize"];
            args = namedargs2cell(rmfield(kwargs, setdiff(fieldnames(kwargs), lut_args)));

            % beamform
            b = bfDASLUT(us, chd, tau_rx, tau_tx, varargin{:}, args{:});
        end
    
        function b = bfDASLUT(us, chd, tau_rx, tau_tx, varargin, kwargs)
            % BFDASLUT - Look-up table delay-and-sum beamformer
            %
            % b = BFDASLUT(us, chd, tau_rx, tau_tx) creates a b-mode 
            % image b from the ChannelData chd using the time delay tables
            % tau_rx and tau_tx. 
            % 
            % The delay tables tau_rx and tau_tx must be broadcastable to
            % size [us.scan.size, chd.N, 1] and [us.scan.size, 1, chd.M]
            % respectively i.e. tau_rx and tau_tx should be in the order of
            % pixels x receives and pixels x transmits respectively.
            % 
            % b = BFDASLUT(us, chd, tau) uses the same delay table for
            % receive and transmit. This is valid primarily for FSA data
            % acquisitions.
            % 
            % b = BFDASLUT(us, chd, tau, A) or
            % b = BFDASLUT(us, chd, tau_rx, tau_tx, A) uses an apodization
            % defined by the ND-array A.
            %
            % b = BFDASLUT(..., 'apod', A) is supported for backwards
            % compatability.
            %
            % b = BFDASLUT(us, chd, tau, A1, A2, ..., An) or
            % b = BFDASLUT(us, chd, tau_rx, tau_tx, A1, A2, ..., An) uses 
            % multiple separable apodization matrices A1, A2, ... An
            % defined by point-wise multiplication. Each matrix Ai must be
            % broadcastable to size I1 x I2 x I3 x N x M where I1 x I2 x I3
            % is the size of the image, N is the number of receivers, and M 
            % is the number of transmits. In other words, for each Ai,
            % ```
            % assert(all( ...
            %     size(Ai,1:5) == 1 | ...
            %     size(Ai,1:5) == [us.scan.size, us.xdc.numel, us.seq.numPulse] ...
            % ));
            %  ```
            % 
            % b = BFDASLUT(..., Name, Value, ...) passes additional Name/Value
            % pair arguments
            % 
            % b = BFDASLUT(..., 'keep_tx', true) preserves the transmit
            % dimension in the output image b.
            %
            % b = BFDASLUT(..., 'keep_rx', true) preserves the receive
            % dimension in the output image b.
            %
            % b = BFDASLUT(..., 'fmod', fc) upmixes the data at a modulation
            % frequency fc. This undoes the effect of
            % demodulation/downmixing at the same frequency.
            %
            % b = BFDASLUT(..., 'interp', method) specifies the method for
            % interpolation. Support is provided by the
            % ChannelData/sample2sep method. The default is 'cubic'.
            %        
            % b = BFEIKONAL(..., 'bsize', B) uses a maximum block size of B
            % when vectorizing computations. A larger block size will run 
            % faster, but use more memory. The default is us.seq.numPulse 
            % or the largest value that restricts the apodization matrix to
            % 1GB.
            %
            % Example:
            % 
            % % Define the setup
            % us = UltrasoundSystem(); % get a default system
            % scat = Scatterers('pos', [0;0;30e-3], 'c0', us.seq.c0); % define a point target
            % 
            % % Compute the image
            % chd = greens(us, scat); % compute the response
            % 
            % % Get rx and tx delay tables for beamforming
            % [~, tau_rx, tau_tx] = bfDAS(us, chd, 'delay_only', true); % beamform the data
            % 
            % % Beamform
            % b = bfDASLUT(us, chd, tau_rx, tau_tx);
            % 
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(us.scan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also DAS BFDAS BFEIKONAL BFADJOINT CHANNELDATA/SAMPLE PARCLUSTER PARPOOL
            arguments
                us (1,1) UltrasoundSystem
                chd ChannelData
                tau_rx (:,:,:,:,1) {mustBeNumeric}
                tau_tx (:,:,:,1,:) {mustBeNumeric} = swapdim(tau_rx,4,5)
            end
            arguments(Repeating)
                varargin (:,:,:,:,:) {mustBeNumericOrLogical} % apodization matrices (I1 x I2 x I3 x N x M)
            end
            arguments
                kwargs.fmod (1,1) {mustBeNumeric} = 0
                kwargs.apod {mustBeNumericOrLogical} = 1
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "lanczos3"])} = 'cubic'
                kwargs.keep_tx (1,1) logical = false
                kwargs.keep_rx (1,1) logical = false
                kwargs.bsize (1,1) {mustBePositive, mustBeInteger} = max(1,min(max(size(tau_tx,5), us.seq.numPulse, 'omitnan'), floor(size(tau_tx,5)/(1*2^-30*8*prod(max(cell2mat(cellfun(@(x)size(x,1:5),varargin','uni',0)),[],1)))))) % heuristic 1GB  limit
            end
            
            % coerce apodization: ensure non-empty + backwards compatibility
            varargin = [{kwargs.apod}, varargin];
            
            % validate / parse receive table sizing
            N = unique([chd.N]);
            if ~isscalar(N)
                error( ...
                    "QUPS:UltrasoundSystem:bfDASLUT:nonUniqueReceiverSize", ...
                    "Expected a single receiver size, but instead they have sizes [" ...
                    + join(string(N),",") + "]." ...
                    )
            end
            if ~all(double([us.scan.size, N]) == size(tau_rx, 1:4))
                D = ndims(tau_rx); % last dim
                [I, L] = deal(prod(size(tau_rx, 1:D-1)), size(tau_rx, D)); % pixel and rx sizing
                if I == us.scan.nPix && L == N % sizing matches
                    tau_rx = reshape(tau_rx, [us.scan.size, N]);
                else
                    error( ...
                        "QUPS:UltrasoundSystem:bfDASLUT:incompatibleReceiveDelayTable", ...
                        "Expected a table with " + us.scan.nPix ...
                        + " pixels and " + N ...
                        + " receives but instead there are " + I ...
                        + " pixels and " + L + " receives." ...
                        );
                end
            end

            % validate / parse transmit table sizing
            M = unique([chd.M]);
            if ~isscalar(M)
                error( ...
                    "QUPS:UltrasoundSystem:bfDASLUT:nonUniqueTransmitSize", ...
                    "Expected a single transmit size, but instead they have sizes [" ...
                    + join(string(M),",") + "]." ...
                    )
            end
            if ~all(double([us.scan.size, 1, M]) == size(tau_tx, 1:5))
                D = ndims(tau_tx); % last dim
                [I, L] = deal(prod(size(tau_tx, 1:D-1)), size(tau_tx, D)); % pixel and rx sizing
                if I == us.scan.nPix && L == N % sizing matches
                    tau_tx = reshape(tau_tx, [us.scan.size, 1, M]);
                else
                    error( ...
                        "QUPS:UltrasoundSystem:bfDASLUT:incompatibleTransmitDelayTable", ...
                        "Expected a table with " + us.scan.nPix ...
                        + " pixels and " + M ...
                        + " transmits but instead there are " + I ...
                        + " pixels and " + L + " receives." ...
                        );
                end
            end

            % identify dimensions to sum after apodization
            sdim = [];
            if ~kwargs.keep_rx, sdim = [sdim, 4]; end
            if ~kwargs.keep_tx, sdim = [sdim, 5]; end

            % max dimension of data
            D = max([3, ndims(chd), cellfun(@ndims, {chd.data}), cellfun(@ndims, {chd.t0})]);

            % sample, apodize, and sum over tx/rx if requested
            for i = numel(chd):-1:1
                if isfinite(kwargs.bsize)
                    [chds, im] = splice(chd(i), chd.mdim, kwargs.bsize); % always splice tx
                    Ma = cellfun(@(x) size(x, 5), varargin); % size in tx dim
                    bi = {0}; % implicit preallocation
                    for m = numel(chds):-1:1 % each set of txs
                        tau_txm = sub(tau_tx, im(m), 5); % index delays per tx
                        a = sub(varargin{end}, unique(min(im{m},Ma(end))), 5); % init apodization per tx
                        for s = 1:numel(varargin)-1, a = a .* sub(varargin{s}, unique(min(im{m},Ma(s))), 5); end % reduce apodization per tx
                        bim = sample2sep(chds(m), tau_txm, tau_rx, kwargs.interp, a, sdim, kwargs.fmod, [4, 5]); % (I1 x I2 x I3 x [1|N] x [1|M] x [F x ... ])
                        if kwargs.keep_tx, bi{m} = bim;
                        else, bi{1} = bi{1} + bim;
                        end
                    end
                    bi = cat(5, bi{:});
                else  
                    a = varargin{end}; for s = 1:numel(varargin)-1, a = a .* varargin{s}; end % reduce apodization
                    bi = sample2sep(chd(i), tau_tx, tau_rx, kwargs.interp, a, sdim, kwargs.fmod, [4, 5]); % (I1 x I2 x I3 x [1|N] x [1|M] x [F x ... ])
                end

                % move aperture dimension to end
                bi = swapdim(bi, 4:5, D+(3:4)); % (I1 x I2 x I3 x 1 x 1 x [F x ... ] x [1|N] x [1|M])
                bi = reshape(bi, size(bi, [1:3, 6:ndims(bi)])); % (I1 x I2 x I3 x [F x ... ] x [1|N] x [1|M])
                b{i} = bi;
            end

            % combine to form output data
            chddim = find(size(chd) > 1, 1, 'first');
            if isempty(chddim), b = b{1}; % scalar
            else, b = cat(chddim, b{:}); end % array
        end

        function [b, bscan] = bfMigration(us, chd, Nfft, kwargs)
            % BFMIGRATION - Plane-Wave Stolt's f-k migration beamformer
            %
            % b = BFMIGRATION(us, chd) creates a b-mode image b from
            % the plane-wave ChannelData chd created from a TransducerArray
            % defined by us.xdc.
            %
            % [...] = BFMIGRATION(us, chd, [F, K]) uses a F-point FFT in time and
            % a K-point FFT laterally. If F < chd.T, the data is truncated temporally
            % and if K < chd.N, the data is truncated laterally. The default is
            % [chd.T, chd.N].
            %
            % [b, bscan] = BFMIGRATION(...) additionally returns a
            % ScanCartesian bscan on which the bmode image is naturally defined.
            %
            % [...] = BFMIGRATION(..., 'c0', c0) sets the sound speed. The
            % default is us.seq.c0.
            %
            % [...] = BFMIGRATION(..., 'keep_tx', true) preserves the
            % transmit dimension in the output image b.
            %
            % [...] = BFMIGRATION(..., 'interp', method) specifies the method for
            % interpolation. The default is 'cubic'.
            %
            % [...] = BFMIGRATION(..., 'bsize', B) uses an block size of B to
            % compute at most B transmits at a time. A larger block size
            % will run faster, but use more memory. The default is chosen
            % heuristically.
            %
            % [...] = BFMIGRATION(..., 'fmod', fc) upmixes the data at a
            % modulation frequency fc. This undoes the effect of
            % demodulation/downmixing at the same frequency.
            %
            % [...] = BFMIGRATION(..., 'jacobian', false) does not apply a
            % jacobian update when mapping the frequncies.
            %
            % References:
            % [1] Garcia D, Le Tarnec L, Muth S, Montagnon E, Porée J, Cloutier G.
            % Stolt's f-k migration for plane wave ultrasound imaging.
            % IEEE Trans Ultrason Ferroelectr Freq Control. 2013 Sep;60(9):1853-67.
            % doi: <a href="matlab:web('https://doi.org/10.1109%2FTUFFC.2013.2771')">10.1109/TUFFC.2013.2771</a>. PMID: 24626107; PMCID: PMC3970982.
            %
            % Example:
            % % A plane-wave transmission is required - migration works
            % % best for small angles
            % c0 = 1500;
            % seq = SequenceRadial('type', 'PW', 'angles', -10 : 0.25 : 10, 'c0', c0);
            %
            % % Define the setup
            % us = UltrasoundSystem('seq', seq); % get a default system
            % scat = Scatterers('pos', 1e-3*[0;0;1].*(5:5:30), 'c0', c0); % define a point target
            %
            % % Compute the response
            % chd = greens(us, scat);
            %
            % % beamform the data, with implicit zero-padding
            % [b, bscan] = bfMigration(us, chd, [2*chd.T, 4*chd.N]);
            %
            % % Display the image
            % bim = mod2db(b); % log-compression
            % figure;
            % imagesc(bscan, bim, [-80 0] + max(bim(:)));
            % colormap gray; colorbar;
            %
            % See also BFDAS BFADJOINT BFEIKONAL

            arguments
                us (1,1) UltrasoundSystem
                chd (1,1) ChannelData
                Nfft (1,2) {mustBeInteger, mustBePositive} = [chd.T, chd.N]; % FFT-lengths
                kwargs.c0 (1,1) {mustBeNumeric} = us.seq.c0
                kwargs.fmod (1,1) {mustBeNumeric} = 0 % modulation frequency
                kwargs.keep_tx (1,1) logical = false % whether to preserve transmit dimension
                kwargs.bsize (1,1) {mustBeNumeric, mustBeInteger, mustBePositive} = max(1,floor(2^(30-4-2-2) / prod([Nfft,prod(size(chd.data,4:max(4,ndims(chd.data))))]))); % vector computation block size - defaults to ~ 1GB
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "next", "previous", "spline", "pchip", "cubic", "makima", "freq", "lanczos3"])} = 'cubic'
                % kwargs.verbose (1,1) logical = true
                kwargs.jacobian (1,1) logical = true
            end

            % This is intended for plane waves and linear arrays
            if us.seq.type ~= "PW"
                warning('Expected a Sequence of type "PW", but instead it was type "' ...
                    + string(us.seq.type) ...
                    + '". Unexpected results may occur.');
            end
            if ~isa(us.xdc, 'TransducerArray')
                warning('Expected a TransducerArray but the Transducer is a ' ...
                    + string(class(us.xdc)) ...
                    + '". Unexpected results may occur.');
            end

            % Choose the FFT size in time/frequency and laterally
            [F, K] = deal(Nfft(1), Nfft(end));

            % Exploding Reflector Model velocity
            cs = kwargs.c0/sqrt(2);

            % get the frequency domains' axes with negative frequencies
            f  = ((0 : F - 1) - floor(F/2)) / F * chd.fs        ; % 1 x T - temporal frequencies
            kx = ((0 : K - 1) - floor(K/2)) / K / us.xdc.pitch; % 1 x K - lateral spatial frequencies

            % move to dimensions aligned to the data
            f  = swapdim(f , 2, chd.tdim);
            kx = swapdim(kx, 2, chd.ndim);

            % get array elements
            pn = us.xdc.positions;
            x0 = pn(1,1); % element lateral start position

            % get transmit mapping in compatible dimensions
            sq = copy(us.seq);
            sq.c0 = kwargs.c0;
            th0 = us.xdc.rot(1); % array orientation (azimuth)
            tau = sq.delays(us.xdc); % N x M
            ord = [chd.ndim, chd.mdim]; % send to these dimensions
            ord = [ord, setdiff(1:max(ord), ord)]; % account for all dimensions
            tau = ipermute(tau, ord); % permute to compatible dimensions
            gamma = swapdim(sind(sq.angles-th0) ./ (2 - cosd(sq.angles-th0)), 2, chd.mdim); % lateral scaling

    	    % splice the data to operate per block of transmits
    	    [chds, ix] = splice(chd, chd.mdim, kwargs.bsize); % split into groups of data
    	    tau   = arrayfun(@(i) {sub(tau  , i, chd.mdim)}, ix);
    	    gamma = arrayfun(@(i) {sub(gamma, i, chd.mdim)}, ix);

    	    chd0 = chd; % save og ChannelData, time delays
    	    if kwargs.keep_tx, bm = cell(1,numel(chds)); else, bm = 0; end % init
            for j = 1:numel(chds) % loop over blocks of transmits

                chd = chds(j); % reference the ChannelData

                % Move data to the temporal frequency domain
                x = chd.data;
                x = x .* exp(2j*pi*kwargs.fmod .* chd.time); % remodulate data
                x = fftshift(fft(x, F, chd.tdim), chd.tdim); % temporal fft

                % re-align time axis to the frequency domain
                x = x .* exp(-2j*pi*f .* chd.t0);

                % align transmits
                x = x .* exp(-2j*pi*f .* tau{j});

                % move to lateral frequency domain
                x = fftshift(fft(x, K, chd.ndim), chd.ndim); % lateral fft

                % get the Stolt's mapping from temporal frequency to spatial
                % (depth) frequency
                fkz = cs*sign(f).*sqrt(kx.^2 + f.^2 / cs^2);
                kkz = (fkz - f(1)) .* F ./ chd.fs; % convert to 0-based frequency index

                % resample using the spatio-temporal mapping
                y = wsinterpd(x, kkz, chd.tdim, 1, [], kwargs.interp, 0);

                % Jacobian (optional)
                if kwargs.jacobian
                    kz = f / cs;
                    y = (y .* kz) ./ (fkz + eps);
                end

                % re-align time axis to the time/space domain
                y = y .* exp(+2j*pi*f .* chd.t0);

                % move to the temporal domain
                b = ifft(ifftshift(y, chd.tdim), F, chd.tdim);
                tb = chd.t0 + (0 : F - 1)' ./ chd.fs;
                % tb = chd.time;

                % get the spatial axes
                zax = sq.c0 ./ 2 .* tb;
                xax = us.xdc.pitch .* (0 : K-1) + x0;

                % align data laterally using Garcia's PWI mapping
                b = b .* exp(2j*pi*kx.*gamma{j}.*zax);

                % move data back to spatial domain
                b = ifft(ifftshift(b, chd.ndim), K, chd.ndim);
                b = sub(b, {1:chd.T,1:chd.N}, [chd.tdim, chd.ndim]);

                % sum or store the transmits
                if kwargs.keep_tx, bm{j} = b; else, bm = bm + sum(b, chd.mdim); end

            end % for j

            % get full image cube
            if kwargs.keep_tx, b = cat(chd0.mdim, bm{:}); else, b = bm; end

            % create the corresponding scan - it aligns with our data
            bscan = ScanCartesian( ...
                'z', us.xdc.offset(3) + double(zax(1:chd0.T)), ...
                'x', xax(1:chd0.N), ...
                'y', us.xdc.offset(2) ...
                );

            % work-around: sometimes the numerical precision is
            % insufficient and interp2 is thrown off: ensure that the data
            % is regularly sampled by recreating the axes
            bscan.z = bscan.z(1) + mean(diff(bscan.z)) .* (0 : chd0.T - 1);

            % resample the data onto the original imaging grid if no
            % output scan was requested (risky)
            if nargout < 2
                warning("QUPS:bfMigration:artefacts", "Resampling a complex image can produce artefacts: request the output Scan to avoid resampling.");
                % resample data onto the given scan (risky)
                [sz, sx] = deal(us.scan.z, us.scan.x); % og scan
                [bz, bx] = ndgrid(bscan.z, bscan.x); % vectors -> matrix
                bint = num2cell(b, [1,2]); % place all transmits/frames in cells
                parfor(j = 1:numel(bint), 0) % interp for each transmit/frame
                    bint{j} = interp2(bx, bz, bint{j}, sx(:)', sz(:), kwargs.interp, 0); %#ok<PFBNS>
                end
                bint = cat(1, bint{:});
                bint = reshape(bint, [us.scan.size([1,2]), size(bint, 3 : max(3, ndims(bint)))]);
                b = bint;
            end
        end
    end

    % apodization functions
    methods
        function apod = apScanline(us, tol)
            % APSCANLINE - Create scanline apodization array
            %
            % apod = APSCANLINE(us) creates an ND-array apod to mask 
            % data using the UltrasoundSystem us in order to form an image
            % using scanlines.
            %
            % Scanline apodization is determined by accepting only
            % transmits and scanlines that are aligned across the transmit
            % aperture, where a scan line is a series of points where only
            % the range varies.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x M where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements, and M is the number of transmits.
            %
            % apod = APSCANLINE(us, tol) uses a tolerance of tol to decide 
            % whether to accept the transmit. The default is the mean
            % difference between foci.
            %
            % Example:
            %
            %   % Define the setup
            %   dx = 0.5e-3; % spacing between foci
            %   xf  = -10e-3 : dx : 10e-3; % focal lateral position
            %   pf  = [1;0;0].*xf + [0;0;1].*50e-3; % focal positions (50mm focal depth)
            %   seq = Sequence('type', 'FC', 'focus', pf, 'c0', 1500);
            %   us = UltrasoundSystem('seq', seq); % get a default system
            %   scat = Scatterers('pos', [0;0;20e-3], 'c0', us.seq.c0); % define a point target
            %
            %   % Compute the image
            %   chd = greens(us, scat); % compute the response
            %   chd = zeropad(singleT(chd), 0, max(0, chd.T - 2^9)); % precondition the data
            %   apod = apScanline(us, dx);
            %   b0 = DAS(us, chd, 'apod',    1); % beamform the data w/o apodization
            %   ba = DAS(us, chd, 'apod', apod); % beamform the data with apodization
            %
            %   % Display the images
            %   figure;
            %   bim0 = mod2db(b0); % log-compression
            %   nexttile(); imagesc(us.scan, bim0, [-80 0] + max(bim0(:)));
            %   colormap gray; colorbar; title("No apodization");
            %
            %   bima = mod2db(ba); % log-compression
            %   nexttile(); imagesc(us.scan, bima, [-80 0] + max(bima(:)));
            %   colormap gray; colorbar; title("Scanline apodization");
            %
            % See also: ULTRASOUNDSYSTEM/APMULTILINE ULTRASOUNDSYSTEM/APTRANSLATINGAPERTURE

            arguments
                us (1,1) UltrasoundSystem
                tol (1,1) {mustBePositive} = us.scan.("d"+scanlat(us.scan)); % numerical tolerance - e.g. us.scan.dx: assumes 2nd index the index of change
            end

            % soft validate the transmit sequence type: it should be focused
            styps = ["VS", "FC", "DV"]; % valid sequence types
            if all(us.seq.type ~= styps), warning(...
                    "Expected a Sequence of type " + join(styps, " or ") + ...
                    ", but instead got a Sequence of type " + us.seq.type + ": This may produce unexpected results." ...
                    );
            end

            % create a mask such that the transmit and pixel lateral
            % positions must 'match'
            if isa(us.scan, 'ScanCartesian')
                xi = swapdim(us.scan.x, 2, us.scan.xdim); % lateral per pixel
                xv = swapdim(sub(us.seq.focus,1,1), 2, 5); % 1 x 1 x 1 x 1 x M
            elseif isa(us.scan, 'ScanPolar')
                xi = swapdim(us.scan.a, 2, us.scan.adim); % angle per pixel
                xv = swapdim(us.seq.angles, 2, 5); % 1 x 1 x 1 x 1 x M
            else
                error("QUPS:UltrasoundSystem:UnsupportedScan", ...
                    "apMultiline does not support a " + class(us.scan) + " - please edit the code here." ...
                    );
            end
            apod = abs(xi - xv) < tol; % create mask
        end

        function apod = apMultiline(us)
            % APMULTILINE - Create a multi-line apodization array
            %
            % apod = APMULTILINE(us) creates an ND-array apod 
            % to mask delayed data using the UltrasoundSystem us in order
            % to form an image using scanlines.
            %
            % Multilne apodization is determined by linearly weighing
            % scan lines by the transmits that straddle each scan line.
            % Scan lines outside of the minimum and maximum transmit are
            % weighted by zero.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x M where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements, and M is the number of transmits.
            %
            % Example:
            %
            %   % Define the setup
            %   dx = 0.5e-3; % spacing between foci
            %   xf  = -10e-3 : dx : 10e-3; % focal lateral position
            %   pf  = [1;0;0].*xf + [0;0;1].*50e-3; % focal positions (50mm focal depth)
            %   seq = Sequence('type', 'FC', 'focus', pf, 'c0', 1500);
            %   us = UltrasoundSystem('seq', seq); % get a default system
            %   scat = Scatterers('pos', [0;0;20e-3], 'c0', us.seq.c0); % define a point target
            %
            %   % Compute the image
            %   chd = greens(us, scat); % compute the response
            %   chd = zeropad(singleT(chd), 0, max(0, chd.T - 2^9)); % precondition the data
            %   if isreal(chd), chd = hilbert(chd); end % complex analog
            %   apod = apMultiline(us);
            %   b0 = DAS(us, chd, 'apod',    1); % beamform the data w/o apodization
            %   ba = DAS(us, chd, 'apod', apod); % beamform the data with apodization
            %
            %   % Display the images
            %   figure;
            %   bim0 = mod2db(b0); % log-compression
            %   nexttile(); imagesc(us.scan, bim0, [-80 0] + max(bim0(:)));
            %   colormap gray; colorbar; title("No apodization");
            %
            %   bima = mod2db(ba); % log-compression
            %   nexttile(); imagesc(us.scan, bima, [-80 0] + max(bima(:)));
            %   colormap gray; colorbar; title("Multiline apodization");
            %
            % See also ULTRASOUNDSYSTEM/APSCANLINE ULTRASOUNDSYSTEM/APTRANSLATINGAPERTURE

            arguments
                us (1,1) UltrasoundSystem
            end

            % soft validate the transmit sequence type: it should be focused
            styps = ["VS", "FC", "DV"]; % valid sequence types
            if all(us.seq.type ~= styps), warning(...
                    "Expected a Sequence of type " + join(styps, " or ") + ...
                    ", but instead got a Sequence of type " + us.seq.type + ": This may produce unexpected results." ...
                    );
            end

            % extract lateral or angle of transmit in order to compare
            if isa(us.scan, 'ScanCartesian')
                xdim = us.scan.xdim; % dimension of change
                xi = swapdim(us.scan.x, 2, us.scan.xdim); % lateral per pixel
                xv = swapdim(sub(us.seq.focus,1,1), 2, 5); % 1 x 1 x 1 x 1 x M
            elseif isa(us.scan, 'ScanPolar')
                xdim = us.scan.adim; % dimension of change
                xi = swapdim(us.scan.a, 2, us.scan.adim); % angle per pixel
                xv = swapdim(us.seq.angles, 2, 5); % 1 x 1 x 1 x 1 x M
            else
                error("QUPS:UltrasoundSystem:UnsupportedScan", ...
                    "apMultiline does not support a " + class(us.scan) + " - please edit the code here." ...
                    );
            end

            % TODO: switch this to accept multiple transmits instead of
            % just left/right transmit

            % get the apodization based on interpolation weights
            % between each transmit angle
            da = xi - xv; % difference in lateral % Rx x Tx
            lind = da >= 0; % tx left  of xi
            rind = da <= 0; % tx right of xi

            % choose right-most left tx and left-most right tx
            lind = cellfun(@(i) find(i,1,'last' ), num2cell(lind,5), 'UniformOutput', false); % X x 1
            rind = cellfun(@(i) find(i,1,'first'), num2cell(rind,5), 'UniformOutput', false); % X x 1
            val = ~cellfun(@isempty, lind) & ~cellfun(@isempty, rind); % X x 1
            [lind, rind] = deal(cell2mat(lind(val)), cell2mat(rind(val))); % X' x 1

            % set the apodization values for these transmits
            da_lr = swapdim(abs(xv(lind) - xv(rind)), xdim,5); % difference in transmit angle (0 for identical left/right)
            a_l = 1 - (abs(swapdim(xv(lind),xdim,5) - xi(val)) ./ da_lr); % left side apodization
            a_r = 1 - (abs(swapdim(xv(rind),xdim,5) - xi(val)) ./ da_lr); % right side apodization
            ind0 = da_lr == 0; % if no difference between left/right transmits ...
            [a_l(ind0), a_r(ind0)] = deal(1, 0); %  set left to 1, right to 0

            % build Tx x Rx apodization matrix
            apod  = zeros([us.scan.size(xdim), us.seq.numPulse]); % build in reduced dimensions
            alind = sub2ind(size(apod), find(val), lind); % left matrix indices
            arind = sub2ind(size(apod), find(val), rind); % right matrix indices
            apod(alind) = apod(alind) + a_l; % add left apod
            apod(arind) = apod(arind) + a_r; % add right apod
            apod = ipermute(apod, [xdim, 5, setdiff(1:5, [xdim,5])]); % send X, M to dimensions xdim, 5
        end

        function apod = apTranslatingAperture(us, tol)
            % APTRANSLATINGAPERTURE - Create translating aperture apodization array
            %
            % apod = APTRANSLATINGAPERTURE(us) creates an ND-array apod to
            % mask data from the UltrasoundSystem us in order to emulate a 
            % translating transmit aperture beamformer configuration. The
            % transmit Sequence us.seq must be a focused.
            % 
            % If us.scan is a ScanCartesian, us.rx must be a
            % TransducerArray and the aperture is limited to receive
            % elements that are within tol of the focus laterally.
            %
            % If us.scan is a ScanPolar, us.rx must be a TransducerConvex
            % and us.seq must be a SequenceRadial. The aperture is
            % limited to receive elements that are within tol of the focus
            % in azimuth.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x M where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements, and M is the number of transmits.
            %
            % Example:
            %
            %   % Define the setup
            %   dx = 0.5e-3; % spacing between foci
            %   xf  = -10e-3 : dx : 10e-3; % focal lateral position
            %   pf  = [1;0;0].*xf + [0;0;1].*50e-3; % focal positions (50mm focal depth)
            %   seq = Sequence('type', 'FC', 'focus', pf, 'c0', 1500);
            %   us = UltrasoundSystem('seq', seq); % get a default system
            %   scat = Scatterers('pos', [0;0;20e-3], 'c0', us.seq.c0); % define a point target
            %
            %   % Compute the image
            %   chd = greens(us, scat); % compute the response
            %   chd = zeropad(singleT(chd), 0, max(0, chd.T - 2^9)); % precondition the data
            %   if isreal(chd), chd = hilbert(chd); end % complex analog
            %   apod = apTranslatingAperture(us, 64*us.xdc.pitch); % 64-element window
            %   b0 = DAS(us, chd, 'apod',    1); % beamform the data w/o apodization
            %   ba = DAS(us, chd, 'apod', apod); % beamform the data with apodization
            %
            %   % Display the images
            %   figure;
            %   bim0 = mod2db(b0); % log-compression
            %   nexttile(); imagesc(us.scan, bim0, [-80 0] + max(bim0(:)));
            %   colormap gray; colorbar; title("No apodization");
            %
            %   bima = mod2db(ba); % log-compression
            %   nexttile(); imagesc(us.scan, bima, [-80 0] + max(bima(:)));
            %   colormap gray; colorbar; title("Translating aperture apodization");
            %
            % See also ULTRASOUNDSYSTEM/APSCANLINE ULTRASOUNDSYSTEM/APMULTILINE

            arguments
                us (1,1) UltrasoundSystem
                tol (1,2) {mustBePositive} = us.scan.("d"+scanlat(us.scan)); % numerical tolerance - e.g. us.scan.dx
            end
            
            % soft validate the transmit sequence type: it should be focused
            styps = ["VS", "FC", "DV"]; % valid sequence types
            if all(us.seq.type ~= styps), warning(...
                    "Expected a Sequence of type " + join(styps, " or ") + ...
                    ", but instead got a Sequence of type " + us.seq.type + ": This may produce unexpected results." ...
                    );
            end

            % extract lateral or angle of transmit in order to compare
            if isa(us.scan, 'ScanCartesian')
                % soft error on transducer type
                if ~isa(us.rx, 'TransducerArray'), warning( ...
                        "Expected a TransducerArray but instead got " + class(us.rx) + ": This may produce unexpected results."...
                        ); end %#ok<ALIGN>
                xi = swapdim(us.scan.x, 2, us.scan.xdim); % lateral per pixel
                xv = swapdim(sub(us.seq.focus,1,1), 2, 5); % lateral per transmit 1 x 1 x 1 x 1 x M
                xn = swapdim(sub(us.rx.positions,1,1), 2,4); % lateral per receiver
            elseif isa(us.scan, 'ScanPolar')
                % soft error on transducer type
                if ~isa(us.rx, 'TransducerConvex'), warning( ...
                        "Expected a TransducerConvex but instead got " + class(us.rx) + ": This may produce unexpected results."...
                        ); end %#ok<ALIGN>
                xi = swapdim(us.scan.a, 2, us.scan.adim); % angle per pixel
                xv = swapdim(us.seq.angles, 2, 5); % angle per transmit (1 x 1 x 1 x 1 x M)
                xn = swapdim(us.rx.orientations,2,4); % angle per receiver (1 x 1 x 1 x N x 1)
            else
                error( ...
                    "QUPS:UltrasoundSystem:UnsupportedScan", ...
                    "UltrasoundSystem.apTranslatingAperture does not support a " + class(us.scan) + ". " ...
                    +newline+"Use a ScanCartesian or ScanPolar instead." ...
                    );
            end
            apod = abs(xi - xv) <= tol(1) & abs(xi - xn) <= tol(end); % create mask
        end

        function apod = apApertureGrowth(us, f, Dmax)
            % APAPERTUREGROWTH - Create an aperture growth aperture apodization array
            %
            % apod = APAPERTUREGROWTH(us) creates an
            % ND-array to mask delayed data using the transmit Sequence seq
            % and the receive Transducer rx in order to form the
            % corresponding image that shrinks the receive aperture at
            % shallower depths in order to maintain a minimum f-number.
            %
            % apod = APAPERTUREGROWTH(us, f) uses an
            % f-number of f to limit the aperture. The default is 1.5.
            %
            % apod = APAPERTUREGROWTH(us, f, Dmax) restricts the
            % maximum size of the aperture to Dmax. The default is Inf.
            %
            % The Transducer us.rx must be a TransducerArray or 
            % TransducerConvex.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x M where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements, and M is the number of transmits.
            %
            %
            % Example:
            %
            %   % Define the setup
            %   dx = 0.5e-3; % spacing between foci
            %   xf  = -10e-3 : dx : 10e-3; % focal lateral position
            %   pf  = [1;0;0].*xf + [0;0;1].*40e-3; % focal positions (50mm focal depth)
            %   seq = Sequence('type', 'FC', 'focus', pf, 'c0', 1500);
            %   us = UltrasoundSystem('seq', seq); % get a default system
            %   scat = Scatterers('pos', [0;0;20e-3], 'c0', us.seq.c0); % define a point target
            %
            %   % Compute the image
            %   chd = greens(us, scat); % compute the response
            %   chd = zeropad(singleT(chd), 0, max(0, chd.T - 2^9)); % precondition the data
            %   if isreal(chd), chd = hilbert(chd); end % complex analog
            %   apod = apApertureGrowth(us, 2); % use a f# of 2
            %   b0 = DAS(us, chd, 'apod',    1); % beamform the data w/o apodization
            %   ba = DAS(us, chd, 'apod', apod); % beamform the data with apodization
            %
            %   % Display the images
            %   figure;
            %   bim0 = mod2db(b0); % log-compression
            %   nexttile(); imagesc(us.scan, bim0, [-80 0] + max(bim0(:)));
            %   colormap gray; colorbar; title("No apodization");
            %
            %   bima = mod2db(ba); % log-compression
            %   nexttile(); imagesc(us.scan, bima, [-80 0] + max(bima(:)));
            %   colormap gray; colorbar; title("Aperture growth apodization");
            %
            %
            % See also APACCEPTANCEANGLE APCOSINEANGLE

            % defaults
            arguments
                us (1,1) UltrasoundSystem
                f (1,1) {mustBePositive} = 1.5
                Dmax (1,1) {mustBePositive} = Inf
            end

            % soft validate the transmit sequence and transducer types.
            if ~(isa(us.rx, 'TransducerArray') || isa(us.rx, 'TransducerConvex')), warning(...
                    "Expected Transducer to be a TransducerArray or TransducerConvex but instead got a " + class(us.rx) + ". This may produce unexpected results."...
                    )
            end

            % get the transmit foci lateral position, M in dim 6
            % Pv = swapdim(sub(seq.focus,1:3,1), 2, 6);

            % get the receiver lateral/axial positions, N in dim 5
            Pn = swapdim(us.rx.positions, 2, 5); % (3 x 1 x 1 x 1 x N)
            [Xn, Zn] = deal(sub(Pn,1,1), sub(Pn,3,1)); % (1 x 1 x 1 x 1 x N)

            % calculate x and z for scan
            if isa(us.scan, 'ScanCartesian') % already in Cartesian
                % get the pixel positions in the proper dimensions
                Xi = swapdim(us.scan.x, 2, 1+us.scan.xdim); % (1 x 1 x I2 x 1)
                Zi = swapdim(us.scan.z, 2, 1+us.scan.zdim); % (1 x I1 x 1 x 1)
            else % generic -> convert to Cartesian
                [Xi, ~, Zi] = us.scan.getImagingGrid(); % (I1 x I2 x 1)
                [Xi, Zi] = dealfun(@(x)reshape(x,[1 size(x)]), Xi, Zi); % (1 x I1 x I2 x 1)
            end 

            % get the equivalent aperture width (one-sided) and pixel depth
            if any(us.rx.orientations) % non-planar array width and depth calculation
                ae = swapdim(us.rx.orientations(), 2, 5); % angles of elements, in degrees
                rp = hypot( Xi - Xn, Zi - Zn); % radii to scan points
                ap = atan2d(Xi - Xn, Zi - Zn); % angles to scan points
                d =     rp .* sind(ap - ae) ; % equivalent one-sided widths for points
                z = abs(rp .* cosd(ap - ae)); % equivalent depth for points; take abs() to use same apodization calculation as linear arrays
            else % planar array width and depth calculation
                d = Xn - Xi; % one-sided width where 0 is aligned with the pixel
                z = Zi - 0; % depth w.r.t. beam origin (for linear xdcs)
            end

            % determine apodization 
            apod = z > f *  abs(2*d) ; % restrict the f-number
            apod = apod .* (abs(2*d) < Dmax); % also restrict the maximum aperture size

            % shift to (I1 x I2 x I3 x N x M)
            apod = reshape(apod, size(apod,2:6));
        end
        
        function apod = apTxParallelogram(us, theta, phi)
            % apTxParallelogram - Plane-wave illumination apodization
            %
            % apod = apTxParallelogram(us, theta) creates a pixel by
            % transmit based apodization mask apod for pixels that lie
            % within the lateral boundaries of the transducer when
            % projected along the transmit angles theta.
            %
            % apod = apTxParallelogram(us, theta, [phi_min, phi_max])
            % expands the accepted angle of projection theta to 
            % theta + [phi_min, phi_max]. The default is [0 0].
            %
            % Example:
            % seq = SequenceRadial('angles', -10 : 5 : 10); % PW sequence
            % us = UltrasoundSystem('seq', seq);
            % ap = us.apTxParallelogram(seq.angles, [-5 5]);
            % figure; imagesc(us.scan, ap);
            % animate(ap,'fs',1,'loop',false,"title","Angle: "+seq.angles);
            % 
            % See also APACCEPTANCEANGLE APCOSINEANGLE
            arguments
                us (1,1) UltrasoundSystem
                theta (1,:) = atan2d(us.seq.focus(1,:), us.seq.focus(3,:));
                phi (1,2) double = 0;
            end
            phi = swapdim(phi,2,7); % min max  difference in angle
            pg = us.scan.positions();
            nv = swapdim([sind(phi+theta); 0*(phi+theta); cosd(phi+theta)],2,6);  % normal vectors (Tx in dim 6)
            pg = pg - nv .* (sub(pg,3) ./ sub(nv,3)); % project to z == 0
            pb = us.xdc.bounds();
            apod = any(pb(1,1) < sub(pg,1,1),7) & any(sub(pg,1,1) <= pb(1,2),7); % in-bounds of xdc
            apod = reshape(apod, size(apod,2:6)); % to I x N x M
        end

        function apod = apAcceptanceAngle(us, theta)
            % APACCEPTANCEANGLE - Create an acceptance angle apodization array
            %
            % apod = APACCEPTANCEANGLE(us) creates an ND-array to
            % mask delayed data from the UltrasoundSystem us which includes
            % receive elements where the pixel to element angle is within
            % an acceptance angle.
            %
            % apod = APACCEPTANCEANGLE(us, theta) uses an acceptance angle
            % of theta in degrees. The default is 45.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x 1 where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements.
            %
            % Example:
            %
            %   % Define the setup
            %   dx = 0.5e-3; % spacing between foci
            %   xf  = -10e-3 : dx : 10e-3; % focal lateral position
            %   pf  = [1;0;0].*xf + [0;0;1].*40e-3; % focal positions (50mm focal depth)
            %   seq = Sequence('type', 'FC', 'focus', pf, 'c0', 1500);
            %   us = UltrasoundSystem('seq', seq); % get a default system
            %   scat = Scatterers('pos', [0;0;20e-3], 'c0', us.seq.c0); % define a point target
            %
            %   % Compute the image
            %   chd = greens(us, scat); % compute the response
            %   chd = zeropad(singleT(chd), 0, max(0, chd.T - 2^9)); % precondition the data
            %   if isreal(chd), chd = hilbert(chd); end % complex analog
            %   apod = apAcceptanceAngle(us, 25); % use a limit of 25   degrees
            %   b0 = DAS(us, chd, 'apod',    1); % beamform the data w/o apodization
            %   ba = DAS(us, chd, 'apod', apod); % beamform the data with apodization
            %
            %   % Display the images
            %   figure;
            %   bim0 = mod2db(b0); % log-compression
            %   nexttile(); imagesc(us.scan, bim0, [-80 0] + max(bim0(:)));
            %   colormap gray; colorbar; title("No apodization");
            %
            %   bima = mod2db(ba); % log-compression
            %   nexttile(); imagesc(us.scan, bima, [-80 0] + max(bima(:)));
            %   colormap gray; colorbar; title("Acceptance angle apodization");
            %
            % See also APAPERTUREGROWTH APCOSINEANGLE

            % defaults
            arguments
                us (1,1) UltrasoundSystem
                theta (1,1) {mustBePositive} = 45
            end

            % get the receiver positions and orientations, N in dim 5
            Pn      = swapdim(us.rx.positions, 2, 5); % (3 x 1 x 1 x 1 x N) % positions
            [~,~,n] = us.rx.orientations;
            n = swapdim(n, 2, 5); % (3 x 1 x 1 x 1 x N) % element normals

            % get the image pixels
            Pi = us.scan.positions(); % (3 x I1 x I2 x I3 x 1)

            % get the normalized difference vector (3 x I1 x I2 x I3 x N)
            r = Pi - Pn;

            % get the cosine of the angle via the inner product of the
            % normalized direction vectors 
            % (1 x I1 x I2 x I3 x N)
            r = r ./ vecnorm(r,2,1); % normalize
            r = pagemtimes(n, 'transpose', r, 'none'); 
            r = reshape(r, size(r,2:6)); % -> (I1 x I2 x I3 x N x 1)

            % accept if greater than the cutoff angle
            apod = r >= cosd(theta);
        end    


        function apod = apCosineAngle(us, theta, kwargs)
            % APCOSINEANGLE - Create an cosine-weighted apodization array
            %
            % apod = APCOSINEANGLE(us) creates an ND-array to weight
            % delayed data from the UltrasoundSystem us by the angle `phi`
            % between the element normal and the element to pixel vector.
            % The weighting is
            % 
            % `cosd(min(90, (phi / theta)))`
            % 
            % where `theta` is the maximum accepted angle.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x 1 where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements.
            %
            % apod = APCOSINEANGLE(us, theta) uses an acceptance angle
            % of theta in degrees. The default is 45.
            %
            % apod = APCOSINEANGLE(..., 'gpu', false) avoids using a GPU.
            % The default is true if a GPU is available and the
            % intermediate memory usage will not exceed 4GB.
            % 
            % Example:
            % us = UltrasoundSystem(); % default
            % a = us.apCosineAngle(single(45)); % apodization ND-array
            % 
            % figure; imagesc(us.scan, a); title("Cosine Apodization"); % display
            % hold on; plot(us.xdc, 'b','Marker', 'square'); legend(); % add Transducer
            % colorbar; colormap hot; clim([0 1]); % format
            % animate(a, 'fn', true, 'loop', false); % animate vs. elements
            % 
            % See also APAPERTUREGROWTH APACCEPTANCEANGLE

            % defaults
            arguments
                us (1,1) UltrasoundSystem
                theta (1,1) {mustBePositive} = 45
                kwargs.gpu (1,1) logical = canUseGPU() && (4 * us.scan.nPix * us.xdc.numel < 4 * 2^30 / 2^3)
            end

            % cosine apodization
            [pg, pn] = deal(us.scan.positions(), us.xdc.positions()); % pixels | elems
            [~,~,nn] = us.xdc.orientations(); % elem normals
            [pn, nn] = deal(swapdim(pn,2,5), swapdim(nn,2,5)); % match dims
            if kwargs.gpu, theta = gpuArray(theta); end % imlpicit casting
            [pn, nn, pg] = dealfun(@(x) cast(x, 'like', theta), pn, nn, pg);
            r  = (pg - pn); % elem -> pixel vector
            r  = r ./ vecnorm(r,2,1); % normalized
            r  = swapdim(pagemtimes(nn, 'transpose', r, 'none'), 2:6, 1:5); % normalized inner product
            r = max(-1,min(1,r)); % coerce in-bounds (numerical precision issues)
            apod = cosd(min(90,(90/theta)*acosd(r))); % cosine gradient with (scaled) angle
        end
    end

    % dependent methods
    methods
        function f = get.fc(us)
            if us.rx == us.tx || isalmostn(us.rx.fc, us.tx.fc)
                f = us.rx.fc;
            else
                f = [us.tx.fc, us.rx.fc];
            end
        end
        function l = get.lambda(us), l = us.seq.c0 ./ us.fc; end
        function x = get.xdc(us)
            if us.tx == us.rx
                x = us.rx;
            else
                error( ...
                    "QUPS:UltrasoundSystem:ambiguousTransducer", ...
                    "Call to 'xdc' ambiguous; transmitter and receiver are not the same." ...
                    );
            end
        end
        function set.xdc(us, xdc), [us.tx, us.rx] = deal(xdc); end
        function s = get.sequence(us), warning("QUPS:UltrasoundSystem:syntaxDeprecated","UltrasoundSystem.sequence is deprecated. Use the .seq property instead."); s = us.seq; end
        function set.sequence(us, s), warning("QUPS:UltrasoundSystem:syntaxDeprecated","UltrasoundSystem.sequence is deprecated. Use the .seq property instead."); us.seq = s; end
    end
    
    % recompilation helper functions
    methods
        function [mcom, nvcom] = recompile(us), mcom = recompileMex(us); if gpuDeviceCount, nvcom = recompileCUDA(us); else, nvcom = string.empty; end, end
        % RECOMPILE - Recompile mex and CUDA files
        %
        % RECOMPILE(us) recompiles all mex binaries and CUDA files and 
        % stores them in us.tmp_folder. If there are no MATLAB compatible
        % GPUs, it does not attempt to recompile CUDA files.
        %
        % Example:
        % us = UltrasoundSystem();
        % us.recompile();
        % ls(us.tmp_folder)
        % 
        % See also ULTRASOUNDSYSTEM.RECOMPILECUDA ULTRASOUNDSYSTEM.RECOMPILEMEX
        function mcom = recompileMex(us, defs)
            % RECOMPILEMEX - Recompile mex files
            %
            % RECOMPILEMEX(us) recompiles all mex binaries and stores
            % them in us.tmp_folder.
            %
            % RECOMPILEMEX(us, defs) compiles for the compiler 
            % definition structs defs. These structures are generated by
            % the UltrasoundSystem class. The default is all the
            % definitions returned by UltrasoundSystem.genMexdefs().
            % 
            % mcom = RECOMPILEMEX(...) returns the commands sent to the mex
            % compiler. This command can be tweaked and resent to debug or 
            % fix compilation errors.
            %
            % Example:
            % % Create a new system and recompile 
            % us = UltrasoundSystem(); % create a new system
            % mcom = us.recompileMex(); % recompile mex files
            %
            % % change the source code ...
            % 
            % % recompile the last mex file manually
            % args = cellstr(mcom(:,end)); % extract separated arguments
            % mex(args{:}); % send to mex to compile
            % 
            % See also ULTRASOUNDSYSTEM.GENMEXDEFS ULTRASOUNDSYSTEM.RECOMPILE

            arguments
                us (1,1) UltrasoundSystem
                defs (1,:) struct = UltrasoundSystem.genMexdefs();
            end

           % compile each definition
            for i = numel(defs):-1:1
                d = defs(i);
                % make full command
                mcom(:,i) = cat(1,...
                    '-outdir', us.tmp_folder, ... place binaries in system's tmp folder
                    join("-I" + d.IncludePath), ...
                    join("-L" + d.Libraries), ...
                    join("-D" + d.DefinedMacros),...
                    fullfile(UltrasoundSystem.getSrcFolder(), 'FMM', 'functions', d.Source) ...
                    );
                
                try
                    arg = cellstr(mcom(:,i));
                    s = mex(arg{:});
                    if s, warning("Unable to recompile " + d.Source); else, disp("Success recompiling " + d.Source); end
                catch
                    warning("Unable to recompile code!");
                end
            end
   
        end
        function com = recompileCUDA(us, defs, arch, kwargs)
            % RECOMPILECUDA - Recompile CUDA ptx files
            %
            % RECOMPILECUDA(us) recompiles all CUDA files to ptx and stores
            % them in us.tmp_folder.
            %
            % RECOMPILECUDA(us, defs) compiles for the compiler 
            % definition structs defs. These structures are generated by
            % the UltrasoundSystem class. The default is all the
            % definitions returned by UltrasoundSystem.genCUDAdefs().
            %
            % RECOMPILECUDA(us, defs, arch) controls the architecture
            % string used in compilation. These must start with 
            % 'compute_' e.g. 'compute_86' for RTX 3000 series or 
            % 'compute_75' for RTX 2000 series. This typically matches the 
            % "compute capability" of the device.
            % 
            % If arch is an empty string, no architecture argument is used.
            % The default is the architecture of the current GPU device.
            %
            % RECOMPILECUDA(..., 'mex', true) uses `mexcuda` rather than
            % `nvcc` to compile ptx files (requires R2023a or later). 
            % Architecture selection, warning selection and suppression,
            % and other compiler options are not available via mexcuda. The
            % default is true if `nvcc` is not found on the path.
            % 
            % nvcom = RECOMPILECUDA(...) returns the command sent to nvcc.
            % This command can be tweaked and resent to debug or fix
            % compilation errors.
            %
            % Example:
            % % Create a new system and recompile
            % us = UltrasoundSystem(); % create a new system
            % nvcom = us.recompileCUDA(); % recompile CUDA files
            %
            % % change the source code ...
            %
            % % recompile the 3rd CUDA file manually
            % system(nvcom(3)); % via nvcc
            % % mexcuda(nvcom{3}{:}); % via mexcuda
            % 
            % See also ULTRASOUNDSYSTEM.RECOMPILE ULTRASOUNDSYSTEM.RECOMPILEMEX 
            % ULTRASOUNDSYSTEM.GENCUDADEFS

            arguments
                us (1,1) UltrasoundSystem
                defs (1,:) struct = UltrasoundSystem.genCUDAdefs();
                arch (:,1) string {mustBeScalarOrEmpty, mustBeArch(arch)} = nvarch()
                kwargs.mex (1,1) logical = isempty(argn(2, @system, 'which nvcc')) && ~isMATLABReleaseOlderThan("R2023a")
            end

            % ensure that the output folder exists
            if ~exist(us.tmp_folder,'dir'), warning("QUPS:UltrasoundSystem:ReinitializingTmpFolder","Recreating a tmp_folder that was deleted (somehow)."); mkdir(us.tmp_folder); addpath(us.tmp_folder); end
            
            % compile each
            for i = numel(defs):-1:1
                d = defs(i);
                % make full command
                if ~kwargs.mex % via nvcc
                com(i) = join(cat(1,...
                    "nvcc ", ...
                    "--ptx " + which(string(d.Source)), ...
                    "-arch=" + arch + " ", ... compile for active gpu
                    "-o " + fullfile(us.tmp_folder, argn(2, @fileparts, d.Source) + ".ptx"), ...
                    join("--" + d.CompileOptions),...
                    join("-I" + d.IncludePath), ...
                    join("-L" + d.Libraries), ...
                    join("-W" + d.Warnings), ...
                    join("-D" + d.DefinedMacros)...
                    ));
                    comp = @system;
                else % via mexcuda
                com{i} = cellstr(cat(2,...
                    "-ptx", ...
                    "-output", fullfile(us.tmp_folder, argn(2, @fileparts, d.Source) + ".ptx"), ...
                    ... ("--" + d.CompileOptions),...
                    ("-I" + d.IncludePath), ...
                    ("-L" + d.Libraries), ...
                    ... ("-W" + d.Warnings), ...
                    ("-D" + d.DefinedMacros),...
                    which(string(d.Source)) ...
                    ... "-arch=" + arch + " ", ... compile for active gpu
                    )');
                    comp = @(s) mexcuda(s{1}{:});
                end                
                try s = comp(com(i));
                    if s, warning( ...
                            "QUPS:recompile:UnableToRecompile", ...
                            "Unable to recompile " + d.Source ...
                            + " (" + s + ")"); 
                    else
                        disp("Success recompiling " + d.Source); 
                    end
                catch
                    warning("QUPS:recompile:UnableToRecompile", "Unable to recompile code!");
                end
            end
        end
        function def = getDASConstCudaDef(us, chd, varargin, kwargs)
            % GETDASCONSTCUDADEF - Constant size compilation definition for DAS
            %
            % def = GETDASCONSTCUDADEF(us) creates a compilation
            % definition for the CUDA executables used by
            % UltrasoundSystem.DAS for the current Scan, Transducer, and
            % Sequence. 
            % 
            % Using a fixed data size triggers compiler-level optimizations
            % that can improve performance for iterative calls. This comes
            % at the cost of compile time and flexibility.
            %
            % After compiling with this definition, calls with a different
            % input or output data size will have unexpected results. Use
            % us.recompileCUDA() to reset the binaries to handle variable sizes. 
            % 
            % def = GETDASCONSTCUDADEF(us, chd) additionally uses a fixed
            % size ChannelData object.
            % 
            % def = GETDASCONSTCUDADEF(us, chd, a1, a2, ..., an)
            % additionally sets the number of apodization matrices. If not
            % included, the default number of apodization matrices is 0.
            % 
            % def = GETDASCONSTCUDADEF(..., 'keep_tx', true) preserves the transmit
            % dimension.
            %
            % def = GETDASCONSTCUDADEF(..., 'keep_rx', true) preserves the receive
            % dimension.
            % 
            % def = GETDASCONSTCUDADEF(..., 'interp', method) specifies the method for
            % interpolation. Must be one of {'nearest', 'linear', 'cubic', 'lanczos3'}.
            % The default is 'cubic'.
            %
            % def = GETDASCONSTCUDADEF(..., 'interp', 'none') avoids specifying the 
            % interpolation method. Use this to also avoid specifying the preservation
            % of the receive and transmit dimensions.
            %
            % def = GETDASCONSTCUDADEF(..., 'no_apod', true) avoids
            % setting the number of apodization matrices. The default is
            % false.
            % 
            % Example:
            % % Setup the data
            % us = UltrasoundSystem();
            % sct = Scatterers();
            % chd = greens(us, sct);
            % apod = us.apAcceptanceAngle(30);
            % 
            % % Compile and extract
            % def = getDASConstCudaDef(us, chd, apod, 'keep_rx', true);
            % us.recompileCUDA(def);
            % [b, k, PRE_ARGS, POST_ARGS] = DAS(us, chd, apod, 'keep_rx', true);
            % 
            % % Compute
            % F = size(chd.data(:,:,:,:),4); % frames
            % for f = 1:F
            %     bf{f} = k.feval(PRE_ARGS{:}, chd.data(:,:,:,f), POST_ARGS{:});
            % end
            % bf = cat(6, bf{:}); % concatenate frames in the 6th dimension
            % 
            % See also ULTRASOUNDSYSTEM.DAS ULTRASOUNDSYSTEM.GENCUDADEFS 
            % ULTRASOUNDSYSTEM.RECOMPILE

            arguments
                us (1,1) UltrasoundSystem
                chd {mustBeScalarOrEmpty} = ChannelData.empty
            end
            arguments(Repeating)
                varargin % apodization matrices
            end
            arguments
                kwargs.no_apod (1,1) = false % don't set apodization sizes
                kwargs.interp (1,1) string {mustBeMember(kwargs.interp, ["linear", "nearest", "cubic", "lanczos3", "none"])} = 'cubic'
                kwargs.keep_tx (1,1) logical = false % whether to preserve transmit dimension
                kwargs.keep_rx (1,1) logical = false % whether to preserve receive dimension
            end
            
            % get the other sizes for beamform.m
            VS = ~(us.seq.type == "PW"); % whether virtual source or plane-wave
            DV = any(us.seq.type == ["DV", "FSA"]); % whether virtual source or plane-wave
            Isz = us.scan.size; % size of the scan
            N = us.rx.numel; % number of receiver elements
            M = us.seq.numPulse; % number of transmits
            assert(~isnan(M) && M > 0); % M must be defined
            
            % get all source code definitions
            def = UltrasoundSystem.genCUDAdefs('beamform');

            % get options
            switch kwargs.interp
                case "nearest", flagnum = 0;
                case "linear",  flagnum = 1;
                case "cubic",   flagnum = 2;
                case "lanczos3",flagnum = 3;
                case "none",    flagnum = nan;
                otherwise, error('Interp option not recognized: ' + string(interp));
            end
            flagnum = flagnum + 8*kwargs.keep_rx + 16*kwargs.keep_tx + 32*(chd.ndim>chd.mdim);
            
            % add the defined macros
            def.DefinedMacros = cat(1, ...
                def.DefinedMacros, ... keep the current defs
                "QUPS_" + {... prepend 'QUPS_'
                "VS="+VS,... virtual source model
                "DV="+DV,... diverging wave model
                "N="+N,... elements
                "M="+M,... transmits
                "I1="+Isz(1),... pixel dim 1
                "I2="+Isz(2),... pixel dim 2
                "I3="+Isz(3) ... pixel dim 3
                }');
            if  kwargs.interp ~= "none"
                def.DefinedMacros(end+1) = "QUPS_BF_FLAG="+flagnum; % interp / integration flag
            end
            if ~kwargs.no_apod
                def.DefinedMacros(end+1) = "QUPS_S="+numel(varargin); % # of apod matrices
            end
            
            % if T is provided, include it
            if isscalar(chd)
                def.DefinedMacros = [def.DefinedMacros; ...
                    {"QUPS_T="+chd.T}; ... number of time samples
                    ];
            end
        end
        function def = getGreensConstCudaDef(us, scat)
            % GETGREENSCONSTCUDADEF - Constant size compilation definition for greens
            %
            % def = GETGREENSCONSTCUDADEF(us) creates a compilation
            % definition for the CUDA executables used by
            % UltrasoundSystem.greens for the current Transducer, Sequence,
            % and sampling frequency.
            % 
            % Using a fixed data size triggers compiler-level optimizations
            % that can improve performance for iterative calls. This comes
            % at the cost of compile time and flexibility.
            %
            % After compiling with this definition, calls with a different
            % size of data will have unexpected results. Use
            % UltrasoundSystem.recompileCUDA to reset the binaries to
            % handle variable sizes. 
            % 
            % def = GETGREENSCONSTCUDADEF(us, scat) additionally uses a
            % fixed data size from the Scatterers scat.
            % 
            % See also ULTRASOUNDSYSTEM.GREENS ULTRASOUNDSYSTEM.GENCUDADEFS 
            % ULTRASOUNDSYSTEM.RECOMPILE

            arguments
                us (1,1) UltrasoundSystem
                scat Scatterers {mustBeScalarOrEmpty} = Scatterers.empty
            end

            % get the Waveform length
            wv = conv(us.rx.impulse, ...
                conv(us.tx.impulse, us.seq.pulse, us.fs), ...
                us.fs); % transmit waveform, convolved at US frequency
            wv.fs = us.fs;
            T = length(wv.samples);

            % get the other sizes for greens.cu
            N = us.rx.numel; % number of receiver elements
            M = us.tx.numel; % number of transmits

            % get all source code definitions
            def = UltrasoundSystem.genCUDAdefs('greens');

            % add the defined macros
            def.DefinedMacros = cat(1, ...
                def.DefinedMacros, ... keep the current defs
                "QUPS_" + {... prepend 'QUPS_'
                "T="+T,... elements
                "N="+N,... elements
                "M="+M,... transmits
                }');

            % if I is provided, include it
            if isscalar(scat)
                def.DefinedMacros = [def.DefinedMacros; ...
                    {"QUPS_I="+scat.numScat}; ... number of time samples
                    ];
            end
        end
    end
   
    % source file recompilation definitions
    methods(Static,Hidden)
        function f = getSrcFolder(), f = fileparts(mfilename('fullpath')); end
    end
    methods(Static)
        function defs = genCUDAdefs(name)
            % GENCUDADEFS - Generate CUDA compilation definitions.
            % 
            % defs = GENCUDADEFS() returns a struct array of definition
            % structs defs to compile all CUDA kernels used by QUPS.
            %
            % defs = GENCUDADEFS(name) returns the kernels specified within
            % the string array name.
            % 
            % defs = GENCUDADEFS("all") returns all available kernels. The 
            % default is "all".
            % 
            % Use with UltrasoundSystem.recompileCUDA to compile the
            % kernels for the UltrasoundSystem.
            % 
            % Example:
            % us = UltrasoundSystem('recompile', false);  
            % defs = UltrasoundSystem.genCUDAdefs(["interpd", "convd"]);
            % us.recompileCUDA(defs);
            % 
            % See also ULTRASOUNDSYSTEM.RECOMPILECUDA ULTRASOUNDSYSTEM.GETDASCONSTCUDADEF 
            % ULTRASOUNDSYSTEM.GETGREENSCONSTCUDADEF 
            arguments
                name (1,:) string {mustBeMember(name, ["all" , ...
                    "beamform", "interpd", "wbilerp", "greens", "convd" ...
                   ])} = "all"
            end
            if isscalar(name) && name == "all"
                name = ["beamform", "interpd", "wbilerp", "greens", "convd"];
            end
            for i = numel(name):-1:1
            switch name(i)
                case "beamform"
                    def.Source = {...
                        'bf.cu', ...
                        }';

                    def.IncludePath = {}; % include folders
                    def.Libraries = {}; % libraries

                    def.CompileOptions = {...  compile options
                        'use_fast_math',...
                        };

                    def.Warnings = {... warnings
                        'no-deprecated-gpu-targets'...
                        };

                    def.DefinedMacros = {};
                case "interpd"
                    def.Source = {...
                        'interpd.cu', ...
                        }';

                    def.IncludePath = {}; % include folders
                    def.Libraries = {}; % libraries

                    def.CompileOptions = {...  compile options
                        'use_fast_math',...
                        };

                    def.Warnings = {... warnings
                        'no-deprecated-gpu-targets'...
                        };

                    def.DefinedMacros = {};
                case "wbilerp"
                    % filename
                    def.Source = {...
                        'wbilerp.cu', ...
                        }';

                    def.IncludePath = {}; % include folders
                    def.Libraries = {}; % libraries

                    def.CompileOptions = {...  compile options
                        };

                    def.Warnings = {... warnings
                        'no-deprecated-gpu-targets'...
                        };

                    def.DefinedMacros = {};
                case "greens"
                    % filename
                    def.Source = {...
                        'greens.cu', ...
                        }';

                    def.IncludePath = {}; % include folders
                    def.Libraries = {}; % libraries

                    def.CompileOptions = {...  compile options
                        'use_fast_math',...
                        };

                    def.Warnings = {... warnings
                        'no-deprecated-gpu-targets'...
                        };

                    def.DefinedMacros = {};
                case "convd"
                    % filename
                    def.Source = {...
                        ... 'conv_cuda.cu', ...
                        'convd.cu', ...
                        }';

                    def.IncludePath = {}; % include folders
                    def.Libraries = {}; % libraries

                    def.CompileOptions = {...  compile options
                        'use_fast_math',...
                        };

                    def.Warnings = {... warnings
                        'no-deprecated-gpu-targets'...
                        };

                    def.DefinedMacros = {};
            end
            defs(i) = def;
            end
        end
        function defs = genMexdefs(name)
            % GENMEXDEFS - Generate mex compilation definitions.
            % 
            % defs = GENMEXDEFS() returns a struct array of definition
            % structs defs to compile all mex kernels used by QUPS.
            %
            % defs = GENMEXDEFS(name) returns only the kernels specified in
            % the string array name.
            % 
            % defs = GENMEXDEFS("all") returns all available kernels. The 
            % default is "all".
            % 
            % Use with UltrasoundSystem.recompileMex to compile the
            % kernels for the UltrasoundSystem.
            % 
            % Example:
            % us = UltrasoundSystem('recompile', false);  
            % defs = UltrasoundSystem.genMexdefs("msfm2d");
            % us.recompileMex(defs);
            %
            % See also ULTRASOUNDSYSTEM.RECOMPILEMEX 
            arguments
                name (1,:) string {mustBeMember(name, ["all" , ...
                    "msfm2d", "msfm3d" ...
                    ])} = "msfm2d"
            end
            if isscalar(name) && name == "all"
                name = ["msfm2d", "msfm3d"];
            end

            for i = numel(name):-1:1
                switch name(i)
                    case {"msfm2d","msfm3d"}
                        def.Source = char(name(i)+'.c');

                        def.IncludePath = {... include folders
                            fullfile(UltrasoundSystem.getSrcFolder(), 'FMM', 'functions'), ... helper_math
                            };

                        def.Libraries = {...
                            ...
                            };

                        def.DefinedMacros = {...
                            };
                end
                defs(i) = def;
            end
        end
    end
end

function arch = nvarch()
% NVARCH - compute capability architecture string for the current gpu
try
    g = gpuDevice();
    arch = "compute_" + replace(g.ComputeCapability,'.',''); % use the current GPU's CC number
catch
    warning("Unable to access GPU.");
    arch = string.empty; % don't use this argument
end
end

function lat_name = scanlat(scan)
arguments, scan (1,1) Scan, end
if     isa(scan, 'ScanCartesian'),  lat_name = "x";
elseif isa(scan, 'ScanPolar'),      lat_name = "a";
elseif isa(scan, 'ScanGeneric'),    lat_name = "v";
elseif isa(scan, 'ScanSpherical'),  lat_name = "a";
end
end

% validator
function mustBeArch(s)
if ~startsWith(s, "compute_")
    error("Nvidia architectures must start with 'compute_'.");
end
end
