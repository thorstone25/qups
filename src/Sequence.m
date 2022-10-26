classdef Sequence < matlab.mixin.Copyable
    % SEQUENCE - Class defining transmit sequences
    %
    % A SEQUENCE object defines the parameters for common transmit
    % sequences and is used to define beamforming delays and apodizations
    % per element per pulse. The same waveform must be sent for all
    % elements and pulses. Delays and apodization matrices are generated
    % when given a Transducer.
    %
    % The interpretation of the time-axis of the generated delays and the 
    % foci depend on the Sequence type property. For type 'PW', the foci
    % are normal vectors and time 0 is when the wavefront passes through
    % the spatial origin (i.e. x=y=z=0). For type 'VS', the foci are
    % spatial positions and time 0 is when the wavefront passes through the
    % foci. For type 'FSA', the foci are ignored and time 0 is when the
    % wavefront passes through each element of the given Transducer.
    % 
    % 
    % 
    % See also: SEQUENCERADIAL WAVEFORM
    
    properties
        % TYPE - Type of pulse sequence definition
        %
        % SEQUENCE.TYPE determines the type of pulse sequence. It must be
        % one of {'FSA', 'PW', 'VS'}.
        %
        % When the type is 'FSA', the pulse sequence represents a full 
        % synthetic aperture acquisition where each pulse has one element 
        % transmitting at a time with no delays applied. When using this
        % type, the numpulse property must be set.
        %
        % When the type is 'PW', the pulse sequence represents a plane-wave
        % acquisition where the time delays are applied such that a planar
        % wavefront forms that travels in the direction of the focal
        % vectors and passes through the origin at time 0.
        %
        % When the type is 'VS', the pulse sequence represents a virtual
        % source acquisition where the time delays are applied such that a
        % wavefront forms that travels radially towards and/or away from
        % each foci and passes through the foci at time 0.
        %
        % See also SEQUENCE.FOCUS SEQUENCE.C0 SAEQUENCE.NUMPULSE
        type (1,1) string {mustBeMember(type, ["FSA", "PW", "VS"])} = 'FSA' 
        % FOCUS - Pulse sequence foci or focal vectors
        %
        % SEQUENCE.FOCUS specifies the foci of pulse sequence.
        %
        % When the Sequence type is 'FSA', the focus has no effect.
        %
        % When the Sequence type is 'PW', the foci are unit normal vectors
        % specifying the propagation direction for each plane wave.
        %
        % When the Sequence type is 'VS', the foci are the positions in
        % space where all wavefronts must converge.
        %
        % See also SEQUENCE.TYPE SEQUENCE.C0
        focus (3,:) {mustBeNumeric} = zeros([3,1]);
        % C0 - Reference sound speed
        %
        % SEQUENCE.C0 specifies the sound speed used for creating the
        % delays.
        %
        % See also DELAYS SEQUENCE.FOCUS SEQUENCE.TYPE
        c0 (1,1) {mustBeNumeric} = 1540         
        % PULSE - Transmit pulse waveform
        %
        % SEQUENCE.PULSE specifies the transmit pulse waveform. The default
        % is a Waveform.Delta().
        %
        % See also SEQUENCE.FOCUS SEQUENCE.TYPE SEQUENCE.C0
        pulse (1,1) Waveform = Waveform.Delta() % transmit Waveform
    end
    
    properties(Dependent)
        % NUMPULSE - Number of pulses
        %
        % SEQUENCE.NUMPULSE is the number of transmit pulses used in the
        % sequence. If the Sequence type is 'FSA', this property must be
        % set directly before being used by any function that queries this 
        % property. Otherwise, the number of pulses is inferred from the
        % number of foci.
        %
        % See also SEQUENCE.TYPE
        numPulse (1,1) double {mustBeInteger}
    end
    
    properties(Dependent,Hidden)
        % APODIZATION_ - Manual apodization matrix definition
        %
        % SEQUENCE.APODIZATION_ overrides the default apodization and 
        % specifies the apodization as either a (N x S) matrix, or a
        % function that takes a Transducer and a Sequence and returns a 
        % (N x S) matrix of weights where N is the number of elements and 
        % S is the number of pulses.
        %
        % Example:
        % % Create a hadamard encoding sequence
        % seq = Sequence('type', 'FSA');
        % seq.apodization_ = @(tx, seq) hadamard(tx.numel);
        %
        % % Get a default system and scatterers
        % us = UltrasoundSystem('sequence', copy(seq));
        % scat = Scatterers('pos', [0;0;30e-3]);
        % 
        % % get the ChannelData
        % chdh = greens(us, scat);
        %
        % % decode the data via refocus
        % chd = refocus(us, chdh, 'gamma', 0); % don't use regularization
        % us.sequence.apodization_ = []; % clear the encoding
        % 
        % % beamform and display the image
        % figure; 
        % imagesc(mod2db(DAS(us, chd)));
        % caxis(max(caxis) + [-60 0]);
        % 
        % See also APODIZATION SEQUENCE.DELAYS_
        apodization_ 
        
        % DELAYS_ - Manual delay matrix definition
        %
        % SEQUENCE.DELAYS_ overrides the default delays and 
        % specifies the delays as either a (N x S) matrix, or a
        % function that takes a Transducer and a Sequence and returns a 
        % (N x S) matrix of delays where N is the number of elements and S 
        % is the number of pulses.
        %
        % Example:
        % % Create a random phase sequence
        % seq = Sequence('type', 'FSA');
        % seq.delays_ = @(tx, seq) (randi([0,3],tx.numel) - 1.5) / 4 / tx.fc;
        %
        % % Get a default system and scatterers
        % us = UltrasoundSystem('sequence', copy(seq));
        % scat = Scatterers('pos', [0;0;30e-3]);
        % 
        % % get the ChannelData
        % chdh = greens(us, scat);
        %
        % % decode the data via refocus
        % chd = refocus(us, chdh, 'gamma', 0); % don't use regularization
        % us.sequence.delays_ = []; % clear the encoding
        % 
        % % beamform and display the image
        % figure; 
        % imagesc(mod2db(DAS(us, chd)));
        % caxis(max(caxis) + [-60 0]);
        % 
        % See also DELAYS SEQUENCE.APODIZATION_
        delays_
    end
    
    properties(Hidden,SetAccess=protected)
        FSA_n_tx (1,1) double = nan % hidden storage of the number of pulse for an FSA sequence
        apodizationv_ (:,:) {mustBeNumericOrLogical} = [] % hidden storage of user entered apodization values
        delaysv_      (:,:) {mustBeNumericOrLogical} = [] % hidden storage of user entered delay       values
        apodizationf_ function_handle {mustBeScalarOrEmpty} % hidden storage of user entered apodization values
        delaysf_      function_handle {mustBeScalarOrEmpty} % hidden storage of user entered delay       values
    end
    
    methods
        % constructor
        function self = Sequence(kwargs)
            % SEQUENCE - Sequence constructor
            %
            % seq = SEQUENCE() constructs a Sequence
            %
            % seq = SEQUENCE(Name, Value, ...) uses Name/Value pair
            % arguments to construct a Sequence
            %
            % seq = SEQUENCE('type', 'FSA', 'numPulse', N) defines a full
            % synthetic aperture (FSA) sequence for a transducer with N 
            % elements. The numPulse property must be set explicitly 
            %
            % seq = SEQUENCE('type', 'PW', 'focus', 1 .* [cosd(theta); 0*theta; sind(theta)]) 
            % defines a plane  wave (PW) sequence at the 1 x S array of 
            % angles theta. The norm of the focus should always be 1.
            %
            % seq = SEQUENCE('type', 'VS', 'focus', FOCI) defines a 
            % focused or diverging virtual source (VS) sequence with 
            % focal point locations at the focal points FOCI. FOCI is a
            % 3 x S array of focal points. 
            %
            % seq = SEQUENCE(..., 'c0', c0) sets the beamforming sound 
            % speed to c0. 
            %
            % seq = SEQUENCE(..., 'pulse', wv) defines the transmitted 
            % pulse to be the Waveform wv.
            %
            % See also SEQUENCERADIAL WAVEFORM
            arguments
                kwargs.type (1,1) string {mustBeMember(kwargs.type, ["FSA", "PW", "VS"])}
                kwargs.focus (3,:) double
                kwargs.c0 (1,1) double
                kwargs.pulse (1,1) Waveform
                kwargs.numPulse (1,1) double {mustBeInteger}
            end

            % initialize
            for s = string(fieldnames(kwargs))', self.(s) = kwargs.(s); end
        end

        % scaling
        function self = scale(self, kwargs)
            % SCALE - Scale units
            %
            % seq = SCALE(seq, 'dist', factor) scales the distance of the
            % properties by factor. This can be used to convert from meters
            % to millimeters for example.
            %
            % seq = SCALE(seq, 'time', factor) scales the temporal
            % properties by factor. This can be used to convert from
            % seconds to microseconds and hertz to megahertz.
            %
            % Example:
            %
            % % Create a Sequence
            % seq = Sequence('type', 'VS', 'c0', 1500, 'focus', [0;0;30e-3]); % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % seq = scale(seq, 'dist', 1e3, 'time', 1e6); % mm, us, MHz
            % seq.focus % in mm
            % seq.c0 % in mm/us
            %
            %
            arguments
                self Sequence
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            self = copy(self);
            cscale = 1;
            if isfield(kwargs, 'dist')
                if self.type ~= "PW" % interpret PW as unit vector
                    self.focus = kwargs.dist * self.focus; % scale distance (e.g. m -> mm)
                end
                cscale = cscale * kwargs.dist;
            end
            if isfield(kwargs, 'time')
                self.pulse = scale(self.pulse, 'time', kwargs.time);
                cscale = cscale / kwargs.time;
            end
            self.c0 = self.c0 * cscale;
        end
    end
    
    % conversion methods
    methods
        function seq = getUSTBSequence(self, xdc, t0)
            % GETUSTBSEQUENCE - Get a USTB/UFF uff.sequence object
            %
            % seq = GETUSTBSEQUENCE(self, xdc, t0) creates a USTB
            % compatible sequence object from the QUPS Sequence object
            % where xdc is a QUPS transducer and t0 is the start time in
            % the QUPS coordinate system.
            %
            % See also TRANSDUCER/GETUSTBPROBE
            arguments
                self (1,1) Sequence
                xdc (1,1) Transducer
                t0 (1,1) double % bulk offset of the data
            end

            % initialize all wave objects
            N = self.numPulse;
            for n = N:-1:1, seq(n) = uff.wave(); end
            
            % set the common settings
            [seq.probe] = deal(xdc.getUSTBProbe());
            [seq.sound_speed] = deal(self.c0);
            
            switch self.type
                case {'PW'}
                    [seq.wavefront] = deal(uff.wavefront.plane);
                    theta = atan2(self.focus(1,:),self.focus(3,:));
                    phi   = atan2(self.focus(2,:),hypot(self.focus(1,:),self.focus(3,:)));
                    for n=1:N, seq(n).source = uff.point(...
                            'azimuth', theta(n), ...
                            'elevation', phi(n), ...
                            'distance', inf ...
                            );
                    end
                    [seq.delay] = deal(t0);
                    
                case {'FSA'}
                    p = xdc.positions();
                    [seq.wavefront] = deal(uff.wavefront.spherical);
                    for n=1:N, seq(n).source.xyz = p(:,n).'; end
                    for n=1:N, seq(n).delay = p(:,n)/self.c0 + t0; end
                    
                case {'VS'}
                    [seq.wavefront] = deal(uff.wavefront.spherical);
                    for n=1:N, seq(n).source.xyz = self.focus(:,n).'; end
                    [seq.delay] = deal(t0);
                    
            end   
        end
    end
    
    % temporal response methods
    methods   
        function tau = delays(self, tx)
            % DELAYS - Transmit Sequence delays
            % 
            % tau = DELAYS(self, tx) returns the delays for the Sequence
            % self and transmitting Transducer tx as an N x S array where N
            % is the number of transmitter elements (tx.numel) and S is the
            % number of transmit pulses (self.numPulse).
            %
            % Delays are computed with respect to a time t=0 based on the
            % Sequence type. 
            % 
            % Type:
            %     'VS' : t = 0 when a wave intersects the focus
            %     'PW' : t = 0 when a wave intersects the point [0;0;0]
            %     'FSA': t = 0 when a wave intersects the transmit element
            %
            % If using the plane wave method, the focus is instead
            % interpreted as a normal unit vector. 
            %
            % See also APODIZATION
            arguments
                self Sequence
                tx Transducer
            end
            
            % element positions (3 x 1 x N)
            p = swapdim(tx.positions(),2,3); 
            
            if isempty(self.delays_)
                switch self.type
                    case 'VS'
                        % TODO: use more robust logic for diverging wave test
                        v = self.focus - p; % element to focus vector (3 x S x N)
                        s = ~all(sub(self.focus,3,1) > sub(p,3,1), 3); % whether in behind of the transducer (1 x S)
                        tau = hypot(hypot(sub(v,1,1), sub(v,2,1)),sub(v,3,1)) ./ self.c0; % delay magnitude (1 x S x N)
                        tau = (-1).^s .* tau; % swap sign for diverging transmit

                    case 'PW'
                        % use inner product of plane-wave vector with the
                        % positions to get plane-aligned distance
                        tau = -sum(self.focus .* p, 1) ./ self.c0; % delay (1 x S x N)
                    case 'FSA'
                        tau = zeros([1 size(p,3) size(p,3)]);
                    otherwise
                        error('Reached an unexpected state :(');
                end

                % reshape for output (N x S)
                tau = permute(tau, [3 2 1]);
            else
                if isa(self.delays_, 'function_handle')
                    tau = self.delays_(tx, self); % call the function on tx
                elseif isnumeric(self.delays_)
                    tau = self.delays_; % return the user supplied values
                else, warning("Unable to interpret delays; not a function handle or numeric type.");
                    tau = self.delays_; % return the user supplied values anyway
                end
            end

            
        end

        function a = apodization(self, tx)
            % APODIZATION - Transmit Sequence apodization
            % 
            % a = APODIZATION(self, tx) returns the apodization values for
            % the Sequence self and transmitting Transducer tx as an N x S 
            % array where N is the number of transmitter elements
            % (tx.numel) and S is the number of transmit pulses
            % (self.numPulse).
            %
            % To use a custom apodization, set the hidden property
            % self.apodization_ to be either an array of the proper size 
            % for the Transducer that you plan to use or a function that
            % accepts a Transducer and a Sequence and returns an N x S
            % array of apodization values.
            %
            % Example:
            % % Make a walking aperture sequence
            % xdc = TransducerArray(); % get a transducer
            % pn = xdc.positions(); % 3 x N
            % xn = pn(1,:)'; % (N x 1) receiver x-positions
            % 
            % % define the focal points
            % xf = xdc.pitch * (-10:1:10); % (1 x S) focal x-positions
            % 
            % % define the walking aperture apodization: use only 32
            % elements nearest to the focus.
            % apod = abs(xn - xf) <= 32/2; % (N x S) array
            % 
            % % construct the Sequence
            % seq = Sequence(...
            % 'type', 'VS', ...
            % 'focus', [0;0;30e-3] + xf .* [1;0;0] ...
            % );
            % 
            % % Define the apodization
            % seq.apodization_ = apod; % set the hidden property
            % 
            % See also DELAYS
            arguments
                self Sequence
                tx Transducer
            end

            if isempty(self.apodization_) % apodization not set by user:
                switch self.type
                    case 'FSA'
                        a = eye(size(tx.positions(),2)); % N x N identity
                    otherwise
                        a = ones([size(tx.positions(),2) self.numPulse]); % N x S
                end
            else
                if isa(self.apodization_, 'function_handle')
                    a = self.apodization_(tx, self); % call the function on tx
                elseif isnumeric(self.apodization_)
                    a = self.apodization_; % return the user supplied values
                else, warning("Unable to interpret apodization; not a function handle or numeric type")
                    a = self.apodization_; % return the user supplied values anyway
                end
            end
        end

        function t0 = t0Offset(self)
            arguments, self Sequence, end
            switch self.type
                case 'VS' % for virtual source, t0 is at the foci
                    t0 = - vecnorm(self.focus, 2,1) ./ self.c0; % (1 x S)
                otherwise % PW - t0 is at origin; FSA - t0 at the element
                    t0 = 0; % (1 x 1)
            end
        end
    end
    
    % get methods
    methods
        % number of transmit pulses 
        function v = get.numPulse(self)
            switch self.type
                case 'FSA'
                    v = self.FSA_n_tx;
                otherwise
                    v = size(self.focus, 2);
            end
        end
        
        function set.numPulse(self, n)
            self.FSA_n_tx = n;
        end

        function set.apodization_(self, apod)
            if isa(apod, 'function_handle')
                self.apodizationf_ = apod;
            elseif isnumeric(apod) || islogical(apod)
                self.apodizationv_ = apod;
            else
                error("Expected a function handle or numeric type; instead got a " + class(apod) + ".");
            end
        end
       
        function set.delays_(self, tau)
            if isa(tau, 'function_handle')
                self.delaysf_ = tau;
            elseif isnumeric(tau) || islogical(tau)
                self.delaysv_ = tau;
            else
                error("Expected a function handle or numeric type; instead got a " + class(tau) + ".");
            end
        end

        function apod = get.apodization_(self)
            if ~isempty(self.apodizationf_), 
                apod = self.apodizationf_;
            else, 
                apod = self.apodizationv_;
            end
        end

        function tau = get.delays_(self)
            if ~isempty(self.delaysf_), 
                tau = self.apodizationf_;
            else, 
                tau = self.delaysv_;
            end
        end
    end

    % plotting methods
    methods
        function h = plot(self, varargin, plot_args)
            % PLOT - overload the plot function
            %
            % PLOT(seq) plots the locations of the foci of the Sequence seq
            % based on it's type. 
            % 
            % If the Sequence type is 'PW', the foci are interpreted as
            % normal vectors and are plotted via quiver. Otherwise, the 
            % foci are plotted as points.
            %
            % PLOT(xdc, ax) uses the axes ax instead of the current axes.
            %
            % PLOT(..., Name, Value, ...) passes name-value pair arguments
            % to the built-in plot function so that name value pairs that
            % are valid for plot are valid here.
            %
            % h = PLOT(...) returns the handle to the plot.
            %
            % Plots only the x-z slice.
            %
            % See also TRANSDUCER/PATCH QUIVER


            arguments
                self (1,1) Sequence
            end
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
                plot_args.DisplayName = 'Sequence'
            end
            
            % extract axis and other non-Name/Value pair arguments
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes')
                hax = varargin{1}; varargin(1) = [];
            else 
                hax = gca;
            end

            plot_args = struct2nvpair(plot_args);
            switch self.type
                case{'PW',}
                    % make a quiver plot, starting at the origin, and
                    % pointing in the vector direction
                    [x, y] = deal(zeros([1, self.numPulse]));
                    [u, v] = deal(self.focus(1,:), self.focus(3,:));
                    h = quiver(hax, x, y, u, v, varargin{:}, plot_args{:});
                otherwise
                    % plot the positions with the options given in the inputs
                    h = plot(hax, self.focus(1,:), self.focus(3,:), varargin{:}, plot_args{:});
            end
        end
    end
end
