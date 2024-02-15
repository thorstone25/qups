% SEQUENCE - Class defining transmit sequences
%
% A SEQUENCE object defines the parameters for common transmit
% sequences and is used to define beamforming delays and apodizations
% per element per pulse. The same waveform must be sent for all
% elements and pulses. Delays and apodization matrices are generated
% when given a Transducer using the delays() method.
%
% The interpretation of the foci (i.e. the `focus` property) and the
% time-axis of the generated delays depend on the Sequence type property.
% 
% For type 'FSA' (full synthetic aperture), the foci are ignored and time 0
% is when the wavefront passes through each element of the given
% Transducer.
%
% For type 'PW' (plane waves), the foci are normal vectors and time 0 is
% when the wavefront passes through the spatial origin (i.e. x=y=z=0).
% 
% For types 'FC' (focused), 'DV' (diverging), the foci are spatial
% positions and time 0 is when the wavefront passes through the foci.
% 
% Use type 'FSA' and set the hidden `delays_` and/or `apodization_`
% properties to use custom transmit delays and apodization. These will be
% compatible with all simulation methods.
% 
% 
% See also: SEQUENCERADIAL SEQUENCEGENERIC WAVEFORM
classdef Sequence < matlab.mixin.Copyable & matlab.mixin.Heterogeneous & matlab.mixin.CustomDisplay
    properties
        % TYPE - Type of pulse sequence definition
        %
        % SEQUENCE.TYPE determines the type of pulse sequence. It must be
        % one of {'FSA', 'PW', 'DV', 'FC', 'VS'}.
        %
        % When the type is 'FSA', the pulse sequence represents a full 
        % synthetic aperture acquisition where each pulse has one element 
        % transmitting at a time, and each element transmits at time 0.
        % When using this type, the numpulse property must be set.
        %
        % When the type is 'PW', the pulse sequence represents a plane-wave
        % acquisition where a planar wavefront forms that travels in the 
        % direction of the focal vectors and passes through the origin of
        % the coordinate system at time 0.
        %
        % When the type is 'DV', the pulse sequence represents a diverging
        % wave acquisition where a radial wavefront forms that travels
        % radially away from each foci, starting from each foci at time 0.
        %
        % When the type is 'FC', the pulse sequence represents a focused
        % transmit acquisition where a radial wavefront forms that travels
        % radially towards, through, then away from each foci, passing
        % through each foci at time 0.
        % 
        % The type 'VS' is a legacy representation of a virtual-source that
        % is either a focused or diverging wave transmit. This can be
        % convenient as a placeholder when importing data. It's usage for
        % beamforming is discouraged, as it can be difficult to
        % disambiguate the sign of the beamforming delays based on the
        % geometry of the transducer and foci alone.
        %
        % See also SEQUENCE.FOCUS SEQUENCE.C0 SAEQUENCE.NUMPULSE
        type (1,1) string {mustBeMember(type, ["FSA", "PW", "FC", "DV", "VS"])} = 'FSA'
        % FOCUS - Pulse sequence foci or focal vectors
        %
        % SEQUENCE.FOCUS specifies the foci of pulse sequence.
        %
        % When the Sequence type is 'FSA', the focus has no effect.
        %
        % When the Sequence type is 'PW', the foci are unit normal vectors
        % specifying the propagation direction for each plane wave.
        %
        % When the Sequence type is 'FC', 'DV', or 'VS', the foci are the
        % positions in space where each wavefront (virtually) converges.
        %
        % See also SEQUENCE.TYPE SEQUENCE.C0
        focus (3,:) {mustBeNumeric} = zeros([3,1]);
        % C0 - Reference sound speed
        %
        % SEQUENCE.C0 specifies the sound speed used for computing delays.
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
        % APD - Apodization function/matrix definition
        %
        % SEQUENCE.APD specifies the apodization weights a as
        % either a (N x S) matrix, or a function that takes a Transducer
        % and a Sequence and returns a (N x S) matrix of weights where N is
        % the number of elements of the corresponding Transducer and S is
        % the number of transmit pulses.
        %
        % Example:
        % % Create a hadamard encoding sequence
        % seq = SequenceGeneric('type', 'FSA');
        % seq.apd = @(tx, seq) hadamard(tx.numel);
        % seq.del = @(tx, seq) zeros(tx.numel);
        %
        % % Get a default system and scatterers
        % us = UltrasoundSystem('seq', seq);
        % us.seq.numPulse = us.tx.numel; % set the number of pulses
        % scat = Scatterers('pos', [0;0;30e-3]);
        %
        % % get the ChannelData
        % chdh = greens(us, scat);
        %
        % % decode the data via refocus
        % chd = refocus(us, chdh, 'gamma', 0); % don't use regularization
        % us.seq.apd = @(tx, seq) eye(tx.numel); % specify FSA apodization
        %
        % % beamform and display the image
        % figure;
        % b = DAS(us, chd);
        % imagesc(us.scan, b);
        % caxis(max(caxis) + [-60 0]);
        %
        % See also APODIZATION SEQUENCE.DEL
        apd  {mustBeA(apd, ["function_handle", "numeric", "logical"])}

        % DEL - Delay function/matrix definition
        %
        % SEQUENCE.DEL specifies the delays as either a (N x S) matrix, or
        % a function that takes a Transducer and a Sequence and returns a
        % (N x S) matrix of delays where N is the number of elements and S
        % is the number of pulses.
        %
        % Example:
        % % Create a random phase sequence
        % del = @(tx, seq) (randi([0,3],tx.numel) - 1.5) / 4 / tx.fc;
        % apd = @(tx, seq) ones(size(tx.numel));
        % seq = SequenceGeneric('del', del, 'apd', apd);
        %
        % % Get a default system and scatterers
        % us = UltrasoundSystem('seq', seq);
        % scat = Scatterers('pos', [0;0;30e-3]);
        %
        % % get the ChannelData
        % chdh = greens(us, scat);
        %
        % % decode the data via refocus
        % chd = refocus(us, chdh, 'gamma', 0); % don't use regularization
        % us.seq.del = zeros(chd.N); % set as a FSA sequence
        % us.seq.apd =   eye(chd.N); % set as a FSA sequence
        %
        % % beamform and display the image
        % figure;
        % b = DAS(us, chd);
        % imagesc(us.scan, b);
        % caxis(max(caxis) + [-60 0]);
        %
        % See also DELAYS SEQUENCE.APD
        del   {mustBeA(del,  ["function_handle", "numeric"])}
    end

    properties(Hidden)
        tscale (1,1) double {mustBeReal, mustBeFinite} = 1;
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
        % us = UltrasoundSystem('seq', seq);
        % us.seq.numPulse = us.tx.numel; % set the number of pulses
        % scat = Scatterers('pos', [0;0;30e-3]);
        % 
        % % get the ChannelData
        % chdh = greens(us, scat);
        %
        % % decode the data via refocus
        % chd = refocus(us, chdh, 'gamma', 0); % don't use regularization
        % us.seq.apodization_ = []; % clear the encoding
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
        % seq.delays_ = (randi([0,3],tx.numel) - 1.5) / 4 / tx.fc;
        %
        % % Get a default system and scatterers
        % us = UltrasoundSystem('seq', seq);
        % us.seq.numPulse = us.tx.numel; % set the number of pulses
        % scat = Scatterers('pos', [0;0;30e-3]);
        % 
        % % get the ChannelData
        % chdh = greens(us, scat);
        %
        % % decode the data via refocus
        % chd = refocus(us, chdh, 'gamma', 0); % don't use regularization
        % us.seq.delays_ = []; % clear the encoding
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
            % seq = SEQUENCE('type', 'FC', 'focus', FOCI) or 
            % seq = SEQUENCE('type', 'DV', 'focus', FOCI) defines a 
            % focused (FC) or diverging (DV) virtual source sequence with 
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
                kwargs.type (1,1) string {mustBeMember(kwargs.type, ["FSA", "PW", "VS", "FC", "DV"])}
                kwargs.focus (3,:) double
                kwargs.c0 (1,1) double
                kwargs.pulse (1,1) Waveform
                kwargs.numPulse (1,1) double {mustBeInteger}
                kwargs.apd {mustBeA(kwargs.apd, ["function_handle", "numeric", "logical"])}
                kwargs.del {mustBeA(kwargs.del, ["function_handle", "numeric"])}
            end

            % initialize
            for s = string(fieldnames(kwargs))', self.(s) = kwargs.(s); end
        end

        function s = obj2struct(seq)
            % OBJ2STRUCT - Convert a QUPS object into a native MATLAB struct
            %
            % seq = OBJ2STRUCT(seq) converts the Sequence seq and all of 
            % it's properties into native MATLAB structs.
            %
            % Example:
            %
            % % Create a Sequence
            % seq = Sequence()
            %
            % % convert to a MATLAB struct
            % seq = obj2struct(seq)
            %
            arguments, seq Sequence {mustBeScalarOrEmpty}; end
            W = warning('off', "MATLAB:structOnObject"); % squash warnings
            s = struct(seq); % convert self
            if ~isempty(s), s.pulse = obj2struct(s.pulse); end % convert pulse
            s.class = class(seq); % append class info
            % remove empty hidden fields
            fds  = string(fieldnames(s))'; % all fieldnames
            ids = endsWith(fds, '_') & arrayfun(@(f) isempty(s.(f)), fds); % empty, hidden props
            s = rmfield(s, [fds(ids), "FSA_n_tx"]); % remove empty, hidden & redundant props
            warning(W); % restore warnings
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
            % seq = Sequence('type', 'FC', 'c0', 1500, 'focus', [0;0;30e-3]); % m, s, Hz
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
                if ~isempty(self.del) % manual delays
                if isnumeric(self.del), self.del = self.del * kwargs.time; % scale values
                else, self.tscale = self.tscale .* kwargs.time; % record functional scaling
                end
                    
                end
            end
            self.c0 = self.c0 * cscale;
        end
    end
    
    % conversion methods
    methods
        function sequence = QUPS2USTB(seq, xdc, t0)
            % QUPS2USTB - Get a USTB compatible uff.wave object array
            %
            % sequence = QUPS2USTB(seq, xdc, t0) creates a USTB compatible
            % uff.wave array from the QUPS Sequence object seq where xdc is
            % a Transducer and t0 is the start time in the QUPS coordinate
            % system.
            %
            % See also TRANSDUCER.QUPS2USTB CHANNELDATA.QUPS2USTB
            arguments
                seq (1,1) Sequence
                xdc (1,1) Transducer
                t0 (1,:) double = 0 % bulk offset of the data
            end

            % initialize all wave objects
            N = seq.numPulse;
            for n = N:-1:1, sequence(n) = uff.wave(); end

            % set the common settings
            [sequence.probe] = deal(xdc.QUPS2USTB());
            [sequence.sound_speed] = deal(seq.c0);
            
            switch seq.type
                case {'PW'}
                    [sequence.wavefront] = deal(uff.wavefront.plane);
                    theta = atan2(seq.focus(1,:),seq.focus(3,:));
                    phi   = atan2(seq.focus(2,:),hypot(seq.focus(1,:),seq.focus(3,:)));
                    for n=1:N, sequence(n).source = uff.point(...
                            'azimuth', theta(n), ...
                            'elevation', phi(n), ...
                            'distance', inf ...
                            );
                    end
                    
                case {'FSA'}
                    p = xdc.positions();
                    [sequence.wavefront] = deal(uff.wavefront.spherical);
                    for n=1:N, sequence(n).source.xyz = p(:,n).'; end
                    t0 = t0 + vecnorm(p,2,1) ./ seq.c0; % delay transform from element to origin for FSA
                    
                case {'VS', 'DV', 'FC'} % focused and diverging wave
                    [sequence.wavefront] = deal(uff.wavefront.spherical);
                    for n=1:N, sequence(n).source.xyz = seq.focus(:,n).'; end
                    switch seq.type % transform for focal point to origin
                        case 'DV', t0 = t0 - vecnorm(seq.focus,2,1) ./ seq.c0;
                        case 'FC', t0 = t0 + vecnorm(seq.focus,2,1) ./ seq.c0;
                        case 'VS', t0 = t0 + vecnorm(seq.focus,2,1) ./ seq.c0;
                            warning("QUPS:QUPS2USTB:ambiguousSequence", ...
                                "A Sequence of type 'VS' (virtual source) is ambiguous and will be treated as a focused transmit: " ...
                                + "set the type to 'FC' or 'DV' to avoid this warning." ...
                                );                                   
                    end
            end   

            % set the start time
            t0 = num2cell(t0);
            [sequence.delay] = deal(t0{:});
        end
    end

    methods(Static)
        function seq = UFF(sequence, c0)
            % UFF - Create a Sequence from a uff.wave object array
            %
            % seq = UFF(sequence) creates a
            % Sequence seq from the
            % uff.wave object array sequence.
            %
            % seq = UFF(sequence, c0) additionally sets the
            % Scan us.scan from the uff.scan uscan.
            %
            % See also SEQUENCE.UFF

            arguments
                sequence (1,:) uff.wave
                c0 {mustBeReal, mustBeScalarOrEmpty} = uniquetol([sequence.sound_speed])
            end

            wvt = unique([sequence.wavefront]);
            if(~isscalar(wvt)), error( ...
                    'QUPS:Sequence:nonUniqueWavefront', ...
                    'The uff.channel_data object must contain a unique wavefront type.'...
                    );
            end

            prb = unique([sequence.probe]);
            if(~isscalar(prb)), error( ...
                    'QUPS:Sequence:nonUniqueProbe', ...
                    'The uff.channel_data must contain a unique probe.'...
                    );
            end
            
            % get the list of sources
            p0 = [sequence.source];

            % get the sequence type
            switch wvt 
                case uff.wavefront.plane
                    type = 'PW'; 
                
                case uff.wavefront.spherical
                    pn = gather(single(prb.geometry(:,1:3)));
                    pm = gather(single(cat(1, p0.xyz)));
                    pd = [p0.distance];
                    if(~(all(isinf(pd)) || all(~isinf(pd)))), error( ...
                        'QUPS:Sequence:nonUniqueWavefront', ...
                        'The uff.channel_data object contains heterogeneous wavefront types.'...
                        );
                    end
                    if all(isinf(pd)), type = 'PW';
                    elseif isalmostn(pn, pm), type = 'FSA';
                    else, type = 'VS';
                    end
                otherwise
                    error( ...
                        'QUPS:Sequence:unknownWavefront', ...
                        "The uff.channel_data object uses unrecognized wavefront type '" + wvt + "'." ...
                        );
            end

            switch type
                case {'FC','DV','VS'}
                    seq = Sequence('type', type, 'c0', c0, 'focus', cat(1,p0.xyz)');
                case 'FSA'
                    seq = Sequence('type', type, 'c0', c0, 'numPulse', numel(sequence));
                case 'PW'
                    th  = [p0.azimuth];
                    phi = [p0.elevation];
                    nf = [sin(th).*cos(phi); cos(0).*sin(phi); cos(th).*cos(phi)];
                    seq = SequenceRadial('type', type, 'c0', c0, 'focus', nf, 'apex', [0;0;0]);     
            end

        end
    
        function [seq, t0] = Verasonics(TX, Trans, TW, kwargs)
            % VERASONICS - Construct a Sequence from Verasonics structs
            %
            % seq = Sequence.Verasonics(TX, Trans) constructs a Sequence
            % from the Verasonics 'TX' and 'Trans' structs.
            %
            % seq = Sequence.Verasonics(TX, Trans, TW) additionally imports
            % the trilevel excitation waveform from the 'TW' struct as a
            % Waveform. If omitted or empty, seq.pulse is a Waveform.Delta
            % instead. 
            % 
            % [seq, t0] = Sequence.Verasonics(...) additionally returns an
            % offset time array t0 between the transmit delay conventions
            % used by QUPS and Verasonics. If the delays cannot be
            % validated, t0 will be NaN.
            % 
            % [...] = Sequence.Verasonics(..., 'tol', tol) sets the numeric
            % threshold for verifying the parsed Verasonics delays are
            % equivalent to the delays used by QUPS. The default is 1e-16.
            % 
            % [...] =  Sequence.Verasonics(..., 'c0', c0) sets the sound
            % speed c0 in m/s. This should match the value of the Versonics
            % variable 'Resource.Parameters.speedOfSound'. The default is
            % 1540.
            % 
            % [...] =  Sequence.Verasonics(..., 'aperture', ap) explicitly
            % provides the element to channel mapping. The default is
            % 'Trans.HVMux.Aperture' if TX has an 'aperture' property or
            % 'Trans.ConnectorES' otherwise.
            % 
            % [...] =  Sequence.Verasonics(..., 'xdc', xdc) provides a
            % Transducer for verifying the parsed delays. This is helpful
            % if you have a custom transducer or if Transducer.Verasonics
            % fails to import the Trans struct properly.
            %
            % Example:
            % 
            % % get the reference sound speed in meters / second
            % c0 = Resource.Parameters.speedOfSound;
            % 
            % % import the Sequence and delay offsets
            % [seq, t0] = Sequence.Verasonics(TX, Trans, TW, 'c0', c0);
            % 
            % See also TRANSDUCER.VERASONICS WAVEFORM.VERASONICS SCAN.VERASONICS
            arguments
                TX struct
                Trans (1,1) struct
                TW struct {mustBeScalarOrEmpty} = struct.empty
                kwargs.c0 (1,1) {mustBeNumeric, mustBePositive, mustBeFloat} = 1540
                kwargs.tol (1,2) {mustBeNumeric, mustBePositive, mustBeFloat} = 1e-16 
                kwargs.aperture (:,:) {mustBeNumeric, mustBeInteger, mustBePositive}
                kwargs.xdc Transducer {mustBeScalarOrEmpty} = TransducerArray.empty
            end
            
            % get channel mapping
            ismux = isfield(TX, 'aperture'); % whether muxed or not
            if isfield(kwargs, 'aperture')
                            ap = kwargs.aperture; % custom muxing
            elseif ismux,   ap = Trans.HVMux.Aperture; % mux
            else,           ap = Trans.ConnectorES; % no muxing
            end

            % constants
            tol = kwargs.tol;
            c0 = kwargs.c0;
            fc = 1e6*Trans.frequency;
            lambda = c0 / fc;

            % tx params
            apd = cat(1,TX.Apod  ); % apodization
            ang = cat(1,TX.Steer ); % angles
            rf  = cat(1,TX.focus );% .* lambda; % focal range (lambda)
            pog = cat(1,TX.Origin);% .* lambda; % beam origin (lambda)
            tau = cat(1,TX.Delay ) ./ fc; % tx delays (s)

            % build the full delay/apodization matrix
            [apdtx, tautx] = deal(zeros(Trans.numelements, numel(TX))); % pre-allocate
            for i = 1 : numel(TX) % for each transmit
                if ismux,   api = ap(:, TX(i).aperture); 
                else,       api = ap(logical(ap)); 
                end %       selected aperture
                j = logical(api); % active elements
                apdtx(j,i) = apd(i,:); % apodization
                tautx(j,i) = tau(i,:); % delays
            end

            % attempt to import the Transducer
            xdc = kwargs.xdc;
            if isempty(xdc)
            try    xdc = Transducer.Verasonics(Trans, c0);
            catch, xdc = TransducerGeneric.empty; 
            end
            end

            % virtual source ambiguous sequence type warning
            [wid, wmsg] = deal( ...
                "QUPS:Verasonics:ambiguousSequenceType", ...
                "Cannot infer whether sequence is focused or diverging." ...
                );

            % create the corresponding Sequence
            if isfield(TX, "FocalPt") % focal points -> VS
                pf = cat(1, TX.FocalPt)' .* lambda; % focal points
                % attempt to infer focused or diverging wave
                if     isa(class(xdc), "TransducerArray" ) ...
                        || isa(class(xdc), "TransducerMatrix")
                    if     all(pf(3,:) < 0), styp = "DV";
                    elseif all(pf(3,:) > 0), styp = "FC";
                    end
                elseif isa(class(xdc), "TransducerConvex")
                    r = vecnorm(pf - xdc.center,2,1);
                    if     all(r < xdc.radius), styp = "DV";
                    elseif all(r > xdc.radius), styp = "FC";
                    end
                end
                if ~exist('styp', 'var') % fallback to virtual source
                    warning(wid, wmsg); % VS ambiguity warning
                    styp = "VS"; % default type
                end
                seq = Sequence("type",styp, "focus", pf);
            
            elseif ~any(tau,'all') % no delays -> FSA
                seq = Sequence("type","FSA", "numPulse",numel(TX));
            
            elseif all(all(pog == 0,2) & all(rf == 0,1) & any(ang,'all'),1) % PW
                az = rad2deg(ang(:,1)'); % azimuth
                el = rad2deg(ang(:,2)'); % elevation
                if any(el)
                    seq = SequenceSpherical("type","PW","angles",[az; el]);
                else
                    seq = SequenceRadial(   "type","PW","angles", az     );
                end
            
            elseif any(rf)
                pf = pog + rf .* [
                    sin(ang(:,1)) .* cos(ang(:,2)), ...
                          1       .* sin(ang(:,2)), ...
                    cos(ang(:,1)) .* cos(ang(:,2)), ...
                    ]; % focal points
                if     all(rf > 0), styp = "FC"; % focused
                elseif all(rf < 0), styp = "DV"; % diverging
                else,               styp = "VS"; % unclear
                                    warning(wid, wmsg); % VS ambiguity warning
                end
                seq = Sequence("type",styp, "focus", lambda * pf.');
            
            else
                warning( ...
                    "QUPS:Verasonics:ambiguousSequenceType", ...
                    "Unable to infer transmit sequence type." ...
                    );
                seq = SequenceGeneric("apd",apdtx, "del",tautx, "numPulse",numel(TX));
            end
            seq.c0 = c0;

            % validate the apodization and override if necessary
            val = ~isempty(xdc) && isalmostn(apdtx, seq.apodization(xdc), tol(end));
            if ~val
                    warning(...
                        "QUPS:Verasonics:overrideSequenceApodization", ...
                        "Overriding QUPS apodization with Vantage defined values." ...
                        );
                    seq.apodization_ = apdtx;
            end

            % validate the delays
            [t0, val] = deal(NaN, false); % no offset / unverified until proven successful
            if ~isempty(xdc) % transducer successfully imported
                act  = logical(apdtx); % whether elements were active
                tauq = seq.delays(xdc); % QUPS delays (N x S)
                tauv = -tautx;          % VSX  delays (N x S)
                if isequal(size(tauq), size(tauv)) % sizing matches (proceed)
                    [tauq(~act), tauv(~act)] = deal(nan); % set inactive delays to NaN
                    t0 = mean(tauv - tauq,1,'omitnan'); % offset time per transmit
                    val = isalmostn(tauv, tauq + t0, tol(1)); % verified with offset
                end
            end

            % override delays if they don't match 
            if ~val
                warning(...
                    "QUPS:Verasonics:overrideSequenceDelays", ...
                    "Overriding QUPS delays with Vantage defined values." ...
                    );
                seq.delays_ = tautx;
            end

            % import waveform
            if isscalar(TW), seq.pulse = Waveform.Verasonics(TW, fc); 
            else,            seq.pulse = Waveform.Delta();
            end
        end
    end
    
    % transmit sequence generator functions
    methods(Static)
        function apd = apWalking(N, sz, str, off)
            % APWALKING - Generate a walking active aperture
            %
            % apd = apWalking(N, sz) generates a walking active aperture
            % apodization array apd for a Transducer with N elements with
            % sz contiguous elements active per transmit. This method is
            % valid for 1D Transducers with contiguous elements e.g. a
            % TransducerArray or TransducerConvex.
            % 
            % apd = apWalking(N, sz, str) additionally sets the stride str
            % between active apertures. The default is 1.
            % 
            % apd = apWalking(N, sz, str, off) additionally sets an offset
            % for the first active aperture. The default is 0.
            %
            % Example:
            % % Create a Transducer
            % xdc = TransducerArray.L12_3v(); % 192 elements
            % 
            % % Generate active apertures of 64 elements moving by 4
            % apd = Sequence.apWalking(xdc.numel, 64, 4);
            %
            % % Create corresponding focal positions
            % zf  = 50e-3; % focal depth
            % pf  = xdc.focActive(apd, zf);
            % 
            % % Create a focused Sequence
            % seq = Sequence('type','FC', 'focus',pf);
            % seq.apodization_ = apd;  % set apodization
            % 
            % % Create and plot the system
            % us = UltrasoundSystem('xdc', xdc, 'seq', seq);
            % plot(us);
            % 
            % See also SEQUENCE.APODIZATION_ TRANSDUCER.FOCACTIVE


            arguments
                N (1,1) % number of total elements
                sz (1,1) {mustBePositive, mustBeInteger} % active size
                str (1,1) {mustBePositive, mustBeInteger} = 1 % stride
                off (1,1) {mustBeNonnegative, mustBeInteger} = 0 % starting offset
            end
            apd = cell2mat(arrayfun(@(i){ ...
                circshift((1:N)' <= sz, i) ...
                }, off : str : max(off,N-sz)));
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
            %     'PW' : t = 0 when a wave intersects the point [0;0;0]
            %     'FSA': t = 0 when a wave intersects the transmit element
            %     'FC' : t = 0 when a wave intersects the focus
            %     'DV' : t = 0 when a wave intersects the focus
            %     'VS' : t = 0 when a wave intersects the focus
            %
            % If using the plane wave method, the focus is instead
            % interpreted as a normal unit vector. 
            %
            % See also APODIZATION
            arguments
                self (1,1) Sequence
                tx (1,1) Transducer
            end
            
            % element positions (3 x 1 x N)
            p = swapdim(tx.positions(),2,3); 
            
            if isempty(self.del)
                switch self.type
                    case {'FC','DV','VS'}
                        v = self.focus - p; % element to focus vector (3 x S x N)
                        tau = hypot(hypot(sub(v,1,1), sub(v,2,1)),sub(v,3,1)) ./ self.c0; % delay magnitude (1 x S x N)
                        switch self.type % get sign swap
                            case 'VS', s = (-1).^(~all(sub(self.focus,3,1) > sub(p,3,1), 3)); % whether behind the transducer (1 x S)
                            case 'FC', s = +1; % positive delays
                            case 'DV', s = -1; % negate delays
                        end
                        tau = tau .* s; % swap sign for diverging transmit
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
                if isa(self.del, 'function_handle')
                    tau = self.del(tx, self); % call the function on tx
                    tau = tau .* self.tscale; % scale 
                elseif isnumeric(self.del)
                    tau = self.del; % return the user supplied values
                else, warning("Unable to interpret delays; not a function handle or numeric type.");
                    tau = self.del; % return the user supplied values anyway
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
            % To use a custom apodization, set the property self.apd to be
            % either an array of the proper size for the Transducer that
            % you plan to use or a function that accepts a Transducer and a
            % Sequence and returns an N x S array of apodization values.
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
            % apd = abs(xn - xf) <= 32/2; % (N x S) array
            % 
            % % construct the Sequence
            % seq = Sequence(...
            % 'type', 'FC', ...
            % 'focus', [0;0;30e-3] + xf .* [1;0;0] ...
            % );
            % 
            % % Define the apodization
            % seq.apodization_ = apd; % set the hidden property
            % 
            % See also DELAYS
            arguments
                self (1,1) Sequence
                tx (1,1) Transducer
            end

            if isempty(self.apd) % apodization not set by user:
                switch self.type
                    case 'FSA'
                        a = eye(size(tx.positions(),2)); % N x N identity
                    otherwise
                        a = ones([size(tx.positions(),2) self.numPulse]); % N x S
                end
            else
                if isa(self.apd, 'function_handle')
                    a = self.apd(tx, self); % call the function on tx
                elseif isnumeric(self.apd) || islogical(self.apd)
                    a = self.apd; % return the user supplied values
                else, warning("Unable to interpret apodization; not a function handle or numeric type")
                    a = self.apd; % return the user supplied values anyway
                end
            end
        end

        function t0 = t0Offset(seq)
            % T0OFFSET - Compute the start time offset to the origin
            %
            % t0 = t0Offset(seq) computes the start time offset t0 for the
            % Sequence seq. For FSA and PW sequences, t0 is always 0. For
            % virtual source sequences, it can be used to shift the spatial
            % location of t0 from the foci to the origin of the coordinate
            % system.
            % 
            % Example: 
            % % get a default system
            % us = UltrasoundSystem();
            % 
            % % create a focused Sequence and a scatterer at the focus
            % us.seq = Sequence('type', 'FC', 'focus', [0,0,30e-3]', 'c0', 1500);
            % scat = Scatterers('pos', us.seq.focus,'c0',us.seq.c0);
            %
            % % get channel data for a scatterrer at the focus
            % % t0 == 0 for this data starts at the focus, at [0,0,30e-3]
            % chd = greens(us, scat);
            %
            % % shift t0 == 0 to start at the origin
            % t0off = swapdim(us.seq.t0Offset(), 2, chd.mdim); % ensure transmits are across the matching dimension
            % chd_og = copy(chd);
            % chd_og.t0 = chd_og.t0 - t0off;
            % 
            % % show the data with each time scale
            % figure;
            % imagesc(chd   , 1, nexttile()); title('Default axis');
            % imagesc(chd_og, 1, nexttile()); title('Shifted axis');
            % 
            % See also SEQUENCE/TYPE

            arguments, seq (1,1) Sequence, end
            switch seq.type
                case {'VS', 'FC'} % for virtual source, t0 is at the foci
                    t0 = - vecnorm(seq.focus, 2,1) ./ seq.c0; % (1 x S)
                case {'DV'} % for virtual source, t0 is at the foci
                    warning("Untested: please verify this code.");
                    t0 = + vecnorm(seq.focus, 2,1) ./ seq.c0; % (1 x S)
                case {'FSA', 'PW'} % PW - t0 is at origin; FSA - t0 at the element
                    t0 = 0; % (1 x 1)
            end
        end
    end
    
    % dependent properties
    methods
        % number of transmit pulses 
        function v = get.numPulse(self)
            [asz, tsz] = deal(size(self.apodizationv_), size(self.delaysv_)); % override data sizing
            if any([asz, tsz])
                if     ~any(tsz), v = asz(2);
                elseif ~any(asz), v = tsz(2);
                else
                    if ~any(asz ~= tsz), v = asz(2);
                    else, error("QUPS:Sequence:matrixSizeMismatch",...
                            "Expected the delay and apodization matrices to have identical sizing."),
                    end
                end
            else
                switch self.type
                    case 'FSA'
                        v = self.FSA_n_tx;
                        if isempty(v) || ~isfinite(v), warning("Number of pulses is unset."); end
                    otherwise
                        v = size(self.focus, 2);
                end
            end
        end
        function set.numPulse(self, n)
            arguments, self (1,1) Sequence, n (1,1) double, end
            self.FSA_n_tx = n;
        end

        function set.apd(seq, apd)
            arguments
                seq (1,1) Sequence
                apd {mustBeA(apd, ["function_handle", "numeric", "logical"])}
            end
            if isa(apd, 'function_handle')
                [seq.apodizationv_, seq.apodizationf_] = deal([], apd);
            else
                [seq.apodizationv_, seq.apodizationf_] = deal(apd, function_handle.empty);
            end
        end
        function apd = get.apd(seq)
            arguments, seq (1,1) Sequence, end
            if     ~isempty(seq.apodizationf_), apd = seq.apodizationf_;
            elseif ~isempty(seq.apodizationv_), apd = seq.apodizationv_;
            else, apd = [];
            end
        end
        function set.del(seq, tau)
            arguments
                seq (1,1) Sequence
                tau {mustBeA(tau, ["function_handle", "numeric"])}
            end
            if isa(tau, 'function_handle')
                [seq.delaysv_, seq.delaysf_] = deal([], tau);
            else
                [seq.delaysv_, seq.delaysf_] = deal(tau, function_handle.empty);
            end
        end
        function tau = get.del(seq)
            arguments, seq (1,1) Sequence, end
            if     ~isempty(seq.delaysf_), tau = seq.delaysf_;
            elseif ~isempty(seq.delaysv_), tau = seq.delaysv_;
            else, tau = [];
            end
        end

        function set.apodization_(self, apd)
            warning("QUPS:Sequence:DeprecatedProperty", ...
                "The 'apodization_' property is deprecated - use 'apd' instead." ...
                );
            self.apd = apd;
        end
        function apd = get.apodization_(self)
            warning("QUPS:Sequence:DeprecatedProperty", ...
                "The 'apodization_' property is deprecated - use 'apd' instead." ...
                );
            apd = self.apd;
        end
        function set.delays_(self, tau)
            warning("QUPS:Sequence:DeprecatedProperty", ...
                "The 'delays_' property is deprecated - use 'del' instead." ...
                );
            self.del = tau;
        end
        function tau = get.delays_(self)
            warning("QUPS:Sequence:DeprecatedProperty", ...
                "The 'delays_' property is deprecated - use 'del' instead." ...
                );
            tau = self.del;
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

    % object display
    methods(Access = protected)
        function propgrp = getPropertyGroups(seq)
            if ~isscalar(seq)
                propgrp = getPropertyGroups@matlab.mixin.CustomDisplay(seq);
            else
                p = string(properties(seq)); % get public properties
                for f = ["apd", "del"] % fields to remove if empty
                if isprop(seq, f) && isempty(seq.(f)), p(p == f) = []; end
                end
                propgrp = matlab.mixin.util.PropertyGroup(p);
            end
        end
    end
end
