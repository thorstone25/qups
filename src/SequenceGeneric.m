classdef SequenceGeneric < Sequence
    % SEQUENCEGENERIC - Class defining transmit sequences with arbitrary delays
    %
    % A SEQUENCEGENERIC object defines a generic pulse sequence where an
    % arbitrary set of apodization weight and pulse delays is defined per
    % element per pulse. The same Waveform must be sent for all elements
    % and pulses. Delays and apodization matrices must have compatiable
    % sizing for a given Transducer when used with other methods.
    %
    % See also SEQUENCE SEQUENCERADIAL WAVEFORM

    properties
        % APOD - Apodization function/matrix definition
        %
        % SEQUENCEGENERIC.APOD specifies the apodization weights a as
        % either a (N x S) matrix, or a function that takes a Transducer
        % and a Sequence and returns a (N x S) matrix of weights where N is
        % the number of elements of the corresponding Transducer and S is
        % the number of pulses.
        %
        % Example:
        % % Create a hadamard encoding sequence
        % seq = SequenceGeneric('type', 'FSA');
        % seq.apod = @(tx, seq) hadamard(tx.numel);
        % seq.del  = @(tx, seq) zeros(tx.numel);
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
        % us.seq.apod = @(tx, seq) eye(tx.numel); % specify FSA apodization
        %
        % % beamform and display the image
        % figure;
        % imagesc(mod2db(DAS(us, chd)));
        % caxis(max(caxis) + [-60 0]);
        %
        % See also APODIZATION SEQUENCEGENERIC.DEL
        apod  {mustBeA(apod, ["function_handle", "numeric", "logical"])}

        % DEL - Delay function/matrix definition
        %
        % SEQUENCEGENERIC.DEL specifies the delays as either a (N x S) matrix, or
        % a function that takes a Transducer and a Sequence and returns a
        % (N x S) matrix of delays where N is the number of elements and S
        % is the number of pulses.
        %
        % Example:
        % % Create a random phase sequence
        % del  = @(tx, seq) (randi([0,3],tx.numel) - 1.5) / 4 / tx.fc;
        % apod = @(tx, seq) ones(size(tx.numel));
        % seq = SequenceGeneric('del', del, 'apod', apod);
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
        % us.seq.del  = zeros(chd.N); % set as a FSA sequence
        % us.seq.apod =   eye(chd.N); % set as a FSA sequence
        %
        % % beamform and display the image
        % figure;
        % imagesc(mod2db(DAS(us, chd)));
        % caxis(max(caxis) + [-60 0]);
        %
        % See also DELAYS SEQUENCEGENERIC.APOD
        del   {mustBeA(del,  ["function_handle", "numeric", "logical"])}
    end
    properties(Hidden)
        tscale = 1;
    end
    methods
        % constructor
        function seq = SequenceGeneric(seq_kwargs, gen_kwargs)
            % SEQUENCEGENERIC - SequenceGeneric constructor
            %
            % seq = SEQUENCEGENERIC() constructs a SequenceGeneric
            %
            % seq = SEQUENCEGENERIC(Name, Value, ...) uses Name/Value pair
            % arguments to construct a SequenceGeneric
            %
            % seq = SEQUENCEGENERIC('numPulse', N) defines a generic
            % sequence for a transducer with N pulses.
            %
            % seq = SEQUENCEGENERIC(..., 'c0', c0) sets the beamforming
            % sound speed to c0.
            %
            % seq = SEQUENCEGENERIC(..., 'pulse', wv) defines the
            % transmitted pulse to be the Waveform wv.
            %
            % See also SEQUENCE SEQUENCERADIAL WAVEFORM
            arguments
                seq_kwargs.c0 (1,1) double
                seq_kwargs.pulse (1,1) Waveform
                seq_kwargs.numPulse (1,1) double {mustBeInteger}
                gen_kwargs.apod {mustBeA(gen_kwargs.apod, ["function_handle", "numeric", "logical"])} = @(tx, sq) eye(tx.numel);
                gen_kwargs.del {mustBeA(gen_kwargs.del, ["function_handle", "numeric", "logical"])} = @(tx, sq) zeros(tx.numel);
            end

            % initialize the Sequence
            seq_kwargs.type = 'FSA'; % always use FSA format.
            seq_args = struct2nvpair(seq_kwargs);
            seq@Sequence(seq_args{:});

            % initialize
            for s = string(fieldnames(gen_kwargs))', seq.(s) = gen_kwargs.(s); end
        end

        function s = obj2struct(seq)
            arguments, seq Sequence {mustBeScalarOrEmpty}; end
            s = obj2struct@Sequence(seq); % call superclass conversion
        end

        function seq = scale(seq, kwargs)
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
            % % Create a SequenceGeneric
            % tau = [-2:1:2; zeros(1,5); 2:-1:-2]'; % 5 x 3 delays
            % seq = SequenceGeneric('c0', 1500, 'del', tau, 'apod', ones(size(tau))); % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % seq = scale(seq, 'dist', 1e3, 'time', 1e6); % mm, us, MHz
            % seq.c0 % in mm/us
            % seq.del % in us
            %
            %
            arguments
                seq Sequence
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end

            % call superclass conversion
            args = struct2nvpair(kwargs);
            seq = scale@Sequence(seq, args{:});

            % record scaling for delays
            seq.tscale = seq.tscale .* kwargs.time;
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

            error('Not implemented.');

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

                case {'VS','DV','FC'}
                    [seq.wavefront] = deal(uff.wavefront.spherical);
                    for n=1:N, seq(n).source.xyz = self.focus(:,n).'; end
                    [seq.delay] = deal(t0);

            end
        end
    end

    % temporal response methods
    methods
        function t0 = t0Offset(seq), arguments, seq SequenceGeneric, end, t0 = 0; end

        function tau = delays(seq, tx)
            arguments, seq SequenceGeneric, tx Transducer, end

            % call superclass for delays, then scale time
            tau = seq.tscale .* delays@Sequence(seq, tx);
        end
    end
end
