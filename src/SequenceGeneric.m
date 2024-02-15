classdef SequenceGeneric < Sequence
    % SEQUENCEGENERIC - Class defining transmit sequences with arbitrary delays
    %
    % A SEQUENCEGENERIC object defines a generic pulse sequence where an
    % arbitrary set of apdization weight and pulse delays is defined per
    % element per pulse. The same Waveform must be sent for all elements
    % and pulses. Delays and apdization matrices must have compatiable
    % sizing for a given Transducer when used with other methods.
    %
    % See also SEQUENCE SEQUENCERADIAL WAVEFORM

    methods
        % constructor
        function seq = SequenceGeneric(kwargs)
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
                kwargs.c0 (1,1) double
                kwargs.pulse (1,1) Waveform
                kwargs.numPulse (1,1) double {mustBeInteger}
                kwargs.apd {mustBeA(kwargs.apd, ["function_handle", "numeric", "logical"])}
                kwargs.del {mustBeA(kwargs.del, ["function_handle", "numeric", "logical"])}
            end

            % set default apdization / delays based on inputs
            hasfld = isfield(kwargs, ["apd", "del"]);
            if     ~hasfld(1) && ~hasfld(2)
                kwargs.apd = @(tx, seq) eye(tx.numel);
                kwargs.del = @(tx, seq) zeros(tx.numel);
            elseif  hasfld(1) && ~hasfld(2)
                kwargs.del = zeros(size(kwargs.apd));
            elseif ~hasfld(1) &&  hasfld(2)
                kwargs.apd = ones(size(kwargs.del));
            end

            % initialize the Sequence
            kwargs.type = 'FSA'; % always use FSA format.
            seq_args = struct2nvpair(kwargs);
            seq@Sequence(seq_args{:});
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
        function t0 = t0Offset(seq), arguments, seq (1,1) SequenceGeneric, end, t0 = 0; end %#ok<MANU>
    end
end
