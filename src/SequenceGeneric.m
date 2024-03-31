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
            % Example:
            % N   = 64; % elems
            % apd = Sequence.apWalking(N, N/2, N/8); % elems x transmits
            % del = 0 * apd; % no delays
            % M   = size(apd, 2); % number of transmits
            % seq = SequenceGeneric('apd', apd, 'del', del, 'numPulse', M)
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
end
