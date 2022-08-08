classdef SequenceRadial < Sequence
    % SEQUENCERADIAL - Sequence defined by range/angle coordinates
    %
    % A SEQUENCERADIAL object expands the Sequence class by offereing 
    % utilities for transmit sequences that can be defined by range and 
    % angle with respect to an apex.
    %
    % Use with a TransducerConvex with the apex at the center or use to 
    % conveniently create plane wave focal sequences with a range of one
    % across all angles.
    %
    % See also: SEQUENCE TRANSDUCERCONVEX
    
    properties
        apex (3,1) {mustBeNumeric} = [0;0;0] % 3 x 1 center of polar coordinate system
    end
    
    properties(Dependent)
        ranges (1,:) {mustBeNumeric} % range from the apex, or length of the vector (m)
        angles (1,:) {mustBeNumeric} % angle with respect to the z-axis (deg)
    end
    
    methods
        % constructor
        function self = SequenceRadial(varargin, seq_args, seqr_args)
            % SEQUENCERADIAL - SequenceRadial constructor
            %
            % seq = SEQUENCERADIAL() constructs a Sequence
            %
            % seq = SEQUENCERADIAL(Name, Value, ...) uses Name/Value pair
            % arguments to construct a Sequence
            %
            % seq = SEQUENCERADIAL('type', 'PW', 'angles', theta) defines a
            % plane wave (PW) sequence at the 1 x S array of angles theta. 
            % The angles theta must be in degrees.
            %
            % seq = SEQUENCERADIAL('type', 'VS', 'apex', apex, 'ranges', r, 'angles', theta)
            % defines a focused or diverging virtual source (VS) sequence 
            % with focal point locations at the ranges r and angles theta 
            % in polar coordinates with respect to an origin at apex. the
            % angles theta are defined in degrees.
            %
            % seq = SEQUENCERADIAL(..., 'c0', c0) sets the beamforming 
            % sound speed to c0. 
            %
            % seq = SEQUENCERADIAL(..., 'pulse', wv) defines the 
            % transmitted pulse to be the Waveform wv.
            %
            % See also SEQUENCERADIAL WAVEFORM

            arguments(Repeating)
                varargin
            end
            arguments % Sequence arguments
                seq_args.type (1,1) string {mustBeMember(seq_args.type, ["PW", "VS"])} = "PW"
                % seq_args.focus (3,:) double % this object sets the focus
                seq_args.c0 (1,1) double
                seq_args.pulse (1,1) Waveform
            end
            arguments % SequenceRadial arguments
                seqr_args.apex (3,1) {mustBeNumeric} = [0;0;0]
                seqr_args.ranges (1,:) {mustBeNumeric} = 1
                seqr_args.angles (1,:) {mustBeNumeric} = 0
            end
            
            % forward unrecognized inputs to Sequence constructor
            for i = 1:2:numel(varargin), seq_args.(varargin{i}) = varargin{i+1}; end

            % initialize Sequence properties
            seq_args = struct2nvpair(seq_args);
            self@Sequence(seq_args{:}); % initialize with all other N/V pairs
            
            % initialize SequenceRadial properties 
            self.apex = seqr_args.apex; 

            % initialize range/angle together
            [r, a] = deal(seqr_args.ranges, seqr_args.angles);
            if isempty(r) && isempty(a) % both empty - do nothing
            elseif ~isempty(r) && ~isempty(a), % both set - init
                self.setPolar(r, a);
            else % error
                error('Both range and angle must be set together.')
            end
        end
        
        % scaling
        function self = scale(self, kwargs)
            arguments
                self SequenceRadial
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            args = struct2nvpair(kwargs); % gather args
            self = scale@Sequence(self, args{:}); % call the superclass method
            if isfield(kwargs, 'dist')
                self.apex = kwargs.dist * self.apex; % scale distance (e.g. m -> mm)
            end
        end
    end

    % get/set methods
    methods
        function setPolar(self, ranges, angles, apex)
            % SETPOLAR - Set the focal points in polar coordinates
            %
            % SETPOLAR(self, ranges, angles) defines the foci given the  
            % ranges and angles.
            %
            % SETPOLAR(self, ranges, angles, apex) additionally redefines 
            % the apex.
            %
            % See also MOVEAPEX
            arguments % SequenceRadial arguments
                self SequenceRadial
                ranges (1,:) {mustBeNumeric}
                angles (1,:) {mustBeNumeric}
                apex (3,1) {mustBeNumeric} = self.apex;
            end

            self.apex = apex(:); % set apex
            [ranges, angles] = deal(ranges(:)', angles(:)'); % 1 x [1|S]
            foci = ranges .* SequenceRadial.vectors_(angles); % 3 x S
            self.focus = foci + self.apex;
        end
        function r = get.ranges(self), r = vecnorm(self.focus - self.apex, 2, 1); end
        function a = get.angles(self), a = atan2d(sub(self.focus - self.apex,1,1), sub(self.focus - self.apex,3,1)); end
        function set.ranges(self, r), self.focus = self.apex + r .* self.vectors(); end
        function set.angles(self, a), self.focus = self.apex + self.ranges .* SequenceRadial.vectors_(a); end
        function moveApex(self, apex)
            % MOVEAPEX - Set a new apex, preserving the ranges and angles
            %
            % MOVEAPEX(self, apex) moves the apex for the SequenceRadial,
            % preserving the ranges / angles of the Sequence
            %
            % See also: SEQUENCERADIAL/SETPOLAR
            arguments
                self SequenceRadial
                apex (3,1)
            end
            
            % to move the apex, if ranges and angles are there, preserve
            % them for a SEQUENCERADIAL
            [r, a] = deal(self.ranges, self.angles);
            self.setPolar(r, a, apex);
        end
    end

    % helper functions
    methods
        function v = vectors(self), v = SequenceRadial.vectors_(self.angles); end
    end

    % helper functions
    methods(Static, Hidden)
        function v = vectors_(angles), v = cat(1, sind(angles), zeros(size(angles)), cosd(angles)); end
    end

    
    % temporal response methods
    methods           
        function a = apodization(self, tx)
            arguments
                self SequenceRadial
                tx Transducer
            end
            switch self.type
                case 'FSA'
                    % todo: use apodization
                    a = eye(size(tx.positions(),2)); % N x N identity
                otherwise
                    a = ones([size(tx.positions(),2) self.numPulse]); % N x S
            end
        end
        function t0 = t0Offset(self)
            arguments, self SequenceRadial, end
            switch self.type
                case 'VS' % for virtual source, t0 is at the foci
                    % transducer intersects 0: offset to that distance
                    t0 = - (self.ranges - vecnorm(self.apex,2,1)) ./ self.c0; % (1 x S)
                otherwise % PW - t0 is at origin; FSA - t0 at the element
                    t0 = 0; % (1 x 1)
            end
        end
    end

    % plotting functions
    methods
        function h = plot(self, varargin)
            if nargin >= 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes') % determine the axis
                hax = varargin{1}; varargin(1) = []; % delete this input
            else
                hax = gca;
            end

            % make a quiver plot, starting at the origin, and
            % pointing in the vector direction
            vecs = self.vectors() .* self.ranges;
            og = repmat(self.apex, [1,size(vecs,2)]);
            [x, y] = deal(og(1,:), og(3,:));
            [u, v] = deal(vecs(1,:), vecs(3,:));
            h = quiver(hax, x, y, u, v, varargin{:});
        end
    end
        
end