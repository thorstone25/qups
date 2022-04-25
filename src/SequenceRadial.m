classdef SequenceRadial < Sequence
    % SEQUENCERADIAL - Sequence defined by range/angle coordinates
    %
    % A SEQUENCERADIAL object expands the Sequence class by offereing 
    % utilities for transmit sequences that can be defined by range and 
    % angle with respect to an apex.
    %
    % Use with a TransducerConvex.
    %
    % See also: SEQUENCERADIAL/SEQUENCERADIAL SEQUENCE TRANSDUCERCONVEX
    
    properties
        apex = [0;0;0] % 3 x 1 center of polar coordinate system
    end
    
    properties(Dependent)
        ranges          % range from the apex, or length of the vector (m)
        angles          % angle with respect to the z-axis (deg)
    end
    
    methods
        % constructor
        function self = SequenceRadial(varargin)
            % initialize Sequence properties
            self@Sequence(varargin{:});

            [r, a] = deal([]); % update ranges and angles together

            % initialize SequenceRadial properties
            for i = 1:2:nargin
                switch varargin{i}
                    case {'apex'}
                        self.(varargin{i}) = varargin{i+1};
                    case 'ranges'
                        r = varargin{i+1};
                    case {'angles'}
                        a = varargin{i+1};
                end
            end

            if isempty(r) && isempty(a)
            elseif ~isempty(r) && ~isempty(a),
                self.setPolar(r, a);
            else
                error('Both range and angle must be set together.')
            end
        end
    end

    % get/set methods
    methods
        function setPolar(self, ranges, angles, apex)
            if nargin >= 4, self.apex = apex(:); end % change apex if requested
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
            switch self.type
                case 'FSA'
                    % todo: use apodization
                    a = eye(size(tx.positions(),2)); % N x N identity
                otherwise
                    a = ones([size(tx.positions(),2) self.numPulse]); % N x S
            end
        end
        function t0 = t0Offset(self)
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