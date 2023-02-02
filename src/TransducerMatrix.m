% TRANSDUCERMATRIX - Matrix Array Transducer class
% 
% A TransducerMatrix defines a matrix array transducer where all elements 
% lie on a plane.
% 
% See also TRANSDUCER TRANSDUCERARRAY TRANSDUCERCONVEX

classdef TransducerMatrix < Transducer

    properties
        pitch (1,2) {mustBeNumeric} = 0.3e-3  % lateral and elevational interelement distance
        numd (1,2) {mustBeNumeric, mustBeInteger}  % number of elements in lateral and elevation dimensions
    end

    properties(Hidden)
        mux_offset (3,:) double = [0;0;0];
    end
        
    % constructor 
    methods(Access=public)
        function self = TransducerMatrix(array_args, xdc_args)
            % TRANSDUCERMATRIX - TransducerMatrix constructor
            %
            % xdc = TRANSDUCERMATRIX(Name, Value, ...) constructs a
            % TransducerMatrix using name/value pair arguments.
            %
            % See also TRANSDUCERARRAY TRANSDUCERCONVEX
            arguments
                array_args.pitch (1,2) double
                array_args.numd (1,2) double
                xdc_args.?Transducer
            end

            % if width/height not set but pitch is, make it the pitch
            if ~isfield(xdc_args, 'width') && isfield(array_args, 'pitch') 
                xdc_args.width = array_args.pitch(1); 
            end
            if ~isfield(xdc_args, 'height') && isfield(array_args, 'pitch') 
                xdc_args.height = array_args.pitch(end); 
            end

            % if numd set, make sure the number of elements 
            % corresponds to the number in each dimension by definition
            if isfield(array_args, 'numd') 
                if isfield(xdc_args, 'numel')
                assert(xdc_args.numel == prod(array_args.numd), ...
                    "The numd (" + (array_args.numd + ",") ...
                    + ") must be factors of the number of elements (" ...
                    + xdc_args.numel + ").");
                else
                    xdc_args.numel = prod(array_args.numd); % set numel to match dimension
                end
            end
            
            % initialize the Transducer
            xdc_args = struct2nvpair(xdc_args);
            self@Transducer(xdc_args{:}) 
            
            % if numd not set, make it a resonably guessed breakdown
            if ~isfield(array_args, 'numd')
                % use product of every other factor -
                % this is a heuristic breakdown that is guaranteed to be
                % the square root if N has a square
                f = factor(self.numel); % factors (sorted)
                array_args.numd = [prod(f(1:2:end)), prod(f(2:2:end))];
            end

            % ensure width/height <= pitch
            % self.width  = min(self.width , self.pitch(1));
            % self.height = min(self.height, self.pitch(1));

            % initialize the TransducerMatrix
            for f = string(fieldnames(array_args))'
                self.(f) = array_args.(f);
            end
        end
    end
    
    % manipulation
    methods
        % scaling
        function self = scale(self, kwargs)
            % SCALE - Scale units
            %
            % xdc = SCALE(xdc, 'dist', factor) scales the distance of the
            % properties by factor. This can be used to convert from meters
            % to millimeters for example.
            %
            % xdc = SCALE(xdc, 'time', factor) scales the temporal
            % properties by factor. This can be used to convert from
            % seconds to microseconds and hertz to megahertz.
            %
            % Example:
            %
            % % Create a TransducerMatrix
            % xdc = TransducerMatrix('fc', 5e6, 'width', 0.3e-3); % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % xdc = scale(xdc, 'dist', 1e3, 'time', 1e6) % mm, us, MHz
            %
            %

            arguments
                self Transducer
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            args = struct2nvpair(kwargs); % get the arguments as Name/Value pairs
            self = scale@Transducer(self, args{:}); % call superclass method
            if isfield(kwargs, 'dist') % scale distance (e.g. m -> mm)
                self.pitch = self.pitch * kwargs.dist;
                self.mux_offset = self.mux_offset * kwargs.dist;
            end
        end
    end
    
    % define position methods
    methods   
        % get methods
        function p = positions(self)
            array_width  = (self.numd(1) - 1) * self.pitch(1);
            array_height = (self.numd(end) - 1) * self.pitch(end);
            x = linspace(-array_width/2,  array_width/2,  self.numd(1));
            y = linspace(-array_height/2, array_height/2, self.numd(end));
            z = 0;
            [x,y,z] = ndgrid(x,y,z);
            p = cat(2, x(:), y(:), z(:))' + self.offset + self.mux_offset;
            % returns a 1 x N vector of the positions of the N elements with 0
            % at the center
        end
        
        function [theta, phi, normal, width, height] = orientations(self)            
            theta = zeros([1, self.numel]);
            phi   = zeros(size(theta));
            ZERO  = zeros(size(theta));
            normal     = [cosd(phi).*sind(theta); sind(phi);  cosd(phi).*cosd(theta)];
            width      = [cosd(theta);            sind(ZERO); -cosd(ZERO).*sind(theta)];
            height     = [sind(phi).*sind(ZERO);  cosd(phi);   sind(phi).*cosd(ZERO)];
        end        
    end

    
    % SIMUS conversion functions
    methods
        % MUST only supports linear arrays or curvilinear arrays
        function p = getSIMUSParam(self), error('MUST does not support matrix arrays.'), end
    end

    % Field II conversion function
    methods
        function aperture = getFieldIIAperture(self, focus, element_sub_divisions)
            if nargin < 2 || isempty(focus), focus = [0 0 realmax('single')]; end % ~ 0-deg plane wave
            if nargin < 3, element_sub_divisions = [1,1]; end % no sub divisions
            % TODO: error if origin not at 0.

            % Field II parameters
            xdc_2d_array_params = { ...
                self.numd(1), ...
                self.numd(end), ...
                self.width, ...
                self.height,...
                abs(self.pitch(1)  ) - self.width, ... kerf in x
                abs(self.pitch(end)) - self.height, ... kerf in y
                ones(self.numd), ...
                element_sub_divisions(1), ...
                element_sub_divisions(end), ...
                reshape(focus, 1, []),...
                };
            
            % ensure double type
            xdc_2d_array_params = cellfun(@double, xdc_2d_array_params, 'UniformOutput', false);
            xdc_2d_array_params = cellfun(@gather, xdc_2d_array_params, 'UniformOutput', false);
            
            % Generate aperture for emission
            try evalc('field_info'); catch, field_init(-1); end
            aperture = xdc_2d_array(xdc_2d_array_params{:});
        end  

        % override
        function p = getFieldIIPositions(xdc)
            % GETFIELDIIPOSITIONS - get an array of element positions
            %
            % p = getFieldIIPositions(xdc) returns a 3 x N vector of the
            % positions of the N elements as represented in FieldII.
            %
            % See also GETFIELDIIAPERTURE FIELD_INFO

            ap = xdc.getFieldIIAperture();
            data = xdc_get(ap, 'rect');
            p = data([24 25 26], :) + xdc.offset;
            xdc_free(ap);
            % p = p(:,1:3:end); % data repeated 3 times for some reason?
        end
    end
    
    % USTB conversion function
    methods
        function probe = getUSTBProbe(self)
            probe = uff.matrix_array(...
                'N_x', self.numd(1), ...
                'N_y', self.numd(end), ...
                'pitch_x', self.pitch(1), ...
                'pitch_y', self.pitch(end), ...
                'element_width', self.width, ...
                'element_height', self.height, ...
                'origin', uff.point('xyz', self.offset(:)') ...
                );
        end
    end

    % Fullwave functions (unsupported)
    methods
        function xdc = getFullwaveTransducer(self, sscan), return; end
    end

    methods
        function set.numd(xdc, numd)
            arguments
                xdc(1,1) TransducerMatrix
                numd (1,2) {mustBeInteger, mustBePositive}
            end
            xdc.numd = numd;
            xdc.numel = prod(numd);
        end
        function n = get.numd(xdc)
            arguments
                xdc (1,1) TransducerMatrix
            end
            assert(isequal(xdc.numel, prod(xdc.numd)), ...
                "The number of elements (" ...
            + xdc.numel ...
            + ") is inconsistent with the number of elements per dimension (" ...
            + xdc.numd(1) + " x " + xdc.numd(end) + " = " + prod(xdc.numd) ...
            + ").");
            n = xdc.numd;
        end
    end

    % dependent variable methods
    methods(Static)
        function xdc = PO192O()
            xdc = TransducerMatrix(...
                'fc', 3.5e6, ...
                'bw_frac', 0.6, ...
                'numd', [32 32], ...
                'width', (0.3e-3), ... 
                'height', (0.3e-3), ...
                'pitch', 0.3e-3, ...
                'el_focus', 20e-3 ...
                );
        end
        function xdc = PO1921()
            xdc = TransducerMatrix(...
                'fc', 7.5e6, ...
                'bw_frac', 0.6, ...
                'numd', [32 32], ...
                'width', (0.3e-3), ... 
                'height', (0.3e-3), ...
                'pitch', 0.3e-3, ...
                'el_focus', 20e-3 ...
                );
        end
        function xdc = Verasonics(Trans, c0)
            arguments
                Trans struct
                c0 (1,1) double = 1540
            end
            % determine the scaling of the properties
            switch Trans.units
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6; % lambda -> m
                case 'mm', scale = 1e-3; % mm -> m
            end

            % infer the number of element in each dimension
            numd = arrayfun(@(i)numel(unique(Trans.ElementPos(:,i))),1:2);
            x = reshape(Trans.ElementPos(:,1), numd);
            z = reshape(Trans.ElementPos(:,2), numd);
            pitch = [mode(mode(diff(x,1,1))), mode(mode(diff(z,1,2)))];

            % for true matrix arrays, create a length property matching the
            % width 
            if ~isfield(Trans, 'elementLength'), Trans.elementLength = Trans.elementWidth; end
            
            % set relevant properties
            xdc = TransducerMatrix(...
                'fc', 1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'bw', 1e6*Trans.Bandwidth([1 end]), ... % bandwidth [Hz]
                'width' , scale*Trans.elementWidth, ... % linear kerf
                'height', scale*Trans.elementLength, ... % Height of element
                'numd', numd, ... % number of elements in each axes
                'pitch', 1e-3*pitch, ... % probe pitch in each axes
                'el_focus', inf ... % elevation focal depth (none)
                );

            % apply mux hardware offset
            pv = scale .* Trans.ElementPos(:,1:2)'; % reported positions
            pv(3,:) = 0; 
            pq = xdc.positions(); % current qups positions
            xdc.mux_offset = pv - pq; % offset from Verasonics definition

        end
        function xdc = UFF(probe)
            arguments
                probe (1,1) {mustBeA(probe, 'uff.matrix_array')}
            end
            if isempty(probe.origin), probe.origin.xyz = [0;0;0]; end % force 0 if empty
            xdc = TransducerMatrix(...
                "numd", [probe.N_x, probe.N_y], ...
                "pitch", [probe.pitch_x, probe.pitch_y], ...
                "width", probe.element_width, ...
                "height", probe.element_height, ...
                "offset", - probe.origin.xyz ...
                );
        end

    end
end
