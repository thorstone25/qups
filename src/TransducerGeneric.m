% TRANSDUCERGENERIC - Generic Transducer Array class
% 
% A TransducerGeneric defines a generic transducer by the positions and
% orientations of its elements.
% 
% See also TRANSDUCER TRANSDUCERARRAY TRANSDUCERCONVEX TRANSDUCERMATRIX

classdef TransducerGeneric < Transducer

    properties
        pos (3,:) % element positions
        az  (1,:) % element azimuth angle (deg)
        el  (1,:) % element elevation angle (deg)
    end
    
    % constructor 
    methods(Access=public)
        function xdc = TransducerGeneric(gen_args, xdc_args)
            % TRANSDUCERGENERIC - TransducerGeneric constructor
            %
            % xdc = TRANSDUCERGENERIC(Name, Value, ...) constructs a
            % TransducerGeneric using name/value pair arguments.
            %
            % See also TRANSDUCERARRAY TRANSDUCERCONVEX
            arguments
                gen_args.pos (3,:) double = [0;0;0]
                gen_args.az (1,:) double = 0
                gen_args.el (1,:) double = 0
                xdc_args.?Transducer
            end

            % --- check the sizing of the inputs --- %
            N = size(gen_args.pos,2); % number of element positions

            % azimuth angles
            Na = numel(gen_args.az); % get current length
            if Na == 1, gen_args.az = repmat(gen_args.az, [1 N]); end % if scalar, make vector
            assert(any(numel(gen_args.az) == N), 'The number of angles must be either scalar, or equal to the number of positions.');

            % elevation angles
            Ne = numel(gen_args.el); % get current length
            if Ne == 1, gen_args.el = repmat(gen_args.el, [1 N]); end % if scalar, make vector
            assert(any(numel(gen_args.el) == N), 'The number of angles must be either scalar, or equal to the number of positions.');

            % initialize the Transducer
            xdc_args = struct2nvpair(xdc_args);
            xdc@Transducer(xdc_args{:}) 
            
            % initialize the TransducerGeneric 
            for f = string(fieldnames(gen_args))'
                xdc.(f) = gen_args.(f);
            end
        end
    end
    
    % manipulation
    methods
        % scaling
        function xdc = scale(xdc, kwargs)
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
            % % Create a TransducerGeneric
            % xdc = TransducerGeneric('fc', 5e6, 'pos', [0.3e-3; 0; 0.5e-3]); % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % xdc = scale(xdc, 'dist', 1e3, 'time', 1e6) % mm, us, MHz
            % xdc.fc % in MHz
            % xdc.pos % in mm
            %
            %

            arguments
                xdc TransducerGeneric
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            args = struct2nvpair(kwargs); % get the arguments as Name/Value pairs
            xdc = scale@Transducer(xdc, args{:}); % call superclass method
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [xdc.pos] = deal(w*xdc.pos);
            end
        end
    end
    
    % define abstract methods
    methods
        function p = positions(xdc), p = xdc.pos; end
        function [theta, phi, normal, width, height] = orientations(xdc)
            theta = xdc.az;
            phi   = xdc.el;
            ZERO  = zeros(size(theta));
            normal     = [cosd(phi).*sind(theta); sind(phi);  cosd(phi).*cosd(theta)];
            width      = [cosd(phi).*cosd(theta); sind(phi); -cosd(phi).*sind(theta)];
            height     = [sind(phi).*sind(theta); cosd(phi);  sind(phi).*cosd(theta)];
        end
    end

    % Fullwave functions
    methods
        function xdc = getFullwaveTransducer(self, sscan)
            error('Not implemented.');
        end
    end

    % SIMUS conversion functions
    methods
        function p = getSIMUSParam(self)
            error("SIMUS does not support a TransducerGeneric.");
            % TODO: error if origin not at 0.
        end
    end

    % Field II conversion function
    methods(Access=public)
        function aperture = getFieldIIAperture(self, focus, element_sub_divisions)
            error('Not implemented.');
            if nargin < 2 || isempty(focus), focus = [0 0 realmax('single')]; end % ~ 0-deg plane wave
            if nargin < 3, element_sub_divisions = [1,3]; end % no sub divisions
            % TODO: error if origin not at 0.
                        
            % Field II parameters
            xdc_lin_array_params = { ...
                self.numel, ...
                self.width, ...
                self.height,...
                self.kerf, ...
                element_sub_divisions(1), ...
                element_sub_divisions(2), ...
                reshape(focus, 1, []),...
                };
            
            % ensure double type
            xdc_lin_array_params = cellfun(@double, xdc_lin_array_params, 'UniformOutput', false);
            xdc_lin_array_params = cellfun(@gather, xdc_lin_array_params, 'UniformOutput', false);
            
            % Generate aperture for emission
            try evalc('field_info'); catch, field_init(-1); end
            aperture = xdc_linear_array(xdc_lin_array_params{:});
        end        
    end
    
    % USTB conversion function
    methods
        function probe = QUPS2USTB(self)
            probe = uff.probe(...
                'geometry', [ ...
                self.pos; ...
                self.az; ...
                self.el; ...
                repmat([ ...
                self.width; ...
                self.height ...
                ], [1, self.numel]) ...
                ]', ...
                'origin', uff.point('xyz', self.offset(:)') ...
                );
        end
    end

    methods(Static)
        function xdc = UFF(probe)
            arguments
                probe uff.probe
            end
            for i = 1:numel(probe)
                if isempty(probe.origin), probe.origin.xyz = [0;0;0]; end % force 0 if empty
            end
            xdc = arrayfun(@(probe) TransducerGeneric(...
                'pos', [probe.x, probe.y, probe.z]', ...
                'az', probe.theta, ...
                'el', probe.phi, ...
                ... "width", probe.width, ...
                ... "height", probe.height, ...
                ... "numel", probe.N_elements, ...
                "offset", - probe.origin.xyz ...
                ), probe);
        end
    end

    methods(Static)
        function xdc = Verasonics(Trans, c0)
            arguments
                Trans struct
                c0 (1,1) double = 1540
            end
            error('Not implemented.')

            % determine the scaling of the properties
            switch Trans.units
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6; % lambda -> m
                case 'mm', scale = 1e-3; % mm -> m
            end

            % parse the impulse response
            h = Trans.IR1wy; % impulse response
            t0 = - (argmax(hilbert(h))-1) / 250e6; % offset to peak time
            wv = Waveform('t', t0 + (0:numel(h)-1) / 250e6, 'samples',h); % impulse response
            
            % set relevant properties
            xdc = TransducerArray(...
                'fc', 1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'bw', 1e6*Trans.Bandwidth([1 end]), ... % bandwidth [Hz]
                'impulse', wv, ... % impulse response function
                'width', scale*Trans.elementWidth, ... % linear kerf
                'height', 1e-3*Trans.elevationApertureMm, ... % Height of element [m]
                'numel', Trans.numelements, ... % number of elements
                'pitch', 1e-3*Trans.spacingMm, ... % probe pitch [m]
                'el_focus', 1e-3*Trans.elevationFocusMm ... % elevation focal depth
                );
        end
    end
end
