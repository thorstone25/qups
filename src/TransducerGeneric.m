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
        function p = positions(xdc), p = xdc.transPos(xdc.pos); end
        function [theta, phi, normal, width, height] = orientations(xdc)
            theta = xdc.az + xdc.rot(1);
            phi   = xdc.el + xdc.rot(2);
            normal = [cosd(phi).*sind(theta); sind(phi);  cosd(phi).*cosd(theta)];
            width  = [cosd(phi).*cosd(theta); sind(phi); -cosd(phi).*sind(theta)];
            height = [sind(phi).*sind(theta); cosd(phi);  sind(phi).*cosd(theta)];
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

            switch Trans.units % determine the scaling of the properties
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6; % lambda -> m
                case 'mm', scale = 1e-3; % mm -> m
            end

            % set relevant properties
            xdc = TransducerGeneric(...
                'fc',       1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'width',    scale*Trans.elementWidth, ... % element width
                'height',   scale*Trans.elementLength, ... % Height of element [m]
                'pos',      scale*Trans.ElementPos(:,1:3)', ... % elements positions
                'az',       rad2deg(Trans.ElementPos(:,4)'), ... % elements aximuth angle
                'el',       rad2deg(Trans.ElementPos(:,5)') ... % elements elevation angle
                );
        end
    end
end
