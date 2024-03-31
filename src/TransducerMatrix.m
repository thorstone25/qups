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
        function xdc = TransducerMatrix(kwargs)
            % TRANSDUCERMATRIX - TransducerMatrix constructor
            %
            % xdc = TRANSDUCERMATRIX(Name, Value, ...) constructs a
            % TransducerMatrix using name/value pair arguments.
            %
            % See also TRANSDUCERARRAY TRANSDUCERCONVEX
            arguments
                kwargs.?TransducerMatrix
            end
            
            % if only width or height set, then make them match
            if      isfield(kwargs,'width') && ~isfield(kwargs, 'height')
                kwargs.height = kwargs.width;
            elseif ~isfield(kwargs,'width') &&  isfield(kwargs, 'height')
                kwargs.width  = kwargs.height;
            end

            % if width/height not set but pitch is, make it the pitch
            if isfield(kwargs,'pitch')
                if ~isfield(kwargs, 'width'), kwargs.width  = kwargs.pitch(1); end
                if ~isfield(kwargs, 'height'),kwargs.height = kwargs.pitch(end); end
            end

            % if width/height is greater than pitch

            % if numd set, make sure the number of elements 
            % corresponds to the number in each dimension by definition
            if isfield(kwargs, 'numd') 
                if isfield(kwargs, 'numel')
                assert(kwargs.numel == prod(kwargs.numd), ...
                    "The numd (" + (kwargs.numd + ",") ...
                    + ") must be factors of the number of elements (" ...
                    + kwargs.numel + ").");
                else
                    kwargs.numel = prod(kwargs.numd); % set numel to match dimension
                end
            end
            
            % initialize the Transducer
            rmflds = intersect(fieldnames(kwargs), ["numd", "pitch"]);
            args = struct2nvpair(rmfield(kwargs, rmflds));
            xdc@Transducer(args{:}) 
            
            % if numd was not set, make it a resonably guessed breakdown
            if ~isfield(kwargs, 'numd')
                % use product of every other factor -
                % this is a heuristic breakdown that is guaranteed to be
                % the square root if N has a square
                f = factor(xdc.numel); % factors (sorted)
                kwargs.numd = [prod(f(1:2:end)), prod(f(2:2:end))]; % heuristic
            end
            xdc.numd = kwargs.numd;

            % if neither width/height set, make width==height==the smaller of the two
            if ~isfield(kwargs,'width') && ~isfield(kwargs,'height')
                [xdc.width, xdc.height] = deal(min(xdc.width, xdc.height));
            end

            % if pitch was not set, ensure width/height <= pitch
            if isfield(kwargs, 'pitch')
                xdc.pitch = kwargs.pitch;
            elseif ~isfield(kwargs, 'pitch') % && (isfield(kwargs, 'width') || isfield(kwargs, 'height'))
                xdc.pitch = max(xdc.pitch, [xdc.width, xdc.height]);
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
            % % Create a TransducerMatrix
            % xdc = TransducerMatrix('fc', 5e6, 'width', 0.3e-3); % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % xdc = scale(xdc, 'dist', 1e3, 'time', 1e6) % mm, us, MHz
            %
            %

            arguments
                xdc Transducer
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            args = struct2nvpair(kwargs); % get the arguments as Name/Value pairs
            xdc = scale@Transducer(xdc, args{:}); % call superclass method
            if isfield(kwargs, 'dist') % scale distance (e.g. m -> mm)
                xdc.pitch = xdc.pitch * kwargs.dist;
                xdc.mux_offset = xdc.mux_offset * kwargs.dist;
            end
        end
    end
    
    % define position methods
    methods   
        % get methods
        function p = positions(xdc)
            array_width  = (xdc.numd(1) - 1) * xdc.pitch(1);
            array_height = (xdc.numd(end) - 1) * xdc.pitch(end);
            x = linspace(-array_width/2,  array_width/2,  xdc.numd(1));
            y = linspace(-array_height/2, array_height/2, xdc.numd(end));
            z = 0;
            [x,y,z] = ndgrid(x,y,z); % make grid
            p = xdc.transPos([x(:), y(:), z(:)]') + xdc.mux_offset;
        end
        
        function [theta, phi, normal, width, height] = orientations(xdc)            
            theta =  xdc.rot(1) + zeros([1, xdc.numel]);
            phi   = -xdc.rot(2) + zeros(size(theta));
            ZERO  = zeros(size(theta));
            normal     = [cosd(phi).*sind(theta); sind(phi );  cosd(phi ).*cosd(theta)];
            width      = [cosd(theta);            sind(ZERO); -cosd(ZERO).*sind(theta)];
            height     = [sind(phi).*sind(ZERO ); cosd(phi );  sind(phi ).*cosd(ZERO )];
        end        
    end

    
    % SIMUS conversion functions
    methods
        % MUST only supports linear arrays or curvilinear arrays
        function p = getSIMUSParam(xdc), 
            arguments, xdc TransducerMatrix, end
            p = arrayfun(@(xdc) struct( ...
                'fc', xdc.fc, ...
                'elements', sub(xdc.positions(),1:2,1), ...
                'width', xdc.width, ...
                'height', xdc.height, ...
                'bandwidth', 100*xdc.bw_frac ...
                ), xdc);
            if isempty(p), p = struct.empty; end
        end
    end

    % Field II conversion function
    methods
        function aperture = getFieldIIAperture(xdc, sub_div, focus)
            arguments
                xdc TransducerMatrix
                sub_div (1,2) double = [1,1]
                focus (3,1) double = [0 0 realmax('single')]
            end

            focus(isinf(focus)) = realmax('single') .* sign(focus(isinf(focus))); % make focus finite

            % Field II parameters
            xdc_2d_array_params = arrayfun(@(xdc) {{ ...
                xdc.numd(1), ...
                xdc.numd(end), ...
                xdc.width, ...
                xdc.height,...
                abs(xdc.pitch(1)  ) - xdc.width, ... kerf in x
                abs(xdc.pitch(end)) - xdc.height, ... kerf in y
                ones(xdc.numd), ...
                sub_div(1), ...
                sub_div(end), ...
                reshape(focus, 1, []),...
                }}, xdc);
            
            % Generate aperture for emission
            try evalc('field_info'); catch, field_init(-1); end
            i = arrayfun(@(xdc) any(xdc.offset) || any(xdc.rot), xdc); % extra translation/rotation
            aperture( i) = getFieldIIAperture@Transducer(xdc(i), sub_div, focus); % call superclass to make rectangles directly
            aperture(~i) = cellfun(@(p)xdc_2d_array(p{:}), xdc_2d_array_params(~i));
        end
    end
    
    % USTB conversion function
    methods
        function probe = QUPS2USTB(xdc)
            arguments, xdc TransducerMatrix, end
            probe = arrayfun(@(xdc) uff.matrix_array(...
                'N_x', xdc.numd(1), ...
                'N_y', xdc.numd(end), ...
                'pitch_x', xdc.pitch( 1 ), ...
                'pitch_y', xdc.pitch(end), ...
                'element_width', xdc.width, ...
                'element_height', xdc.height, ...
                'origin', uff.point('xyz', xdc.offset(:)') ...
                ), xdc);
            if isempty(probe), probe = reshape(uff.matrix_array.empty, size(probe)); end
        end
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
                Trans (1,1) struct
                c0 (1,1) double = 1540
            end
           
            switch Trans.units % determine the scaling of the properties
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6; % lambda -> m
                case 'mm', scale = 1e-3; % mm -> m
            end

            % infer the number of element in each dimension
            numd = arrayfun(@(i)numel(unique(Trans.ElementPos(:,i))),1:2);
            x = reshape(Trans.ElementPos(:,1), numd);
            z = reshape(Trans.ElementPos(:,2), numd);
            pitch = [mode(mode(diff(x,1,1))), mode(mode(diff(z,1,2)))];

            % for true matrix arrays, create a length property matching the width 
            if ~isfield(Trans, 'elementLength'), Trans.elementLength = Trans.elementWidth; end

            % set relevant properties
            xdc = TransducerMatrix(...
                'fc',       1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'width' ,   scale*Trans.elementWidth, ... % linear kerf
                'height',   scale*Trans.elementLength, ... % Height of element
                'numd',     numd, ... % number of elements in each axes
                'pitch',    scale*pitch ... % probe pitch in each axes
                );

            % apply mux hardware offset
            pv = scale .* Trans.ElementPos(:,1:2)'; % reported positions
            pv(3,:) = 0; 
            pq = xdc.positions(); % current qups positions
            xdc.mux_offset = pv - pq; % offset from Verasonics definition

        end
        function xdc = UFF(probe)
            arguments
                probe uff.matrix_array
            end
            for i = 1:numel(probe)
                if isempty(probe.origin), probe.origin.xyz = [0;0;0]; end % force 0 if empty
            end
            xdc = arrayfun(@(probe) TransducerMatrix(...
                "numd", [probe.N_x, probe.N_y], ...
                "pitch", [probe.pitch_x, probe.pitch_y], ...
                "width", probe.element_width, ...
                "height", probe.element_height, ...
                "offset", - probe.origin.xyz ...
                ), probe);
        end
    end
end
