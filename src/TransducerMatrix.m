% TRANSDUCERMATRIX - Matrix Array Transducer class
% 
% A TransducerMatrix defines a matrix array transducer where all elements 
% lie on a plane.
% 
% See also TRANSDUCER TRANSDUCERARRAY TRANSDUCERCONVEX TRANSDUCERPISTON

classdef TransducerMatrix < Transducer

    properties
        pitch (1,2) {mustBeNumeric} = 0.3e-3  % lateral and elevational interelement distance
        numd (1,2) {mustBeNumeric, mustBeInteger}  % number of elements in lateral and elevation dimensions
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

            % if numd set and numel set, make sure the number of elements 
            % corresponds to the number in each dimension by definition
            if isfield(xdc_args, 'numel') && isfield(array_args, 'numd') 
                assert(xdc_args.numel == prod(array_args.numd), ...
                    "The numd (" + (array_args.numd + ",") ...
                    + ") must be factors of the number of elements (" ...
                    + xdc_args.numel + ").");
            end
            
            % initialize the Transducer
            xdc_args = struct2nvpair(xdc_args);
            self@Transducer(xdc_args{:}) 
            
            % if numd not set, make it a resonably guessed breakdown
            if ~isfield(xdc_args, 'numd')
                % use product of every other factor -
                % this is a heuristic breakdown that is guaranteed to be
                % the square root if N has a square
                f = factor(self.numel); % factors (sorted)
                array_args.numd = [prod(f(1:2:end)), prod(f(2:2:end))];
            end

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
            p = cat(2, x(:), y(:), z(:))' + self.offset;
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
        % MUST does not support arrays that aren't linear arrays or
        % curvilinear arrays
        function p = getSIMUSParam(self), error('MUST does not support matrix arrays.'), end
    end

    % Field II conversion function
    methods
        function aperture = getFieldIIAperture(self, focus, element_sub_divisions)
            if nargin < 2 || isempty(focus), focus = [0 0 realmax('single')]; end % ~ 0-deg plane wave
            if nargin < 3, element_sub_divisions = [1,1]; end % no sub divisions
                        
            % Field II parameters
            xdc_2d_array_params = { ...
                self.numd(1), ...
                self.numd(end), ...
                self.width, ...
                self.height,...
                0, ... kerf in x
                0, ... kerf in y
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
            p = p(:,1:3:end); % data repeated 3 times for some reason?
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
        function xdc = getFullwaveTransducer(self, sscan)
            return;

            [dX, dY] = deal(sscan.dx, sscan.dz); % simulation grid step size
            [X0, Y0]= deal(sscan.x(1), sscan.z(1)); % first value
            nX = sscan.size('X' == sscan.order); % grid size
            nY = sscan.size('Z' == sscan.order);
            % map (X, Z) -> (X, Y)
            
            % define variables
            xdc_.npx     = self.numel; % number of elements
            % xdc.thetas  = self.orientations(); % xdc.dTheta*((-(xdc.npx-1)/2):((xdc.npx-1)/2)); % the thetas defining the transmit elements

            % legacy
            % zero_offset = 12.4e-3;      % (deprecated) offset of transducer face, how much it comes into the grid (m)
            % xdc.ptch    = self.pitch; % sind(self.angular_pitch) * self.radius; % angular pitch of the transducer (pixels)
            % xdc.cen     = [(self.offset(1) - X0)/dX, (self.offset(3) - self.radius - Y0)/dY]; % center of the transducer in grid indices

            %% Make incoords and outcoords curves

            % define the thetas at the center of each element
            % evenly space, centered at 0deg
            % for n=1:xdc.npx, xdc.thetas(n)=n*xdc.dTheta; end
            % xdc.thetas = xdc.thetas-mean(xdc.thetas);

            % get x-axis and y-axis
            x = X0 + (0:nX-1) * dX; % 1 x X
            y = Y0 + (0:nY-1) * dY; % 1 x Y

            % Make a rectangle that defines the transducer surface
            pb = self.bounds;
            inmap  = pb(1,1) < x' & x' < pb(1,2) & y < pb(3,2);
            outmap = zeros(nX,nY);

            % Grab the coords on edge of the rectangle - deeper rect for outcoords
            for i=1:nX
                % find inmap coords
                j = find(inmap(i,:)==0);
                j = j(1);
                inmap(i,1:max([j-8 0]))=0; % make a depth of 8-1=7 pixels in y

                % find outmap coords
                j = find(inmap(i,:)==1);
                if(~isempty(j))
                    j = j(end)+2; % offset by 2 indices in y - this is important!
                    outmap(i,j)= 1; % make a depth of 1 pixel in y
                end
            end

            % convert incoords binary map to a vector of coordinates
            xdc_.inmap     = inmap;
            xdc_.incoords  = mapToCoords(double(inmap));
            [~, idcr]     = sort(xdc_.incoords(:,1)); % sort by x instead of y
            xdc_.incoords  = xdc_.incoords(idcr,:);

            % convert outcoords binary map to a vector of coordinates
            xdc_.outcoords = mapToCoords(outmap);
            [~, idcr]     = sort(xdc_.outcoords(:,1)); % sort by x instead of y
            xdc_.outcoords = xdc_.outcoords(idcr,:);

            %% assign which transducer number each incoord is assigned to
            xn = self.positions();
            xn = xn(1,:); % x-positions of the transmit elements

            % get location of the center of each element (in pixels)
            xdc_.outcoords2 = zeros(0,2);
            xdc_.incoords2  = zeros(0,2);

            xdc_.outcoords(:,3) = 1; % this helps with reading genout
            xdc_.outcoords(:,4) = 0; % This labels which tranducer element each subelement is assigned to
            xdc_.incoords (:,4) = 0; % This labels which tranducer element each subelement is assigned to

            for tt=1:xdc_.npx

                % find which incoords are assigned to tt
                % less_than_max    = xdc.thetas_in < (xdc.thetas(tt) + xdc.dTheta/2);
                % greater_than_min = xdc.thetas_in > (xdc.thetas(tt) - xdc.dTheta/2);
                % idtheta = find( less_than_max & greater_than_min);
                idxn = abs(x(xdc_.incoords(:,1)) - xn(tt)) < self.pitch/2; % x is in-bounds
                xdc_.incoords(idxn,4) = tt; % assign

                % find center of tt tx element - do each dim separate cause sometimes idtheta is just one value
                % xdc.incoords2(tt,1) = mean(xdc.incoords(idtheta,1));
                % xdc.incoords2(tt,2) = mean(xdc.incoords(idtheta,2));
                xdc_.incoords2(tt,1:2) = mean(xdc_.incoords(idxn,1:2),1);

                % find which outcoords are assigned to tt
                % less_than_max    = xdc.thetas_out < (xdc.thetas(tt) + xdc.dTheta/2);
                % greater_than_min = xdc.thetas_out > (xdc.thetas(tt) - xdc.dTheta/2);
                % idtheta = find( less_than_max & greater_than_min);
                idxn = abs(x(xdc_.outcoords(:,1)) - xn(tt)) < self.pitch/2; % x is in-bounds
                xdc_.outcoords(idxn,4) = tt; % assign

                % find center of tt rx element - do each dim separate cause sometimes
                % xdc.outcoords2(tt,1) = mean(xdc.outcoords(idtheta,1));
                % xdc.outcoords2(tt,2) = mean(xdc.outcoords(idtheta,2));
                xdc_.outcoords2(tt,1:2) = mean(xdc_.outcoords(idxn,1:2),1);

            end

            xdc_.nOutPx = size(xdc_.outcoords,1);
            xdc_.nInPx  = size(xdc_.incoords,1);

            %     figure(2); clf;
            %     plot(xdc.incoords(:,1),xdc.incoords(:,2),'.'), hold on
            %     plot(xdc.incoords2(:,1),xdc.incoords2(:,2),'.')
            %     plot(xdc.outcoords(:,1),xdc.outcoords(:,2),'.')
            %     plot(xdc.outcoords2(:,1),xdc.outcoords2(:,2),'.'), hold off

            % make vector which labels where the transducer surface is in pixels in
            % y across x
            xdc_.surf = zeros(1,nX);
            for i = 1:nX

                % find where the transducer surface is
                j = find(xdc_.inmap(i,:)==1);
                if(~isempty(j))
                    j = j(end);
                else
                    j = 1;
                end
                xdc_.surf(i) = j + 6; % round(ppw/2); % buffer ?????????????????????????
            end

            % output a struct with only the required fields (for
            % consistency)
            args = cellstr(["npx", "inmap", "nInPx", "nOutPx", "incoords", "outcoords", "incoords2", "outcoords2", "surf"]);
            for a = 1:numel(args), args{2,a} = xdc_.(args{1,a}); end
            xdc = struct(args{:});

        end
    end

    methods
        function set.numd(xdc, numd)
            arguments
                xdc TransducerMatrix
                numd (1,2) {mustBeInteger, mustBePositive}
            end
            xdc.numd = numd;
            xdc.numel = prod(xdc.numd);
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
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6;
                otherwise
                    error('Conversion from Verasonics trans to TransducerMatrix not implemented for units not in wavelengths.');
            end

            % set relevant properties
            xdc = TransducerMatrix(...
                'fc', 1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'bw', 1e6*Trans.Bandwidth([1 end]), ... % bandwidth [Hz]
                'width', scale*Trans.elementWidth, ... % linear kerf
                'height', 1e-3*Trans.elevationApertureMm, ... % Height of element [m]
                'numd', sqrt(Trans.numelements) * [1 1], ... % number of elements in each axes
                'pitch', 1e-3*Trans.spacingMm * [1 1], ... % probe pitch in each axes [m]
                'el_focus', 1e-3*Trans.elevationFocusMm ... % elevation focal depth
                );
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
