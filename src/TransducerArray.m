% TRANSDUCERARRAY - Linear Array Transducer class
% 
% A TransducerArray defines a linear transducer by it's pitch (a.k.a.
% element spacing) and number of elements.
% 
% See also TRANSDUCER TRANSDUCERCONVEX TRANSDUCERPISTON

classdef TransducerArray < Transducer

    properties
        pitch               % the interelement distance (m)
    end
    
    properties(Dependent)
        kerf                % the spacing between elements (m)
        aperture_size       % size of the aperture
    end
    
    % constructor and get/set methods
    methods(Access=public)
        % constructor: accept name/value pairs
        function self = TransducerArray(array_args, xdc_args)
            % TRANSDUCERARRAY - TransducerArray constructor
            %
            % xdc = TRANSDUCERARRAY(Name, Value, ...) constructs a
            % TransducerArray using name/value pair arguments.
            %
            % xdc = TRANSDUCERARRAY(uff_probe) construct a linear array
            % from the uff.linear_array uff_probe.
            %
            % See also TRANSDUCERCONVEX
            arguments
                array_args.pitch double {mustBeScalarOrEmpty} 
                array_args.kerf double {mustBeScalarOrEmpty} 
                xdc_args.?Transducer
            end
            
            % initialize the Transducer
            xdc_args = struct2nvpair(xdc_args);
            self@Transducer(xdc_args{:}) 
            
            % initialize the TransducerArray 
            for f = string(fieldnames(array_args))'
                self.(f) = array_args.(f);
            end

            % if kerf not set, default it to 30% pitch
            if isempty(self.pitch), self.kerf = 0.3 * self.width; end
        end
    end
    
    % manipulation
    methods
        % scaling
        function self = scale(self, kwargs)
            arguments
                self Transducer
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            args = struct2nvpair(kwargs); % get the arguments as Name/Value pairs
            self = scale@Transducer(self, args{:}); % call superclass method
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [self.pitch] = deal(w*self.pitch);
            end
        end
    end
    
    % define abstract methods
    methods
        function p = positions(self), p = findPositions(self); end
        function [theta, phi, normal, width, height] = orientations(self)
            [theta, phi, normal, width, height] = getOrientations(self);
        end
        function pb = bounds(self), pb = getBounds(self); end
        function pch = patches(self,sub_div), pch = getFieldIIPatches(self,sub_div); end
    end

    
    % define position methods
    methods    
        % get methods
        function p = findPositions(self)
            % returns a 1 x N vector of the positions of the N elements with 0
            % at the center
            array_width = (self.numel - 1) * self.pitch;
            x = linspace(-array_width/2, array_width/2, self.numel);
            p = cat(1, x, zeros(2, numel(x))) + self.offset;
        end
        
        function [theta, phi, normal, width, height] = getOrientations(self)            
            theta = zeros([1, self.numel]);
            phi   = zeros(size(theta));
            ZERO  = zeros(size(theta));
            normal     = [cosd(phi).*sind(theta); sind(phi);  cosd(phi).*cosd(theta)];
            width      = [cosd(theta);            sind(ZERO); -cosd(ZERO).*sind(theta)];
            height     = [sind(phi).*sind(ZERO);  cosd(phi);   sind(phi).*cosd(ZERO)];
        end
        
        function pb = getBounds(self)
            % returns a 3 x 2 matrix of min / max values in x/y/z
            
            % transducer patches of {x,y,z,c} bound tuples
            pch = self.patches([1,1]);
            
            % get min/max bounds of the tx by iterating over each patch
            pb = [inf(3,1), -inf(3,1)];
            for i = 1:self.numel
                pchi = pch{i}(1:3);
                pb(:,1) = min(pb(:,1), cellfun(@(pch) min(pch, [], 'all'), pchi(:)));
                pb(:,2) = max(pb(:,2), cellfun(@(pch) max(pch, [], 'all'), pchi(:)));
            end
        end

        function [p0, dp, N] = getPosDescriptor(self)
            % returns a descriptor of the array's positions in the form
            % {p0, dp, N} where p0 is the first position, dp is the
            % step to the next position, and N is the number of positions
            
            array_width = (self.numel - 1)* self.pitch;
            p0 = [-array_width/2;0;0]; % position
            dp = [self.pitch; 0; 0;]; % sequential change in position
            N = self.numel; % number of positions
        end        
    end
    
    % SIMUS conversion functions
    methods
        function p = getSIMUSParam(self)
            p = struct( ...
                'fc', self.fc, ...
                'pitch', self.pitch, ...
                'width', self.width, ...
                'height', self.height, ...
                'Nelements', self.numel, ...
                'radius', inf, ...
                'bandwidth', 100*self.bw_frac, ... 2-way 6dB fractional bandwidth in % 
                'focus', self.el_focus ... elevation focus
                );
        end
    end

    % Field II conversion function
    methods(Access=public)
        function aperture = getFieldIIAperture(self, focus, element_sub_divisions)
            if nargin < 2 || isempty(focus), focus = [0 0 realmax('single')]; end % ~ 0-deg plane wave
            if nargin < 3, element_sub_divisions = [1,3]; end % no sub divisions
                        
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
        function probe = getUSTBProbe(self)
            probe = uff.linear_array(...
                'N', self.numel, ...
                'pitch', self.pitch, ...
                'element_width', self.width, ...
                'element_height', self.height, ...
                'origin', uff.point('xyz', self.offset(:)') ...
                );
        end
    end

    % Fullwave functions
    methods
        function xdc = getFullwaveTransducer(self, sscan)

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

    % dependent variable methods
    methods
        % get the kerf
        function k = get.kerf(self)
            k = self.pitch - self.width;
        end
        
        % set the kerf
        function set.kerf(self, k)
           self.pitch = self.width + k; 
           if self.pitch < 0
               warning('The interelement spacing is less than 0!');
           end
        end
        
        % get the aperture size
        function a = get.aperture_size(self)
            a = self.numel * self.pitch;            
        end
    end
    
    methods(Static)
        function xdc = L12_3V()
            xdc = TransducerArray(...
                'fc', mean([4e6 11e6]), ...
                'bandwidth', ([4e6 11e6]), ...
                'width', (0.18e-3), ... placeholder @ 90% pitch
                'height', (2e-3), ... placeholder @ 10x picth
                'numel', 192, ...
                'pitch', 0.2e-3, ...
                'focus', 20e-3 ...
                );
        end
        function xdc = L11_5V()
            xdc = TransducerArray(...
                'fc', mean([4.5e6 10e6]), ...
                'bandwidth', ([4.5e6 10e6]), ...
                'width', (0.27e-3), ... placeholder @ 90% pitch
                'height', (3e-3), ... placehoder @ 10x picth
                'numel', 128, ...
                'pitch', 0.3e-3, ...
                'focus', 18e-3 ...
                );
        end
        function xdc = L11_2V()
            xdc = TransducerArray(...
                'fc', 5.1333e+06, ... % Transducer center frequency [Hz]
                'bandwidth', 5.1333e6 + 3e6*[-1 1]/2, ... % bandwidth [Hz]
                'width', 0.270e-3, ... % linear kerf
                'height', 5e-3, ... % Height of element [m]
                'numel', 128, ... % number of elements
                'pitch', 0.300e-3, ... % probe pitch [m]
                'focus', 20e-3 ... % elevation focal depth
                );

        end
        function xdc = L12_5V()
            % Transducer parameters pulled from verasonics L12-5 50mm probe
            xdc = TransducerArray(...
                'fc', 7.5e6, ... % Transducer center frequency [Hz]
                'bandwidth', ([5 11])*1e6, ... % bandwidth [Hz]
                'width', hex2num('3f265251dc6ba641'), ... % element width [m]
                'height', 7.5e-3, ... % Height of element [m]
                'numel', 256, ... % number of elements
                'pitch', hex2num('3f29992e39cf2ea7'), ... % probe pitch [m]
                'focus', 20e-3 ... % elevation focal depth
                );
        end
        function xdc = Verasonics(Trans, c0)
            % xdc = VERASONICS(Trans)
            %
            % Create a Linear Array from the properties defined a
            % Verasonics 'Trans' struct.
            %
            % xdc = VERASONICS(Trans, c0) uses c0 as the sound speed
            % instead of 1540. This is typicaly set by the Verasonics
            % property 'Resource.Parameters.speedOfSound'. Be sure to
            % explicitly set this if other than 1540.

            if nargin < 2, c0 = 1540; end

            % determine the scaling of the properties
            switch Trans.units
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6;
                otherwise
                    error('Conversion from Verasonics trans to TransducerArray not implemented when units not in wavelengths.');
            end

            % set relevant properties
            xdc = TransducerArray(...
                'fc', 1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'bandwidth', 1e6*Trans.Bandwidth([1 end]), ... % bandwidth [Hz]
                'width', scale*Trans.elementWidth, ... % linear kerf
                'height', 1e-3*Trans.elevationApertureMm, ... % Height of element [m]
                'numel', Trans.numelements, ... % number of elements
                'pitch', 1e-3*Trans.spacingMm, ... % probe pitch [m]
                'focus', 1e-3*Trans.elevationFocusMm ... % elevation focal depth
                );
        end
        function xdc = UFF(probe)
            % TRANSDUCERARRAY.UFF - TransducerArray conversion
            %
            % xdc = TRANSDUCERARRAY.UFF(probe) converts the 
            % uff.linear_array probe to a TransducerArray xdc.
            %
            %
            arguments
                probe (1,1) {mustBeA(probe, 'uff.linear_array')}
            end
            xdc = TransducerArray(...
                "pitch", probe.pitch, ...
                "width", probe.element_width, ...
                "height", probe.element_height, ...
                "numel", probe.N_elements, ...
                "offset", - probe.origin.xyz ...
                );
        end

    end
end
