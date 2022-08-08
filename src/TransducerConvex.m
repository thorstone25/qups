
classdef TransducerConvex < Transducer
    
    properties
        % verasonics C5-2V defaults
        radius = 50e-3          % inner radius of curvature (m)
        angular_pitch = 0.5872  % the interelement angular distance (deg)
    end
       
    properties(Dependent)
        pitch                   % interelement spacing along the arc (m)
        kerf                    % the spacing between elements along the arc (m)
        angular_aperture_size   % size of the aperture (deg)
        center                  % center of the circle defining the transducer (m)
    end
    
    % constructor and get/set methods
    methods(Access=public)
        % constructor: accept name/value pairs
        function self = TransducerConvex(varargin)
            
            % setup the transducer args
            if nargin == 1 && isa(varargin{1}, 'struct'), varargin = struct2nvpair(varargin{1}); end

            % initialize the (inherited) Transducer
            self@Transducer(varargin{:}) % assume we don't error on bad inputs
            
            % initialize the TransducerConvex
            if nargin == 1 && isa(varargin{1}, 'uff.curvilinear_array') % (uff object)
                probe = varargin{1};
                props = fieldnames(probe)';
                for p = props
                    switch p{1}
                        case 'radius', self.radius = probe.(p{1});
                    end
                end
                for p = props  % must be set last due to depenedent variable property ordering
                    switch p{1}
                        case 'pitch',  self.pitch  = probe.(p{1}); 
                    end
                end
            else % (name-value)
                for i = 1:2:numel(varargin)
                    switch varargin{i}
                        case 'angular_pitch'
                            self.angular_pitch = (varargin{i+1});
                        case 'radius'
                            self.radius = (varargin{i+1});
                        case 'pitch'
                            self.pitch = varargin{i+1};
                        case 'kerf'
                            if ~isempty(self.width) % we have the width
                                self.kerf = varargin{i+1};
                            elseif ~isempty(self.pitch) % no width: we have pitch
                                kerf = varargin{i+1};
                                self.width = self.pitch - kerf;
                            else % no pitch nor width already:
                                % create an error if no width nor pitch
                                % provided already
                                msgid = 'UltrasoundSimulation:TransducerArray:InvalidArgument';
                                msgtext = 'Either width or pitch must be specified before kerf';
                                ME = MException(msgid, msgtext);
                                throw(ME);
                            end
                    end
                end
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
                [self.radius] = deal(w*self.radius);
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
            % returns a 3 x N vector of the positions of the N elements with
            % the center element at the origin
            
            array_angular_width = (self.numel - 1)* self.angular_pitch;
            theta = linspace(-array_angular_width/2, array_angular_width/2, self.numel);
            z = self.radius * cosd(theta);
            x = self.radius * sind(theta);
            y = zeros(size(theta));
            p = cat(1, x, y, z) + self.center;
        end
        
        function [theta, phi, normal, width, height] = getOrientations(self)
            % Outputs:
            % 
            %   - theta:    the azimuthal angle in cylindrical coordinates
            % 
            %   - phi:      the elevation angle in cylindrical coordinates
            % 
            %   - normal:   a 3 x N array of the element normals
            % 
            %   - width:    a 3 x N array of element width vectors
            % 
            %   - height:   a 3 x N array of element height vectors
            
            array_angular_width = (self.numel - 1)* self.angular_pitch;
            theta = linspace(-array_angular_width/2, array_angular_width/2, self.numel);
            phi   = zeros(size(theta));
            ZERO  = zeros(size(theta));
            normal     = [cosd(phi).*sind(theta); sind(phi); cosd(phi).*cosd(theta)];
            width      = [cosd(theta);           sind(ZERO); -cosd(ZERO).*sind(theta)];
            height     = [sind(phi).*sind(ZERO);  cosd(phi); sind(phi).*cosd(ZERO)];
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
                'radius', self.radius, ...
                'bandwidth', 100*self.bw_frac, ... 2-way 6dB fractional bandwidth in %
                'focus', self.el_focus ... elevation focus
                );
        end
    end


    % Field II conversion function
    methods(Access=public)
        function aperture = getFieldIIAperture(self, focus, element_sub_divisions)
            % creates a transmit / recieve aperture in Field II from the
            % ransducer array
            
            % set defaults
            if nargin < 2 || isempty(focus), focus = [0 0 1e-3]; end % ~ 0-deg plane wave
            if nargin < 3 || isempty(element_sub_divisions), element_sub_divisions = [1,1]; end % compute
            
            % Field II parameters
            xdc_convex_params = { ...
                self.numel,     ... no of elements in x direciton
                self.width,     ... size of elements in x-arc direction
                self.height,    ... size of elements in y direciton
                self.kerf,      ... kerf in x-arc direction
                self.radius,    ... inner convex radius (circumscribed)
                element_sub_divisions(1), ... x sub-divisions
                element_sub_divisions(2), ... y sub-divisions
                focus           ... focal depth
                };
            
            % Generate aperture for emission
            try evalc('field_info'); catch, field_init(-1); end
            aperture = xdc_convex_array (xdc_convex_params{:});
        end
    end
    
    % USTB conversion function
    methods
        function probe = getUSTBProbe(self)
            probe = uff.curvilinear_array(...
                'N', self.numel, ...
                'pitch', self.pitch, ...
                'radius', self.radius, ...
                'element_width', self.width, ...
                'element_height', self.height, ...
                'origin', uff.point('xyz', self.offset(:)') ...
                );
        end
    end

    % Fullwave functions (in-house)
    methods
        function xdc = getFullwaveTransducer(self, sscan)
            
            [dX, dY] = deal(sscan.dx, sscan.dz); % simulation grid step size
            [X0, Y0]= deal(sscan.x(1), sscan.z(1)); % first value
            nX = sscan.size('X' == sscan.order); % grid size
            nY = sscan.size('Z' == sscan.order);
            % map (X, Z) -> (X, Y)
            

            % define variables
            xdc_.npx     = self.numel; % number of elements
            xdc_.dTheta  = self.angular_pitch; % atan2(xdc.ptch,xdc.rad); % angular pitch (radians)
            xdc_.thetas  = self.orientations(); % xdc.dTheta*((-(xdc.npx-1)/2):((xdc.npx-1)/2)); % the thetas defining the transmit elements
            
            % legacy
            % xdc.rad     = self.radius / dY; % 4.957e-2;  % radius of the transducer (pixels)
            % zero_offset = 12.4e-3;      % (deprecated) offset of transducer face, how much it comes into the grid (m)
            % xdc.ptch    = sind(self.angular_pitch) * self.radius / dY; % spatial pitch of the transducer (pixels)
            % xdc.cen     = [(self.offset(1) - X0)/dX, (self.offset(3) - self.radius - Y0)/dY]; % center of the transducer in grid indices
            
            %% Make incoords and outcoords curves

            % define the thetas at the center of each element
            % evenly space, centered at 0deg
            % for n=1:xdc.npx, xdc.thetas(n)=n*xdc.dTheta; end
            % xdc.thetas = xdc.thetas-mean(xdc.thetas);

            % get x-axis and y-axis
            x = X0 + (0:nX-1) * dX; % 1 x X
            y = Y0 + (0:nY-1) * dY; % 1 x Y

            % Make a circle that defines the transducer surface
            inmap = (hypot(x' - self.center(1), y - self.center(3))) < self.radius; % X x Y
            outmap = zeros(nX,nY);
            % inmap(hypot(x+self.offset(1),y+self.offset(3)-self.radius) < self.radius) = true;
            % inmap(circleIdx(size(inmap),xdc.cen,xdc.rad/dY)) = 1;

            % Grab the coords on edge of the circle - larger circle for outcoords
            for i=1:nX

                % find inmap coords
                j = find(inmap(i,:)==0);
                j = j(1);
                inmap(i,1:max([j-8 0]))=0; % make a depth of 8-1=7 pixels in y

                % find outmap coords
                j = find(inmap(i,:)==1);
                if(~isempty(j))
                    j = j(end)+2; % offset by 2 indices in y - this is important!
                    outmap(i,j)=1; % make a depth of 1 pixel in y
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
            % xdc.thetas_in  = atan2(xdc.incoords(:,1)-xdc.cen(1),xdc.incoords(:,2)-xdc.cen(2));
            % xdc.thetas_out = atan2(xdc.outcoords(:,1)-xdc.cen(1),xdc.outcoords(:,2)-xdc.cen(2));

            xdc_.thetas_in   = atan2d(x(1+xdc_.incoords (:,1)) - self.center(1),y(1+xdc_.incoords (:,2))-self.center(3));
            xdc_.thetas_out  = atan2d(x(1+xdc_.outcoords(:,1)) - self.center(1),y(1+xdc_.outcoords(:,2))-self.center(3));

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
                idtheta = abs(xdc_.thetas_in - xdc_.thetas(tt)) < xdc_.dTheta/2;
                xdc_.incoords(idtheta,4) = tt;

                % find center of tt tx element - do each dim separate cause sometimes idtheta is just one value
                % xdc.incoords2(tt,1) = mean(xdc.incoords(idtheta,1));
                % xdc.incoords2(tt,2) = mean(xdc.incoords(idtheta,2));
                xdc_.incoords2(tt,1:2) = mean(xdc_.incoords(idtheta,1:2),1);

                % find which outcoords are assigned to tt
                % less_than_max    = xdc.thetas_out < (xdc.thetas(tt) + xdc.dTheta/2);
                % greater_than_min = xdc.thetas_out > (xdc.thetas(tt) - xdc.dTheta/2);
                % idtheta = find( less_than_max & greater_than_min);
                idtheta = abs(xdc_.thetas_out - xdc_.thetas(tt)) < xdc_.dTheta/2;
                xdc_.outcoords(idtheta,4) = tt;

                % find center of tt rx element - do each dim separate cause sometimes
                % xdc.outcoords2(tt,1) = mean(xdc.outcoords(idtheta,1));
                % xdc.outcoords2(tt,2) = mean(xdc.outcoords(idtheta,2));
                xdc_.outcoords2(tt,1:2) = mean(xdc_.outcoords(idtheta,1:2),1);

            end

            % remove coordinates with no matching element - they are unused
            % this was commented out and it worked
            xdc_.incoords (xdc_.incoords (:,4)==0,:) = [];
            xdc_.outcoords(xdc_.outcoords(:,4)==0,:) = [];

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
        % get the pitch
        function p = get.pitch(self)
            p = 2 * self.radius * sind(self.angular_pitch / 2);
        end
        
        % set the pitch
        function set.pitch(self, p)
            self.angular_pitch = 2 * asind(p / 2 / self.radius);
        end
        
        % get the kerf
        function k = get.kerf(self)
            k = self.pitch - self.width;
        end
        
        % set the kerf
        function set.kerf(self, k)
           self.pitch = self.width + k; 
           if self.angular_pitch < 0
               warning('The interelement angular spacing is less than 0!');
           end
        end
        
        % get the aperture size
        function a = get.angular_aperture_size(self)
            a = (self.numel - 1) * self.angular_pitch;            
        end

        % get the center of the transducer
        function p = get.center(self)
            p = - [0; 0; self.radius] + self.offset;
        end

        % set the center of the transducer
        function set.center(self, p)
            self.offset(1:numel(p)) = p + sub([0; 0; self.radius], 1:numel(p));
        end
    end
    
    methods(Static)
        function xdc = C5_2V()
            xdc = TransducerConvex(...
            'fc', 1e6*mean([2.4 5]), ... % in verasonics doc, this is 4
            'bandwidth', 1e6*([2.4 5]), ...
            'width', 0.46e-3, ...
            'height', 13.5e-3, ...
            'numel', 128, ...
            'radius', 49.57e-3, ...
            'angular_pitch', hex2num('3fe2ca22dae81311'), ...
            'focus', mean([50e-3 70e-3]) ...
            );
        end
        function xdc = Verasonics(Trans, c0)
            % xdc = VERASONICS(Trans)
            %
            % Create a Convex Array from the properties defined a
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
                    warning('Conversion from Verasonics Trans to TransducerArray not supported when units not in wavelengths; values may be incorrect.');
            end

            % set relevant properties
            xdc = TransducerArray(...
                'fc', 1e6*Trans.frequency, ... % Transducer center frequency [Hz]
                'bandwidth', 1e6*Trans.Bandwidth([1 end]), ... % bandwidth [Hz]
                'width', scale*Trans.elementWidth, ... % linear kerf
                'height', 1e-3*Trans.elevationApertureMm, ... % Height of element [m]
                'numel', Trans.numelements, ... % number of elements
                'radius', Trans.radiusMm, ... % radius [m]
                'pitch', 1e-3*Trans.spacingMm, ... % probe pitch [m]
                'focus', 1e-3*Trans.elevationFocusMm ... % elevation focal depth
                );
        end
    end    
end
