classdef Scan < matlab.mixin.Copyable
    properties(Abstract)
        order
    end
    
    % dependent parameters
    properties(Abstract, Dependent)
        size
        nPix
    end
    
    % get/set & constructor
    methods(Abstract)
        
        setImagingGrid(self, xb, yb, zb)
        
        % GETUSTBSCAN - Return a USTB/UFF compatible uff.scan object
        %
        % scan = GETUSTBSCAN(self) returns a uff.scan object
        %
        % See also UFF.SCAN
        scan = getUSTBScan(self)
        
    end
        
    % imaging computations
    methods(Abstract)
        % SCAN/GETIMAGINGGRID - get the multidimensional grid for imaging
        %
        % [X, Y, Z, sz] = GETIMAGINGGRID(self) returns the
        % multidimensional arrays corresponding to the positions of the
        % imaging grid and the size of the arrays
        %
        % The dimensions of change for each variable is given by the
        % 'order' property of the Scan
        %
        % outputs:
        %   - X:    x coordinate (m)
        %   - Y:    y coordinate (m)
        %   - Z:    z coordinate (m)
        %   - sz:   size of the X, Y, and Z multidimensional arrays
        [X, Y, Z, sz] = getImagingGrid(self)
        
        % SETIMAGEGRIDONTARGET - Set the image grid for a target
        %
        % SETIMAGEGRIDONTARGET(self, target) sets the image grid based on 
        % the boundaries of target.
        %
        % SETIMAGEGRIDONTARGET(self, target, margin) additionally adds a
        % margin around the boundaries of target.
        % 
        % This method leaves unchanged the number of points on the grid, so
        % the resolution may change whenever this function is called.
        %
        % Inputs:
        %   target:     Target object
        %   margin:     a 3 x 2 matrix of x/y/z  min/max bounds for the
        %               imaging grid
        setImageGridOnTarget(self, target, margin)

    end


    % apodization functions
    methods
        function apod = scanlineApodization(self, seq, tol)
            % SCANLINEAPODIZATION Create scanline apodization array
            %
            % apod = SCANLINEAPODIZATION(self, seq) creates an ND-array
            % to mask delayed data using the transmit sequence seq in order
            % to form an image using scanlines.
            % 
            % Scanline apodization is determined by accepting only
            % transmits and scanlines that are aligned across the transmit 
            % aperture, where a scan line is a series of points where only
            % the range varies.
            %
            % The output apod has dimensions I1 x I2 x I3 x N x M where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements, and M is the number of transmits.
            %
            % apod = SCANLINEAPODIZATION(self, seq, tol) uses a tolerance
            % of tol to decided whether to accept the transmit.
            %
            % See also: SCAN/MULTILINEAPODIZATION


            % soft validate the transmit sequence type: it should be focused
            if seq.type ~= "VS", warning(...
                    "Expected sequence type to be VS but instead got " + seq.type + ": This may produce unexpected results."...
                    );
            end

            % create a mask such that the transmit and pixel lateral
            % positions must 'match'
            if nargin < 3, tol = min(diff(sub(seq.focus,1,1))); end % numerical tolerance
            if isa(self, 'ScanCartesian')
                xdim = find(self.order == 'X'); % dimension of change in lateral
                xi = shiftdim(self.x(:), xdim-1); % lateral per pixel
                xv = swapdim(sub(seq.focus,1,1), 2, 5); % 1 x 1 x 1 x 1 x M
            elseif isa(self, 'ScanPolar')
                xdim = find(self.order == 'A'); % dimension of change in azimuth
                xi = shiftdim(self.a(:), xdim-1); % angle per pixel
                xv = swapdim(seq.angles, 2, 5); % 1 x 1 x 1 x 1 x M
            end
            apod = abs(xi - xv) < tol; % create mask
        end

        function apod = multilineApodization(self, seq)
            % MULTILINEAPODIZATION Create multi-line apodization array
            %
            % apod = MULTILINEAPODIZATION(self, seq) creates an ND-array
            % to mask delayed data using the transmit sequence seq in order
            % to form an image using scanlines.
            % 
            % Multilne apodization is determined by linearly weighing
            % scan lines by the transmits that straddle each scan line.
            % Scan lines outside of the minimum and maximum transmit are 
            % weighted by zero.  
            %
            % The output apod has dimensions I1 x I2 x I3 x N x M where
            % I1 x I2 x I3 are the dimensions of the scan, N is the number
            % of receive elements, and M is the number of transmits.
            %
            % See also SCAN/SCANLINEAPODIZATION


            % soft validate the transmit sequence type: it should be focused
            if seq.type ~= "VS", warning(...
                    "Expected sequence type to be VS but instead got " + seq.type + ": This may produce unexpected results."...
                    );
            end

            % extract lateral or angle of transmit in order to compare
            if isa(self, 'ScanCartesian')
                xdim = find(self.order == 'X'); % dimension of change in lateral
                xi = shiftdim(self.x(:), xdim-1); % lateral per pixel
                xv = swapdim(sub(seq.focus,1,1), 2, 5); % 1 x 1 x 1 x 1 x M
            elseif isa(self, 'ScanPolar')
                xdim = find(self.order == 'A'); % dimension of change in azimuth
                xi = shiftdim(self.a(:), xdim-1); % angle per pixel
                xv = swapdim(seq.angles, 2, 5); % 1 x 1 x 1 x 1 x M
            end
            X = self.size(xdim);

            % TODO: switch this to accept multiple transmits instead of
            % just left/right transmit

            % get the apodization based on interpolation weights
            % between each transmit angle
            da = xi - xv; % difference in lateral % Rx x Tx
            lind = da >= 0; % tx left  of xi
            rind = da <= 0; % tx right of xi

            % choose right-most left tx and left-most right tx
            lind = cellfun(@(i) find(i,1,'last' ), num2cell(lind,5), 'UniformOutput', false); % X x 1
            rind = cellfun(@(i) find(i,1,'first'), num2cell(rind,5), 'UniformOutput', false); % X x 1
            val = ~cellfun(@isempty, lind) & ~cellfun(@isempty, rind); % X x 1
            [lind, rind] = deal(cell2mat(lind(val)), cell2mat(rind(val))); % X' x 1

            % set the apodization values for these transmits
            da_lr = swapdim(abs(xv(lind) - xv(rind)), xdim,5); % difference in transmit angle (0 for identical left/right)
            a_l = 1 - (abs(swapdim(xv(lind),xdim,5) - xi(val)) ./ da_lr); % left side apodization
            a_r = 1 - (abs(swapdim(xv(rind),xdim,5) - xi(val)) ./ da_lr); % right side apodization
            ind0 = da_lr == 0; % if no difference between left/right transmits ...
            [a_l(ind0), a_r(ind0)] = deal(1, 0); %  set left to 1, right to 0

            % build Tx x Rx apodization matrix
            apod  = zeros([X, seq.numPulse]); % build in reduced dimensions
            alind = sub2ind(size(apod), find(val), lind); %#ok<CPROPLC> % left matrix indices
            arind = sub2ind(size(apod), find(val), rind); %#ok<CPROPLC> % right matrix indices
            apod(alind) = apod(alind) + a_l; % add left apod
            apod(arind) = apod(arind) + a_r; % add right apod
            apod = ipermute(apod, [xdim, 5, setdiff(1:5, [xdim,5])]); % send X, M to dimensions xdim, 5
        end
    
        function apod = translatingApertureApodization(self, seq, rx, tol)
            % no halp

            % extract lateral or angle of transmit in order to compare
            if isa(self, 'ScanCartesian')
                xdim = find(self.order == 'X'); % dimension of change in lateral
                xi = shiftdim(self.x(:), xdim-1); % lateral per pixel
                xv = swapdim(sub(seq.focus,1,1), 2, 5); % lateral per transmit 1 x 1 x 1 x 1 x M
                % soft error on transducer type
                if isa(rx, 'TransducerArray'), else, warning( ...
                        "Expected a TransducerArray but instead got " + class(rx) + ": This may produce unexpected results."...
                        ); end %#ok<ALIGN> 
                xn = swapdim(sub(rx.positions,1,1), 2,4); % lateral per receiver 
            elseif isa(self, 'ScanPolar')
                xdim = find(self.order == 'A'); % dimension of change in azimuth
                xi = shiftdim(self.a(:), xdim-1); % angle per pixel
                xv = swapdim(seq.angles, 2, 5); % angle per transmit 1 x 1 x 1 x 1 x M
                % soft error on transducer type
                if isa(rx, 'TransducerConvex'), else, warning( ...
                        "Expected a TransducerArray but instead got " + class(rx) + ": This may produce unexpected results."...
                        ); end %#ok<ALIGN> 
                xn = swapdim(rx.orientations,2,4); % lateral per receiver 
            end
            apod = abs(xi - xv) <= tol(1)/2 & abs(xi - xn) <= tol(end)/2; % create mask

        end

        function apod = apertureGrowthApodization(self, seq, rx, f, Dmax)
            % no halp

            % defaults
            if nargin < 5, Dmax = Inf; end
            if nargin < 4, f = 1.5; end

            % soft validate the transmit sequence type: it should be focused
            if seq.type ~= "VS", warning(...
                    "Expected sequence type to be VS but instead got " + seq.type + ". This may produce unexpected results."...
                    );
            end

            % TODO: generalize to polar
            % get the transmit foci lateral position, M in dim 6
            Pv = swapdim(sub(seq.focus,1:3,1), 2, 6); 

            % get the receiver lateral positions, N in dim 5
            Pn = swapdim(sub(rx.positions,1:3,1), 2, 5); % (3 x 1 x 1 x 1 x N)

            % get the points as a vairance in depth
            zdim = find(self.order == 'Z');
            Zi = shiftdim(self.z(:), zdim-1+1); % (1 x I1 x I2 x I3)
            % Zi = [0;0;1] .* Zi; % (3 x I1 x I2 x I3)

            % get the equivalent aperture width (one-sided) and pixel depth
            d = sub(Pn - Pv, 1, 1); % one-sided width where 0 is aligned with transmit
            z = Zi - 0; % depth w.r.t. beam origin
            apod = z > f * abs(2*d) ; % restrict the f-number
            apod = apod .* (abs(2*d) < Dmax); % also restrict the maximum aperture size

            % shift to (I1 x I2 x I3 x N x M)
            apod = shiftdim(apod, 1);
        end

        function apod = acceptanceAngleApodization(self, seq, rx, theta)
            % no halp

            % defaults
            if nargin < 4, theta = 45; end

            % soft validate the transmit sequence type: it should be focused
            if seq.type ~= "VS", warning(...
                    "Expected sequence type to be VS but instead got " + seq.type + ". This may produce unexpected results."...
                    );
            end

            % TODO: generalize to polar
            % get the receiver positions and orientations, N in dim 5
            Pn2 = swapdim(rx.positions   , 2, 4); % (3 x 1 x 1 x N)
            thn = swapdim(rx.orientations, 2, 4); % (3 x 1 x 1 x N)

            % get the points as a variance in depth and lateral
            % xdim = find(self.order == 'X');
            % zdim = find(self.order == 'Z');
            % Xi = shiftdim(self.x(:), xdim-1+1); % (I1 x I2 x I3)
            % Zi = shiftdim(self.z(:), zdim-1+1); % (I1 x I2 x I3)
            [Xi, ~, Zi] = self.getImagingGrid(); % (I1 x I2 x I3)

            % restrict to points where the angel is less than theta at the
            % receiver
            thi = atan2d(Xi - sub(Pn2,1,1), Zi - sub(Pn2,3,1)); % angle at which the ray hits the element
            apod = abs(thi - thn) <= theta; % (I1 x I2 x I3 x N x M)
        end

    end


    methods
        function h = imagesc(self, b, varargin)
            % SCAN/IMAGESC - Overload of imagesc
            %
            % h = IMAGESC(b) displays the data b on the Scan.
            %
            % h = IMAGESC(b, ax) uses the axes ax instead of the current
            % axes.
            %
            % h = IMAGESC(..., Name, Value, ...) forwards arguments to
            % MATLAB's IMAGESC function.
            %
            % See also IMAGESC

            if nargin >= 3 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                ax = varargin{1}; varargin(1) = [];
            else
                ax = gca;
            end
            if isa(self, 'ScanCartesian')
                % assumes that data is in 'ZXY' order
                if self.order ~= "ZXY", warning('Scan axes are incorrect - this function may give unexpected results.'); end
                h = imagesc(ax, self.x, self.z, b, varargin{:});
                xlabel(ax,'Lateral (m)');
                ylabel(ax,'Axial (m)');
                axis(ax, 'image');

            elseif isa (self, 'ScanPolar')
                % assumes that data is in 'RAY' order
                if self.order ~= "RAY", warning('Scan axes are incorrect - this function may give unexpected results.'); end
                h = imagesc(ax, self.a, self.r, b, varargin{:});
                xlabel(ax,'Angle (^o)');
                ylabel(ax,'Range (m)');
                axis(ax, 'tight')
            end
        end

        function h = plot(self, varargin)
            % determine the axis
            if nargin >= 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                hax = varargin{1}; varargin(1) = []; % delete this input
            else
                hax = gca;
            end

            % get the imaging grid
            [X, ~, Z] = getImagingGrid(self);

            % plot the positions with the options given in the inputs
            h = plot(hax, X(:), Z(:), varargin{:});
        end
    end
end

