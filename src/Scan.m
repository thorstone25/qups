classdef Scan < handle
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

