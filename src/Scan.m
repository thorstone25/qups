% SCAN - Imaging region definition
%
% The Scan class stores definitions for the imaging region. A Scan
% provides a method to return the image pixel coordinates as an ND-array of
% up to 3 dimensions as well as the row vectors for each individual
% dimension. It also provides convenience methods for defining apodization
% array defined for the Scan.
% 
% Scan is an abstract class. To instantiate a Scan, create a
% ScanCartesian or a ScanPolar.
% 
% See also SCANCARTESIAN SCANPOLAR

classdef Scan < matlab.mixin.Copyable
    properties(Abstract)
        % dimension of change for the pixels
        %
        % order = 'ABC' represents that scan.a varies in dimension 1,
        % scan.b varies in dimension 2, etc. 
        %
        % This is relevant for calls to scan.getImagingGrid() which returns
        % the cartesian coordinates each as a 3D-array.
        order   % dimension of change for the pixels
    end
    
    % dependent parameters
    properties(Abstract, Dependent)
        size        % size of each dimension of the generated pixels
        nPix        % total number of pixels
    end
    
    % get/set & constructor
    methods(Abstract)
        % GETUSTBSCAN - Return a USTB/UFF compatible uff.scan object
        %
        % scan = GETUSTBSCAN(self) returns a uff.scan object
        %
        % See also UFF.SCAN
        scan = getUSTBScan(self)
    end
        
    % imaging computations
    methods(Abstract)
        % GETIMAGINGGRID - get the multidimensional grid for imaging
        %
        % [X, Y, Z, sz] = GETIMAGINGGRID(self) returns the
        % multidimensional arrays corresponding to the positions of the
        % imaging grid and the size of the arrays.
        %
        % G = GETIMAGINGGRID(self) returns the X/Y/Z arrays in a 1 x 3
        % cell (as in G = {X, Y, Z}).
        %
        % P = GETIMAGINGGRID(self, 'vector', true) returns an NDarray whos
        % first dimension is a vector in X,Y,Z order and with the image
        % dimensions raised by one i.e. the following are true:
        %
        % * isequal(size(P)   , [3, self.size])
        % * isequal(P(1,:,:,:), shiftdim(X,-1))
        % * isequal(P(2,:,:,:), shiftdim(Y,-1))
        % * isequal(P(3,:,:,:), shiftdim(Z,-1))
        %
        % The dimension of change for each variable is given by the
        % 'order' property of the Scan
        %
        % Outputs:
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


        % SCALE - Scale units
        %
        % scan = SCALE(scan, 'dist', factor) scales the distance of the
        % properties by factor. This can be used to convert from meters to
        % millimeters for example.
        %
        % Example:
        %
        % % Create a scan
        % scan = ScanCartesian('xb', [-2e-3, 2e-3]); % in meters
        %
        % % convert from meters to millimeters
        % scan = scale(scan, 'dist', 1e3); % in millimeters
        % scan.xb
        %
        % 
        scale(self, kwargs)
    end

    % conversion
    methods
        function s = obj2struct(scan)
            arguments, scan Scan {mustBeScalarOrEmpty}; end
            s = struct(scan); % convert self
        end
    end

    % plotting
    methods
        function gif(self, b_im, filename, h, varargin, kwargs)
            % GIF - Write the series of images to a GIF file
            %
            % GIF(self, b_im, filename) writes the series of images b_im
            % to the file filename. b_im must be a numeric array with
            % images in the first two dimensions.
            %
            % GIF(self, b_im, filename, h) updates the image handle h 
            % rather than creating a new image. Use imagesc to create an 
            % image handle. You can then format the figure prior to 
            % calling this function.
            %
            % GIF(..., Name, Value, ...) forwards Name/Value pairs to
            % imwrite.
            %
            % Example:
            % 
            % % create some data
            % scan = ScanCartesian( ...
            %   'x', 1e-3*(-10:0.2:10), ...
            %   'z', 1e-3*(0:0.2:30) ...
            %  );
            % sz = [scan.size,100];
            % b = complex(rand(sz), rand(sz)) - (0.5 + 0.5i);
            % b_im = rad2deg(angle(b)); % display phase in degrees
            % 
            % % display with imagesc
            % figure;
            % h = imagesc(scan, b_im(:,:,1));
            % colormap hsv;
            % colorbar;
            % title('Random phase');
            % 
            % % make a gif
            % gif(scan, b_im, 'noise.gif', h, 'LoopCount', Inf);
            %
            % See also IMAGESC PLOT

            arguments
                self (1,1) Scan
                b_im {mustBeNumeric, mustBeReal}
                filename (1,1) string
                h (1,1) matlab.graphics.primitive.Image = imagesc(self, b_im, 1)
            end
            arguments(Repeating)
                varargin
            end
            arguments
                kwargs.map (:,:) double = 2^8
                kwargs.LoopCount (1,1) double {mustBePositive} = Inf
                kwargs.DelayTime (1,1) double {mustBePositive} = 1/15
            end

            % parse inputs
            for i = 1:2:numel(varargin), kwargs.(varargin{1}) = varargin{i+1}; end

            % number of slices
            M = prod(size(b_im, 3:max(3,ndims(b_im)))); %#ok<CPROPLC> 

            % get image frames
            for m = M:-1:1, h.CData(:) = b_im(:,:,m); fr{m} = getframe(h.Parent.Parent); end

            % get color space for the image
            [~, map] = rgb2ind(fr{1}.cdata, kwargs.map, 'nodither');

            % get all image data
            im = cellfun(@(fr) {rgb2ind(fr.cdata,map,'nodither')}, fr);
            im = cat(4, im{:});

            % forward options to imwrite (except 'map' option)
            badkwargs = {'map'};
            nvkwargs = struct2nvpair(rmfield(kwargs, badkwargs));
            imwrite(im, map, filename, nvkwargs{:});
        end

        function h = imagesc(self, b, varargin, im_args)
            % IMAGESC - Overload of imagesc
            %
            % h = IMAGESC(self, b) displays the data b on the Scan self.
            %
            % h = IMAGESC(self, b, ax) uses the axes ax instead of the 
            % current axes.
            %
            % h = IMAGESC(..., Name, Value) forwards arguments to MATLAB's 
            % built-in IMAGESC function.
            %
            % See also IMAGESC GIF PLOT

            arguments
                self (1,1) Scan
                b {mustBeNumeric}
            end
            arguments(Repeating)
                varargin
            end
            arguments
                im_args.?matlab.graphics.primitive.Image
            end
            
            % find axis varargs
            if numel(varargin) >= 1 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                ax = varargin{1}; varargin(1) = [];
            else
                ax = gca;
            end

            % imagesc arguments
            im_args = struct2nvpair(im_args);

            % make sure the data size and image size are the same
            assert(all(size(b, 1:3) == self.size, "all"), ...
                "The image size does not match the scan.") %#ok<CPROPLC> 

            % get the axis argumentts
            ax_args = arrayfun(@(c) {self.(c)}, lower(self.order)); % get axis in order
            ax_sing = find(self.size == 1, 1, 'last'); % find singleton dimension
            ax_args = ax_args(~ismember(1:3, ax_sing)); % strip the singleton dimension
            ax_args = ax_args([2 1]); % swap x,y axis arguments
 
            % plot
            if isa(self, 'ScanCartesian')
                h = imagesc(ax, ax_args{:}, squeeze(b), varargin{:},  im_args{:});
                xlabel(ax,'Lateral');
                ylabel(ax,'Axial');
                axis(ax, 'image');

            elseif isa (self, 'ScanPolar')
                h = imagesc(ax, ax_args{:}, squeeze(b), varargin{:}, im_args{:});
                xlabel(ax,'Angle (^o)');
                ylabel(ax,'Range');
                axis(ax, 'tight')
            end
        end

        function h = plot(self, varargin, plot_args)
            % PLOT - Overload of plot
            %
            % h = PLOT(self) plots the pixels of the Scan onto the current 
            % axes and returns the plot handle(s) in h.
            %
            % h = PLOT(self, ax) plot on the axes handle ax.
            %
            % h = PLOT(..., Name, Value, ...) passes the following
            % name/value pairs to the built-in plot function
            %
            % See also PLOT IMAGESC GIF

            arguments
                self (1,1) Scan
            end
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
            end

            % extract axis and other non-Name/Value pair arguments
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes')
                hax = varargin{1}; varargin(1) = [];
            else, hax = gca;
            end

            % get the imaging grid
            [X, ~, Z] = getImagingGrid(self);

            % plot the positions with the options given in the inputs
            plot_args = struct2nvpair(plot_args);
            h = plot(hax, X(:), Z(:), varargin{:}, plot_args{:});
        end
    end
end

