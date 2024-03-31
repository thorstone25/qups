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

classdef Scan < matlab.mixin.Copyable & matlab.mixin.Heterogeneous & matlab.mixin.CustomDisplay
    properties(Abstract)
        % ORDER - ordering of the dimension of the imaging grid
        %
        % scan.order = 'ABC' represents that scan.a is aligned with 
        % dimension 1, scan.b with dimension 2, etc. 
        %
        % This is relevant for calls to scan.getImagingGrid() or 
        % scan.positions() which return the 3D Cartesian coordinates 
        % each as a 3 x A x B x C ND-array.
        order (1,3) char  % dimension of change for the pixels
    end
    properties
        origin (3,1) double % Cartesian origin
    end
    
    % dependent parameters
    properties(Dependent)
        size        % size of each dimension of the generated pixels
        nPix        % total number of pixels
    end

    % USTB interop
    methods
        % QUPS2USTB - Convert a Scan to a USTB/UFF compatible uff.scan
        %
        % uscan = QUPS2USTB(scan) returns a uff.scan object.
        %
        % Example:
        % uscan = QUPS2USTB(ScanCartesian());
        %
        % See also UFF.SCAN
        % USTB interface methods
        function uscan = QUPS2USTB(scan)
            uscan = uff.scan('xyz', reshape(scan.positions(),3,[])');
        end

    end

    % UFF constructor
    methods(Static)
        function scan = UFF(uscan)
        % UFF - Construct a Scan from a uff.scan
        %
        % scan = Scan.UFF(uscan) converts the uff.scan uscan to a Scan.
        %
        % See also SCAN.QUPS2USTB
            arguments, uscan uff.scan; end
            switch class(uscan)
                case 'uff.linear_scan',         scan = ScanCartesian.UFF(uscan);
                case 'uff.sector_scan',         scan = ScanPolar.UFF(uscan);
                case 'uff.linear_3D_scan',      scan = ScanGeneric.UFF(uscan);
                case 'uff.linear_scan_rotated', scan = ScanGeneric.UFF(uscan);
                case 'uff.scan',                scan = ScanGeneric.UFF(uscan);
                otherwise,                      scan = ScanGeneric.UFF(uscan);
            end
        end
    end

    % Verasonics import
    methods(Static)
        function scan = Verasonics(PData, scale)
            % VERASONICS - Create a Scan from a Verasonics PData struct
            %
            % scan = Scan.Verasonics(PData) creates a Scan scan from the
            % pixel data struct PData.
            % 
            % scan = Scan.Verasonics(PData, scale) scales the pixels by
            % scale. The default is 1.
            % 
            % Example:
            %
            % % get wavelengths in meters
            % lambda = Resource.Parameters.speedOfSound / (Trans.frequency * 1e6);
            % 
            % % import in meters
            % scan = Scan.Verasonics(PData, lambda); 
            % 
            % See also TRANSDUCER.VERASONICS
            arguments
                PData (1,1) struct
                scale (1,1) double = 1;
            end
            if ~isfield(PData, 'Coord'), PData.Coord = 'rectangular'; end
            switch PData.Coord
                case "rectangular"
                    ax = arrayfun(@(N,dx) {(0:N-1)*dx}, PData.Size([2,3,1]), PData.PDelta([1,2,3])); % axes
                    ax{1} = ax{1} + 0.5 * PData.PDelta(1); % x-axis offset
                    og = num2cell(PData.Origin); 
                    ax = cellfun(@plus, ax, og, 'UniformOutput',false); % offset by origin
                    scan = ScanCartesian('x', ax{1}, 'y', ax{2}, 'z', ax{3});

                case "polar"
                    warning("Unverified import definition.");
                    ax = arrayfun(@(N,dx) {(0:N-1)*dx}, PData.Size, PData.PDelta([2,1,3])); % axes
                    [r, az, y] = deal(ax{:});
                    apex = PData.Origin;
                    az = rad2deg(az - mean(az));
                    scan = ScanPolar('origin', apex, 'r', r, 'a', az, 'y', y);

                case "spherical"
                    warning("Unverified import definition.");
                    ax = arrayfun(@(N,dx) {(0:N-1)*dx}, PData.Size, PData.PDelta); % axes
                    [r, az, el] = deal(ax{:});
                    apex = PData.Origin;
                    az = rad2deg(az - mean(az));
                    el = rad2deg(el - mean(el));
                    scan = ScanSpherical('origin', apex, 'r', r, 'a', az, 'e', el);
            end
            scan = scan.scale('dist', scale);
        end
    end

    % imaging computations
    methods(Abstract)
        % GETIMAGINGGRID - get the multi-dimensional grid for imaging
        %
        % [X, Y, Z] = GETIMAGINGGRID(scan) returns the 3D arrays 
        % corresponding to the positions of the imaging grid and the size 
        % of the arrays.
        %
        % G = GETIMAGINGGRID(scan) returns the 3D X/Y/Z arrays in a 1 x 3
        % cell (as in G = {X, Y, Z}).
        %
        % P = GETIMAGINGGRID(scan, 'vector', true) returns a single 4D 
        % array whos first dimension is a vector in X,Y,Z order and with 
        % the image dimensions raised by one i.e. the following are true:
        %
        % * isequal(size(P)   , [3, scan.size])
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
        %
        % See also: SCAN/POSITIONS
        [X, Y, Z] = getImagingGrid(scan)
    end
    methods

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
        function scan = scale(scan, kwargs)
            arguments
                scan ScanGeneric
                kwargs.dist (1,1) double
            end
            if kwargs.dist ~= 1
                warning("QUPS:Scan:noOverride", ...
                    "You used 'scale' by " + kwargs.dist + " on a " + class(scan) + "!" ...
                    +newline+"It's not very effective ..." ...
                    );
            end
            scan = copy(scan);
        end
    end

    % positions
    methods
        function  p = positions(scan), p = scan.getImagingGrid("vector", true); end
        % POSITIONS - get the multi-dimensional grid positions
        %
        % P = POSITIONS(scan) returns an NDarray whos first dimension is a 
        % vector in X,Y,Z order and with the image dimensions raised by one
        % i.e. the following is true:
        %
        % * isequal(size(P)   , [3, scan.size])
        %
        % The dimension of change for each variable is given by the
        % 'order' property of the Scan.
        %
        % See also: GETIMAGINGGRID
    end

    % heterogeneous support
    methods (Static,Sealed,Access = protected)
        function scan = getDefaultScalarElement()
            scan = ScanCartesian(); % default heterogeneous instance
        end
    end

    % object display (must be Sealed)
    methods (Sealed, Access = protected, Hidden)
        function groups = getPropertyGroups(scan)
            if ~isscalar(scan)
                groups = getPropertyGroups@matlab.mixin.CustomDisplay(scan);
            else % scan is scalar - show info in groups
                axs = arrayfun(@string, lower(scan.order)); % axes (in order)
                j = scan.size > 1; % filter to non-singular dimensions
                ax = matlab.mixin.util.PropertyGroup(["order","size","origin",axs], "Axes");
                % sz = matlab.mixin.util.PropertyGroup([ "nPix", "n"+axs], "Shape");
                rs = matlab.mixin.util.PropertyGroup(["d"+axs(j)], "Resolution");
                bd = matlab.mixin.util.PropertyGroup([axs(j)+"b"], "Bounds (min/max)");
                lb = matlab.mixin.util.PropertyGroup([axs(j)+"label"],"Labels");
                groups = [ax,rs,bd,lb];

                other_props = setdiff(properties(scan), [groups.PropertyList, "nPix"]);
                other_props = other_props(~contains(other_props, ("d" + axs | axs + ("b" | "label"))));
                if ~isempty(other_props)
                    groups(end+1) = matlab.mixin.util.PropertyGroup(other_props, "Other");
                end
            end
        end

        % Heterogenous display support functions - must be sealed manually
        function header = getHeader(obj)
            header = getHeader@matlab.mixin.CustomDisplay(obj);
        end
        function footer = getFooter(obj)
            footer = getFooter@matlab.mixin.CustomDisplay(obj);
        end
        function displayNonScalarObject(obj)
            displayNonScalarObject@matlab.mixin.CustomDisplay(obj);
        end
        % Do not override this method: a 'Scan' is Abstract and therefore
        % cannot be instanstiated.
        function displayScalarObject(obj)
            displayScalarObject@matlab.mixin.CustomDisplay(obj);
        end
        function displayEmptyObject(obj)
            displayEmptyObject@matlab.mixin.CustomDisplay(obj);
        end
        function displayScalarHandleToDeletedObject(obj)
            displayScalarHandleToDeletedObject@matlab.mixin.CustomDisplay(obj);
        end
    end

    % conversion
    methods
        function s = obj2struct(scan)
            % OBJ2STRUCT - Convert a QUPS object into a native MATLAB struct
            %
            % scan = OBJ2STRUCT(scan) converts the Scan scan and all of 
            % it's properties into native MATLAB structs.
            %
            % Example:
            %
            % % Create a Scan
            % scan = ScanCartesian()
            %
            % % convert to a MATLAB struct
            % scan = obj2struct(scan)
            %
            arguments, scan Scan {mustBeScalarOrEmpty}; end
            W = warning('off', "MATLAB:structOnObject"); % squash warnings
            s = struct(scan); % convert scan
            s.class = class(scan); % append class info
            warning(W); % restore warnings
        end
    end

    % plotting
    methods
        function gif(scan, b_im, filename, h, varargin, kwargs)
            % GIF - Write the series of images to a GIF file
            %
            % GIF(scan, b_im, filename) writes the series of images b_im
            % to the file filename. b_im must be a numeric array with
            % images in the first two dimensions.
            %
            % GIF(scan, b_im, filename, h) updates the image handle h 
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
                scan (1,1) Scan
                b_im {mustBeNumeric, mustBeReal}
                filename (1,1) string
                h (1,1) matlab.graphics.primitive.Image = imagesc(scan, b_im, 1)
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

        function h = imagesc(scan, b, varargin, im_args, kwargs)
            % IMAGESC - Overload of imagesc
            %
            % h = IMAGESC(scan, b) displays the data b on the Scan scan.
            %
            % h = IMAGESC(scan, b, ax) uses the axes ax instead of the 
            % current axes.
            %
            % h = IMAGESC(..., Name, Value) forwards arguments to MATLAB's 
            % built-in IMAGESC function.
            %
            % See also IMAGESC GIF PLOT VOL3D

            arguments
                scan (1,1) Scan
                b {mustBeNumericOrLogical}
            end
            arguments(Repeating)
                varargin
            end
            arguments
                im_args.?matlab.graphics.primitive.Image
                kwargs.slice (1,1) char
                kwargs.index (1,1) {mustBeInteger, mustBeNonnegative} = 1
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
            assert(all(size(b, 1:3) == scan.size, "all"), ...
                "The image size does not match the scan.") %#ok<CPROPLC> 

            % identify the dimension to slice
            ind     = 1:3;
            if isfield(kwargs, 'slice')
                sdim = find(lower(kwargs.slice) == lower(scan.order)); % slicing index
            else
                sdim = find(scan.size == 1, 1, 'last'); % slicing index
            end
            
            % if no empty dim identified, take median of 3rd dim by default
            if isempty(sdim), sdim = 3; kwargs.index = ceil(size(b,sdim)/2); end %#ok<CPROPLC> 
            
            % identify the axes
            ind(sdim) = []; % delete this index
            ax_sym  = lower(scan.order(ind)); % axis symbol
            ax_args = arrayfun(@(c) {scan.(c)}, ax_sym); % get axis in order
            ax_args = ax_args([2 1]); % swap x,y axis arguments ( for imagesc )

            % plot
            if ismatrix(b)
                bpl = b;
            else
                bpl = squeeze(sub(b, kwargs.index, sdim)); % may have higher dimensions
                bpl = bpl(:,:,1);
            end
            if ~isreal(b), bpl = mod2db(bpl); end % convert complex to modulus in dB implicitly
            h = imagesc(ax, ax_args{:}, bpl, varargin{:},  im_args{:}); % plot

            % axes labels
            xlabel(ax,scan.(ax_sym(2) + "label")); % e.g. 'scan.xlabel'
            ylabel(ax,scan.(ax_sym(1) + "label")); % e.g. 'scan.zlabel'

            % default axis setting
            if isa(scan, 'ScanCartesian')
                axis(ax, 'image');
            elseif any(arrayfun(@(s)isa(scan,s), ["ScanPolar", "ScanGeneric", "ScanSpherical"]))
                axis(ax, 'tight')
            else
                warning('QUPS:Scan:UnrecognizedScan', "Unrecognized Scan of class " + class(scan) + ".");
                axis(ax, 'tight')
            end
        end

        function h = plot(scan, varargin, plot_args, kwargs)
            % PLOT - Overload of plot
            %
            % h = PLOT(scan) plots the pixels of the Scan onto the current 
            % axes and returns the plot handle(s) in h.
            %
            % h = PLOT(scan, ax) plot on the axes handle ax.
            %
            % h = PLOT(..., Name, Value, ...) passes the following
            % name/value pairs to the built-in plot function
            %
            % See also PLOT IMAGESC GIF

            arguments
                scan (1,1) Scan
            end
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
                plot_args.DisplayName = 'Grid'
                kwargs.slice (1,1) char
                kwargs.index (1,1) {mustBeInteger, mustBeNonnegative} = 1
            end

            % extract axis and other non-Name/Value pair arguments
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes')
                hax = varargin{1}; varargin(1) = [];
            else, hax = gca;
            end

            % identify the dimension to slice
            if isfield(kwargs, 'slice')
                sdim = find(lower(kwargs.slice) == lower(scan.order)); % slicing index
            else
                sdim = find(scan.size == 1, 1, 'last'); % slicing index
            end

            % if no empty dim identified, take median of 3rd dim by default
            if isempty(sdim), sdim = 3; kwargs.index = ceil(scan.size(sdim)/2); end

            % get the imaging grid
            pts = cell(1,3);
            [pts{:}] = getImagingGrid(scan);

            % reduce to plotting plane
            pts = cellfun(@(axs) {sub(axs, kwargs.index, sdim)}, pts); % slice at slice index
            pts = cellfun(@(x) reshape(x,[],1), pts, 'UniformOutput', false); % vectorize

            % plot the positions with the options given in the inputs
            plot_args = struct2nvpair(plot_args);
            h = plot3(hax, pts{[1 3 2]}, varargin{:}, plot_args{:});

            % set a flat view if possible 
            % TODO: check if default view already set instead
            % TODO: complete for other dims
            if all(diff(h.ZData,1) < 1e-12), view(2); end

            % axes labels
            ax_sym  = lower(scan.order); % axis symbol
            ylabel(hax,scan.(ax_sym(1) + "label")); % e.g. 'scan.zlabel'
            xlabel(hax,scan.(ax_sym(2) + "label")); % e.g. 'scan.xlabel'
            zlabel(hax,scan.(ax_sym(3) + "label")); % e.g. 'scan.ylabel'

        end

        function h = vol3d(scan, b, varargin, im_args)
            % VOL3D - Overload of vol3d
            %
            % h = VOL3D(scan, b) displays the data b on the Scan scan.
            %
            % h = VOL3D(scan, b, ax) uses the axes ax instead of the 
            % current axes.
            %
            % h = VOL3D(..., Name, Value) forwards arguments to vol3d
            %
            % Note: vol3d.m is a MATLAB file exchange file. It can be
            % downloaded <a href="matlab:web('https://www.mathworks.com/matlabcentral/fileexchange/22940-vol3d-v2')">here</a>.
            % 
            % See also VOL3D IMAGESC GIF PLOT

            arguments
                scan (1,1) Scan
                b {mustBeNumeric}
            end
            arguments(Repeating)
                varargin
            end
            arguments
                im_args.texture string {mustBeScalarOrEmpty, mustBeMember(im_args.texture, ["2D", "3D"])}
                im_args.Alpha {mustBeInRange(im_args.Alpha, 0, 1)}
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
            assert(all(size(b, 1:3) == scan.size, "all"), ...
                "The image size does not match the scan.") %#ok<CPROPLC> 

            % identify the axes
            ax_sym  = lower(scan.order); % axis symbol
            ax_args = arrayfun(@(c) {scan.(c+"b")}, ax_sym); % get axis in order
            ax_args = ax_args([2 1 3]); % swap x,y axis arguments ( for vol3d )
            
            % plot
            bpl = b(:,:,:,1);
            if ~isreal(b), bpl = mod2db(bpl); end % convert complex to modulus in dB implicitly
            d_args = [cellstr(["C","X","Y","Z"]+"Data"); [{bpl}, ax_args]];
            h = vol3d('Parent', ax, d_args{:}, varargin{:},  im_args{:}); % plot

            % axes labels
            xlabel(ax,scan.(ax_sym(2) + "label")); % e.g. 'scan.xlabel'
            ylabel(ax,scan.(ax_sym(1) + "label")); % e.g. 'scan.zlabel'
            zlabel(ax,scan.(ax_sym(3) + "label")); % e.g. 'scan.ylabel'

            % default axis setting
            if isa(scan, 'ScanCartesian')
                axis(ax, 'image');
            elseif isa (scan, 'ScanPolar') || isa(scan, 'ScanGeneric')
                axis(ax, 'tight')
            else
                warning('QUPS:Scan:UnrecognizedScan', "Unrecognized Scan of class " + class(scan) + ".");
                axis(ax, 'tight')
            end
        end
    end

    % Dependent properties
    methods
        function sz = get.size(scan)
            sz = arrayfun(@(c) scan.("n"+c), lower(scan.order));
        end
        function set.size(scan, sz)
            sz(numel(sz)+1:3) = 1; % send omitted dimensions to size 1
            for i = 1:3, scan.("n"+lower(scan.order(i))) = sz(i); end
        end
        function n = get.nPix(scan), n = prod(scan.size); end
    end
end

