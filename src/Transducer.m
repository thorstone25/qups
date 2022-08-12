% TRANSDUCER Abstract Transducer class
% 
% Superclass for a medical ultrasound transducer. Any TRANSDUCER has
% methods for the positions and orientations of it's elements as well as
% the characteristics of the elements (width/height, impulse response, 
% central frequency, etc.). This class offers definitions for common
% transducers and conversion functions between real transducers (mainly
% from Verasonics) and simulation programs (k-Wave, Fullwave, FieldII).
% 
% See also TRANSDUCERARRAY TRANSDUCERCONVEX TRANSDUCERPISTON

classdef (Abstract) Transducer < matlab.mixin.Copyable
    properties
        fc = 5e6        % center frequency
        bw = [3.5e6 6.5e6] % bandwidth
        width = 1.19e-4 % width of an element (m)
        height = 6e-3   % height of an element (m)
        numel = 128     % number of elements
        offset = [0;0;0]% the offset from the origin compared to fieldII's definitions
        impulse = repmat(Waveform(),[0,0]) % the impulse response function of the element
    end
    
    properties
        el_focus = realmax('single') % elevation focal depth
    end
    
    properties(GetAccess=public, SetAccess=protected, Dependent)
        area            % area of an element (m^2)
    end
    
    properties(Dependent)
        bw_frac         % fractional bandwidth
    end
    properties(Dependent, Hidden)
        origin
    end
    
    % constructor
    methods
        function self = Transducer(varargin)
            if nargin == 1 && isa(varargin{1}, 'Transducer') % Transducer initialization object
                % set all non-dependent properties of the object
                xdc = varargin{1};
                mc = metaclass(xdc);
                props = mc.PropertyList;
                for p = props(~[props.Dependent] & ~[props.NonCopyable])'
                    self.(p.Name) = xdc.(p.Name);
                end

            elseif nargin == 1 && isa(varargin{1}, 'uff.probe') % uff Constructor
                probe = varargin{1};
                props = fieldnames(probe)';
                for p = props
                    switch p{1}
                        case 'N',  self.numel  = probe.(p{1});
                        case 'N_elements',  self.numel  = probe.(p{1});
                        case 'element_width',       self.width  = probe.(p{1});
                        case 'element_height',      self.height = probe.(p{1});
                        case 'origin',      self.origin(:) = probe.origin.xyz;
                    end
                end
            else % name-value pair initialization
                for i = 1:2:nargin
                    switch lower(varargin{i})
                        case {'width', 'height', 'fc', 'bw', 'el_focus', 'numel', 'offset', 'impulse'}
                            self.(lower(varargin{i})) = (varargin{i+1});
                        case {'bandwidth'}
                            self.bw = (varargin{i+1});
                        case 'fractional_bandwidth'
                            self.setBandwidthFractional(varargin{i+1});
                        case {'focus'}
                            self.el_focus = varargin{i+1};
                    end
                end
            end

            % regardless of input, if impulse is empty, initialize it
            if isempty(self.impulse)
                self.impulse = self.ultrasoundTransducerImpulse(); 
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
            % % Create a Transducer
            % xdc = TransducerArray() % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % xdc = scale(xdc, 'dist', 1e3, 'time', 1e6); % mm, us, MHz
            % xdc.width
            % xdc.fc
            %

            arguments
                self Transducer
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            self = copy(self);
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [self.width, self.height, self.offset] = deal(w*self.width, w*self.height, w*self.offset);
            end
            if isfield(kwargs, 'time')
                w = kwargs.time;
                % scale time (e.g. s -> us / Hz -> MHz)
                [self.fc, self.bw] = deal(self.fc/w, self.bw/w);
                self.impulse = scale(self.impulse, 'time', w);
            end
        end
    end


    % transducer specific methods
    methods (Abstract)
        % inferred position methods
        p = positions(self); % compute the positions of the elements
        % TRANSDUCER/ORIENTATIONS Orientation of the elements
        %
        % [theta, phi] = ORIENTATIONS(self) returns the azimuth and
        % elevation angles of the elements of the transducer. 
        %
        % See also TRANSDUCER/POSITIONS
        [theta, phi, normal, width, height] = orientations(self); % compute the orientations
        pb = bounds(self); % compute the 3 x 2 array of x/y/z by min/max 
        pch = patches(self,sub_div); % computes a [Ndiv x Nel] cell array of
        % {X,Y,Z,C} tuples of 2x2 matrices specifying four corners of each 
        % patch that is a subdivision of the elements
    end
    
    % toolbox conversion functions
    methods (Abstract)
        aperture = getFieldIIAperture(self, focus, element_sub_divisions); % get the fieldII aperture object
        probe = getUSTBProbe(self); % get the USTB probe object 
    end

    % Verasonics conversion functions
    methods (Abstract, Static)
        % xdc = VERASONICS(Trans)
        %
        % Construct a Transducer from the properties defined a
        % Verasonics 'Trans' struct.
        %
        % xdc = VERASONICS(Trans, c0) uses c0 as the sound speed
        % instead of 1540. This is typicaly set by the Verasonics
        % property 'Resource.Parameters.speedOfSound'. Be sure to
        % explicitly set this if other than 1540.
        xdc = Verasonics(Trans, c0)

    end

    % SIMUS functions
    methods (Abstract)
        p = getSIMUSParam(self)
    end
    
    % kWave functions
    methods
        function [ksensor, ksensor_ind, sens_map] = getKWaveSensor(self, kgrid, kgrid_origin, el_sub_div)
            % TODO: doc this function: it gives you a k-Wave ready sensor
            % structure
            %             arguments
            %                 self (1,1)
            %                 kgrid (1,1) kWaveGrid
            %                 kgrid_origin
            %                 el_sub_div (1,2) double = [1,1];
            %             end
            if nargin < 4, el_sub_div = [1,1]; end
            
            % get the sensor and source masks. This is the hard part: how do I do this
            % for a convex probe on a grid surface?
            
            % cast grid points to single type for efficiency and base grid on the origin
            [gxo, gyo, gzo] = deal(kgrid_origin(1), kgrid_origin(2), kgrid_origin(3)); % grid {dim} origin
            [gxv, gyv, gzv] = deal(single(gxo + kgrid.x_vec), single(gyo + kgrid.y_vec), single(gzo + kgrid.z_vec)); % grid {dim} vector
            [gx0, gy0, gz0] = deal(gxv(1), gyv(1), gzv(1)); % grid {dim} first point
            [gxd, gyd, gzd] = deal(single(kgrid.dx), single(kgrid.dy), single(kgrid.dz)); % grid {dim} delta
            
            % get array sizing - make size '1' and step 'inf' for sliced dimensions
            [Nx, Ny, Nz] = dealfun(@(n) n + (n==0), kgrid.Nx, kgrid.Ny, kgrid.Nz);
            if gxd == 0, gxd = inf; end
            if gyd == 0, gyd = inf; end
            if gzd == 0, gzd = inf; end
            
            % get local element size reference
            [width_, height_] = deal(self.width, self.height);
            
            % get regions for each receive transducer element
            el_cen = self.positions(); % center of each element
            [theta, phi, el_dir, el_wid, el_ht] = self.orientations(); % orientations vectors for each element
            
            % reorder to map between kwave and conventional ultrasound coordinates
            [el_dir, el_cen, el_wid, el_ht] = dealfun(@(v)v([3 1 2],:), el_dir, el_cen, el_wid, el_ht);
            
            % get edges in x/y/z (i.e. patches) for each sub element
            patches = self.patches(el_sub_div)'; % (E x N)
            nSub = size(patches, 1); % number of subelements
            
            % convert to cell array of points with x/y/z in 1st dimension
            % in kgrid order
            p_patches = cellfun(@(pch) vec(cellfun(@(p)mean(p(:)), pch([3 1 2]))), patches, 'UniformOutput', false);
            
            % send grid data to workers if we have a parallel pool
            if(isvalid(gcp('nocreate'))), sendDataToWorkers = @parallel.pool.Constant;
            else, sendDataToWorkers = @(x) struct('Value', x);
            end
            [gxv, gyv, gzv] = dealfun(sendDataToWorkers, gxv, gyv, gzv);
            
            % initialize
            ksensor_all = false([Nx, Ny, Nz]);
            
            % get sensing map for all patches for all elements
            fprintf('\n');
            parfor el = 1:self.numel
                tt = tic;
                % set variables for the element
                [el_dir_, el_cen_, el_wid_, el_ht_, p_patches_] = deal(el_dir(:,el), el_cen(:,el), el_wid(:,el), el_ht(:,el), p_patches(:,el));
                
                % initialize
                ksensor{el} = struct('mask', false([Nx, Ny, Nz]));
                
                for i = nSub:-1:1
                    % get the center of the subelement
                    pcen = p_patches_{i};
                    
                    % get zero crossing indices
                    xind = 1 + (pcen(1) - gx0) / gxd;
                    yind = 1 + (pcen(2) - gy0) / gyd;
                    zind = 1 + (pcen(3) - gz0) / gzd;
                    
                    % get index on both sides
                    [xind, yind, zind] = dealfun(@(n) vec(unique([floor(n); ceil(n)])), xind, yind, zind);
                    
                    % shift to appropriate dimensions
                    xind = shiftdim(xind,  0); % in first dimension
                    yind = shiftdim(yind, -1); % in second dimension
                    zind = shiftdim(zind, -2); % in third dimension
                    
                    % get outer product expansion for all sets of indices
                    [xind, yind, zind] = ndgrid(xind, yind, zind);
                    
                    % vectorize
                    [xind, yind, zind] = dealfun(@vec, xind, yind, zind);
                    
                    % [xind, yind, zind] = msk_fun_neighbors(pcen);
                    ind_msk = sub2ind([Nx, Ny, Nz], xind, yind, zind);
                    
                    % get vector from element to the grid pixels
                    vec_ = gather([gxv.Value(xind), gyv.Value(yind), gzv.Value(zind)]') - pcen; %#ok<PFBNS> % 3 x J;
                    
                    % get plane wave phase shift distance as the inner product
                    % sign is whether in front or behind
                    d = (el_dir_' * vec_);
                    
                    % get subelement apodization accounting for cosine
                    % distribution along the transducer surface
                    % I don't ... actually know how to do this ...
                    a = cosd(90 * 2 * el_wid_' * (pcen - el_cen_) / width_ ) ...
                        * cosd(90 * 2 * el_ht_'  * (pcen - el_cen_) / height_);
                    a = 1; %%% DEBUG %%%
                    
                    % save indices, distances, and normal to the sensitivity map
                    sens_map(i,el) = struct(...
                        'amp', a, ...
                        'dist', vec(d),...
                        'mask_indices', ind_msk,...
                        'element_dir', el_dir_ ...
                        );
                    
                    % set the overall mask for each element
                    ksensor{el}.mask(ind_msk) = true;
                end
                
                % reduce to make full recording sensor mask
                ksensor_all = ksensor_all | ksensor{el}.mask;
                
                % report computation timing
                fprintf('Processed %i subelements for element %i in %0.5f seconds.\n', nSub, el, toc(tt));
            end
            
            % get the translation map for sensor to receiver element
            ind_full_sensor = find(ksensor_all);
            parfor (el = 1:self.numel, 0)
                [~, ksensor_ind{el}] = ismember(find(ksensor{el}.mask), ind_full_sensor);
            end
            
            % ensure results on the host
            ksensor_ind  = cellfun(@gather, ksensor_ind, 'UniformOutput', false);
            ksensor = cellfun(@(k)struct('mask',gather(k.mask)), ksensor, 'UniformOutput', false);
        end
    end

    % FDTD functions
    methods
        function [mask, el_weight, el_dist, el_ind] = elem2grid(self, scan, el_sub_div)
            % ELEM2GRID - Transducer element to grid mapping
            % 
            % [mask, el_ind, el_weight, el_dist] = ELEM2GRID(self, scan)
            % returns a binary mask and a set of indices, weights, and 
            % distances defined for each element give a Transducer self and
            % a ScanCartesian scan. The element indicies are defined on the
            % vectorized non-zero indices of the mask i.e. on the indices 
            % of `find(mask)`. 
            
            % arguments
            %     self (1,1) Transducer
            %     scan (1,1) ScanCartesian
            %     el_sub_div (1,2) double = [1,1]
            % end
            if nargin < 3, el_sub_div = [1,1]; end

            % get the sensor and source masks. This is the hard part: how do I do this
            % for a convex probe on a grid surface?

            vec = @(x)x(:); % define locally for compatibility
            
            % cast grid points to single type for efficiency and base grid on the origin
            [gxv, gyv, gzv] = deal(single(scan.x), single(scan.y), single(scan.z)); % grid {dim} vector
            [gx0, gy0, gz0] = deal(scan.x(1), scan.y(1), scan.z(1)); % grid {dim} first point
            [gxd, gyd, gzd] = deal(single(scan.dx), single(scan.dy), single(scan.dz)); % grid {dim} delta
            pg = scan.getImagingGrid();
            pg = cell2mat(cellfun(@(x){shiftdim(x,-1)}, pg(:))); % 3 x Z x X x Y
            pdims = arrayfun(@(c){find(c == scan.order)}, 'XYZ'); % dimensions mapping
            [xdim, ydim, zdim] = deal(pdims{:});
            iord = arrayfun(@(d)find([pdims{:}] == d), [1,2,3]); % inverse mapping
            
            % get array sizing - size '1' and step 'inf' for sliced dimensions
            [Nx, Ny, Nz] = deal(scan.nx, scan.ny, scan.nz);
            
            % get local element size reference
            [width_, height_, sz_] = deal(self.width, self.height, scan.size);
            
            % get regions for each receive transducer element
            el_cen = self.positions(); % center of each element
            [theta, phi, el_dir, el_wid, el_ht] = self.orientations(); % orientations vectors for each element
            
            % reorder to map between kwave and conventional ultrasound coordinates
            % [el_dir, el_cen, el_wid, el_ht] = dealfun(@(v)v([3 1 2],:), el_dir, el_cen, el_wid, el_ht);
            
            % get edges in x/y/z (i.e. patches) for each sub element
            patches = self.patches(el_sub_div)'; % (E x N)
            nSub = size(patches, 1); % number of subelements
            
            % convert to cell array of points with x/y/z in 1st dimension
            p_patches = cellfun(@(pch) vec(cellfun(@(p)mean(p(:)), pch(1:3))), patches, 'UniformOutput', false);
            
            % send grid data to workers if we have a process pool
            if(~isempty(gcp('nocreate')) && isa(gcp('nocreate'), 'parallel.ProcessPool')), 
                  sendDataToWorkers = @parallel.pool.Constant;
            else, sendDataToWorkers = @(x) struct('Value', x);
            end
            [gxv, gyv, gzv] = dealfun(sendDataToWorkers, gxv, gyv, gzv);
            
            % get sensing map for all patches for all elements
            parfor el = 1:self.numel
                % set variables for the element
                [el_dir_, el_cen_, el_wid_, el_ht_, p_patches_] = deal(el_dir(:,el), el_cen(:,el), el_wid(:,el), el_ht(:,el), p_patches(:,el));
                
                % initialize
                % mask{el} = false(sz_);
                
                for i = nSub:-1:1
                    % get the center of the subelement
                    pcen = p_patches_{i}; % (3 x 1)
                    
                    % get zero crossing in indices
                    xind = 1 + (pcen(1) - gx0) / gxd;
                    yind = 1 + (pcen(2) - gy0) / gyd;
                    zind = 1 + (pcen(3) - gz0) / gzd;
                    
                    % get integer index on both sides of zero crossing
                    [xind, yind, zind] = dealfun(@(n) vec((floor(n) + [0; 1])), xind, yind, zind);

                    % set index of sliced dimensions
                    if isinf(gxd), xind(:) = 1; end
                    if isinf(gyd), yind(:) = 1; end
                    if isinf(gzd), zind(:) = 1; end
                    
                    % shift to appropriate dimensions for the scan
                    pind = {xind, yind, zind}; 
                    pind = pind(iord); % shift to proper order given by the scan 
                    [pind{:}] = ndgrid(pind{:}); % outer product expansion 
                    [pind{:}] = dealfun(vec, pind{:}); % vectorize
                    ind_msk = sub2ind(sz_, pind{:}); % get linear indices of the scan (J x 1)
                    pind_ = [pind{:}]; % combine
                    [xind, yind, zind] = deal(pind_(:,xdim), pind_(:,ydim), pind_(:,zdim)); % split
                    
                    % get vector from element to the grid pixels
                    pgrd  = [gxv.Value(xind'); gyv.Value(yind'); gzv.Value(zind')]; %#ok<PFBNS> all data passed is necessary
                    vec_ = gather(pgrd) - pcen; % 3 x J
                    
                    % get plane wave phase shift distance as the inner product
                    % sign is whether in front or behind
                    d = (el_dir_' * vec_).'; % J x 1
                    
                    % get subelement apodization accounting for cosine
                    % distribution along the transducer surface
                    % I don't ... actually know how to do this ...

                    % a = cosd(90 * 2 * el_wid_' * (pcen - el_cen_) / width_ ) ...
                    %   * cosd(90 * 2 * el_ht_'  * (pcen - el_cen_) / height_);
                    % a = a + false(size(d)); 

                    % get subelement apodization accounting for grid
                    % interpolation
                    a = prod(1 - (abs(pgrd - pcen) ./ [gxd; gyd; gzd]),1)'; % J x 1
                    % a = ones(size(d)); %%% DEBUG %%%
                    
                    % save outputs
                    mask_ind{i,el} = ind_msk; % indices of the scan
                    weights {i,el} = a; % amplitudes
                    dists   {i,el} = d; % distances for phase shifting response
                end
            end
            
            % check that each (sub-)element has some an associated grid index
            assert(~any(cellfun(@isempty, mask_ind), 'all'), 'Unable to assign all (sub-)elements to the grid.');
            
            % reduce to make full recording sensor mask
            mask = false(scan.size);
            mask(unique(cat(1, mask_ind{:}))) = true;

            % get the translation map for grid mask index to receiver element
            ind_el  = cellfun(@(m) gather(argn(2, @ismember, m, find(mask))), mask_ind, 'UniformOutput', false);

            % make all into arrays, defined by mask indices
            for el = self.numel:-1:1
                el_weight{el} = cat(1, weights{:,el});
                el_dist  {el} = cat(1, dists  {:,el});
                el_ind   {el} = cat(1, ind_el {:,el});
            end

            % make numeric arrays (J x N)
            el_weight = cat(2, el_weight{:});
            el_dist   = cat(2, el_dist  {:});
            el_ind    = cat(2, el_ind   {:});
        end
    end

    methods %(Abstract)
        % karray = kWaveArray(self, og, varargin)
        function karray = kWaveArray(self, dim, og, varargin)
        % KWAVEARRAY - Create a kWaveArray object from the Transducer
        %
        % karray = KWAVEARRAY(self, dim, og) creates a kWaveArray from the
        % Transducer self defined on a kWaveGrid with dim dimensions and 
        % origin og.
        %
        % karray = KWAVEARRAY(..., Name, Value, ...) forwards Name/Value
        % pair arguments to the kWaveArrray constructor.
        %
        % Inputs:
        %   Axisymmetric
        %   BLITolerance
        %   BLIType
        %   UpsamplingRate
        % 
        % See also KWAVEARRAY

        %{
    end

    methods
        %}
        % function karray = kWaveArray(self, og, varargin)
            
            % TODO: set axisymmetric if it makes sense
            % TODO: allow input argument to specify dimensions
            
            % initialize
            karray = kWaveArray(varargin{:}); 

            % get positions and orientations
            p = self.positions; % TODO: translate the entire array?
            [th, phi] = self.orientations;
            n = self.numel;

            % build the rotation vectors
            % v = [cosd(phi); cosd(phi); sind(phi)] .* [cosd(th); sind(th); 1];
            rot = [phi; th; th] .* [1;1;0];

            % convert to kWaveGrid coorindates: axial x lateral x elevation
            switch dim
                case 2, % we only rotate in x-z plane
                    [p, rot] = deal(p([3,1  ],:), rot(2,:)); 
                    [arr_off, arr_rot] = deal(-og(1:2), 0);
                case 3, 
                    [p, rot] = deal(p([3,1,2],:), rot([3,1,2],:)); % I think this is right?
                    [arr_off, arr_rot] = deal(-og(1:3), [0;0;0]);
                otherwise, error("Wrong number of dimensions (" + dim + "). An array requires 2 or 3 dimensions.");
            end

            % add each element to the array
            for i = 1:n
                karray.addRectElement(p(:,i), self.width, self.height, rot(:,i));
            end

            karray.setArrayPosition(arr_off, arr_rot);
        end
    end

    % fullwave functions
    methods (Abstract)
        % GETFULLWAVETRANSDUCER - define a fullwave transducer structure
        %
        % xdc = GETFULLWAVETRANSDUCER(self, scan) creates a fullwave 
        % compatible structure xdc defined on the ScanCartesian scan.
        %
        % Outputs:
        %   - xdc: a Fullwave compatible struct with (at least) the 
        %           following fields:
        %       - npx           number of elements
        %       - inmap         mask of the input pixels
        %       - nInPx         number of input pixels
        %       - nOutPx        number of output pixels
        %       - incoords      (x,y,1,el) coordinate pairs of the input pixels
        %       - outcoords     (x,y,1,el) coordinate pairs of the output pixels
        xdc = getFullwaveTransducer(self, scan) 
    end

    methods
        % get positions
        function p = getPositions(self)
            % returns a 3 x N vector of the positions of the N elements with
            % the center element at the origin
            p = self.positions();
        end

        % get number of elements
        function n = nTx(self), n = self.numel; end
        function n = nRx(self), n = self.numel; end
    end
    
 
    % Field II calls - rely on being able to get a FieldII aperture
    methods(Access=public)
        function p = getFieldIIPositions(self)
            % returns a 3 x N vector of the positions of the N elements as
            % represented in FieldII
            
            ap = self.getFieldIIAperture();
            data = xdc_get(ap, 'rect');
            p = data([24 25 26], :) + self.offset;
            xdc_free(ap);
        end
        
        function el_patches = getFieldIIPatches(self, el_sub_div)
            % GETFIELDIIPATCHES- get cell array of patches
            % 
            % el_patches = GETFIELDIIPATCHES(self, el_sub_div) returns the 
            % corners of the mathematical elements represented in FieldII. 
            % Each element is represented by a {x,y,z,c} tuples (a.k.a. a 
            % "patch") which holds [x1, x2; x3 x4], which are the corners 
            % of the elements
            % 
            % Inputs:
            %   - el_sub_div    2 x 1 vector of integers specifying the 
            %                   number of mathematical elements per
            %                   physical element
            %
            % Outputs:
            %   - el_patches       array of {x,y,z,c} tuples
            %
            % See also GETFIELDIIBARYCENTERS
            
            % parse inputs
            if nargin < 2, el_sub_div = [1 1]; end
            
            % get the aperture
            ap = self.getFieldIIAperture([0 0 0], el_sub_div);
            
            % extract corners
            data = xdc_get(ap, 'rect');
            ind = 1:size(data,2);
            
            % release the aperture
            xdc_free(ap);
            
            % extract the corners into {x,y,z,c} tuples (a.k.a. a patch)
            % where x/y/z/c are each 2 x 2 arrays of points
            corner_map = @(i){...
                self.offset(1) + ...
                [   data(11,i), data(20,i); ...
                    data(14,i), data(17,i)], ...
                self.offset(2) + ...
                [   data(12,i), data(21,i); ...
                    data(15,i), data(18,i)], ...
                self.offset(3) + ...
                [   data(13,i), data(22,i);...
                    data(16,i), data(19,i)], ...
                    data(5, i) * ones(2,2),...
                };
            el_patches = arrayfun(corner_map, ind, 'UniformOutput', false); % array of {x,y,z,c} tuples

            % reshape by element (N x E), E is subelements
            el_patches = reshape(el_patches, prod(el_sub_div), [])';
        end
        
        function p = getFieldIIBaryCenters(self, el_sub_div)
            % computes the positions of the transducer elements by
            % computing the barycenters (geometric mean) of the corners of
            % the mathematical elements computed by FieldII.
            if nargin < 2, el_sub_div = [1,1]; end
            
            el_patches = self.getFieldIIPatches(el_sub_div);
            
            p = cellfun(@(patch) cellfun(@(comp) mean(comp(:)),...
                         patch(1:3)'), el_patches, 'UniformOutput', false);
            p = reshape(cat(1, p{:}), [3, size(p)]);
        end
    end
    
    % internal subroutines
    methods(Access=public)
        function impulse = ultrasoundTransducerImpulse(self, opt)
            % impulse response function for an ultrasound transducer at time t
            % where 0 is defined at the center of the pulse
            % optionally, choose a waveform as one of {'pulse'* | 'delta'}
            
            % defaults
            if nargin < 2, opt = 'pulse'; end
            
            switch opt
                case 'pulse'
                    % get energy strength cutoff time
                    bwr = -6;
                    tpr = -80;
                    
                    tc = gauspuls('cutoff',self.fc, self.bw_frac, bwr, tpr);
                    
                    % get the function
                    impulse_fun = Transducer.cgauspulsfun(self.fc, self.bw_frac, bwr);
                    
                    % make a Waveform object
                    impulse = Waveform('fun', impulse_fun, 't0', -tc, 'tend', tc);
                case 'delta'
                    impulse = Waveform('fun', @(t) complex(ones(size(t))), 't0', 0, 'tend', 0);
            end
        end
    end
    
    % dependent methods
    methods
        function set.bw(self, b),
            if(length(b) ~= 2),
                error("bw must be a length 2 vector of the passband cutoff frequencies."); 
            end, 
            self.bw = b;
        end

        function a = get.area(self), a = self.width * self.height; end
        
        function b = get.bw_frac(self), b = range(self.bw) ./ self.fc; end
        
        function set.bw_frac(self, bf), self.bw = mean(self.fc) + mean(self.fc) * bf * ([-1 1] ./ 2); end
        
        function o = get.origin(self), o = - self.offset; end
        
        function set.origin(self, o), self.offset = -o; end
    end
    
    % plot functions
    methods(Access=public)
        function varargout = patch(self, el_sub_div, varargin)
            % PATCH - Overload of patch for a transducer 
            %
            % PATCH(self) plots the aperture from its patch representation.
            % It overloads the patch function so that arguments valid for 
            % patch are valid here. 
            % 
            % PATCH(self, el_sub_div) specifices element subdivisions. 
            %
            % PATCH(self, el_sub_div, ax) plots on the axes ax.
            %
            % PATCH(self, el_sub_div, ..., Name, Value, ...) passes 
            % name-value pair arguments to patch.
            %
            % hp = PATCH(self, ...) returns a patch object handle
            %
            % [hp, hax] = patch(self, ...) additionally returns an axes
            % handle
            %
            % Inputs
            %   - el_sub_div:    2 x 1 vector of element subdivisions
            %
            % See also PATCH
            
            if nargin < 2, el_sub_div = [1,1]; end

            % parse inputs
            if nargin > 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                axs = varargin{1}; varargin(1) = [];
            else
                axs = gca;
            end

            % get the aperture as a cell array of {x,y,z,c} pairs
            el_patches = self.patches(el_sub_div);
            
            % map point ordering to generate a square with patch; the
            % corners must be specified in sequential order
            % make column vector, drop apodization,
            ind_map = [1 2 4 3]';
            el_patches = cellfun(@(patch) cell2mat(cellfun(@(component) ...
                reshape(component(ind_map),[],1), ... .* 1e3, ...
                patch, 'UniformOutput', false)), ...
                el_patches, 'UniformOutput', false);
            
            % convert to get vertex x index x X/Y/Z array
            el_patches = permute(cat(3, el_patches{:}), [1,3,2]);
            
            % get X/Y/Z/C components
            Xp = el_patches(:,:,1);
            Yp = el_patches(:,:,2);
            Zp = el_patches(:,:,3);
            % Cp = patches(:,:,4); % apodization - not used
            
            % plot a patch: use depth as the color
            hp = patch(axs, 'XData', Xp, 'YData', Yp, 'ZData', Zp, 'CData', Zp, varargin{:});
            shading(axs, 'flat');
            
            % set default axis arguments
            xlabel(axs, 'x [m]');
            ylabel(axs, 'y [m]');
            zlabel(axs, 'z [m]');
            grid  (axs, 'on');
            zlim  (axs, [-1e-3, 1e-3] + [min(Zp(:)), max(Zp(:))])
            axis  (axs, 'equal')
            % axis  (axs, 'image');
            view  (axs, 3);
            
            % return
            if nargout >= 1, varargout{1} = hp; end
            if nargout >= 2, varargout{2} = axs; end
            
        end
        
        function varargout = plot(self, varargin)
            % parse inputs
            if nargin > 2 && isa(varargin{1}, 'matlab.graphics.axis.Axes')
                axs = varargin{1}; varargin(1) = [];
            else
                axs = gca;
            end

            % get positions
            p = self.getPositions();
            
            % plot
            hp = plot(axs, p(1,:), p(3,:), varargin{:});
            
            % return
            if nargout > 0, varargout{1} = hp; end
        end
        
        function varargout = surf(self, varargin)
            % see Transducer.patch
            varargout = cell([1, nargout]);
            [varargout{:}] = patch(self, varargin{:});
        end
    end
    
    % to be deprecated functions
    methods(Static, Access=private)
         function [y] = cgauspuls(t,fc,bw,bwr)
             %{
             % implement this, but with less checking and shaping
            [yi, yq] = gauspuls(varargin{:});
            y = complex(yi, yq);
             %}
             
             % defaults
             [bwrI, bwI, fcI] = deal(-6, 0.5, 1e3);

             % set inputs
             if nargin >= 4 && ~isempty(bwr), bwrI = bwr; end
             if nargin >= 3 && ~isempty(bw),  bwI = bw; end
             if nargin >= 2 && ~isempty(fc),  fcI = fc; end
             
             % Determine Gaussian mean and variance in the
             % frequency domain to match specifications:
             % r = 10.^(bwrI/20);             % Ref level (fraction of max peak)
             % fv = -bwI*bwI*fcI*fcI/(8*log(r)); % variance is fv, mean is fc
             
             % Determine corresponding time-domain parameters:
             % tv = 1./(4*pi*pi*fv);  % variance is tv, mean is 0

             % Compute time-domain pulse envelope, normalized by sqrt(2*pi*tv):
             % ye = exp(-t.*t./(2*tv));
             
             % Modulate envelope to form in-phase and quadrature components:
             % yc = ye .* cospi(2*fcI*t);  % In-phase
             % ys = ye .* sinpi(2*fcI*t);  % Quadrature

             % tv = ;  % variance is tv, mean is 0
             y = exp(-t.*t./(2*(1./(4*pi*pi*(-bwI*bwI*fcI*fcI/(8*log(10.^(bwrI/20)))))))) ...
                 .* (cospi(2*fcI*t) + 1i .* sinpi(2*fcI*t));
             
             % gauspuls converges to 0 at Inf
             y(isinf(t)) = 0;
             % idxInf = isinf(t);
             % yc(idxInf) = 0;
             % ys(idxInf) = 0;             
             
         end
         
         function f = cgauspulsfun(fc, bw, bwr)             
             isig = ((4*pi*pi*(-bw*bw*fc*fc / (8*log(10.^(bwr/20))))) ./ 2);
             f = @(t) exp(-t .* t .* isig) .* (cospi(2*fc*t) + 1i .* sinpi(2*fc*t));
         end         
         function f =  gauspulsfun(fc, bw, bwr)             
             isig = ((4*pi*pi*(-bw*bw*fc*fc / (8*log(10.^(bwr/20))))) ./ 2);
             f = @(t) exp(-t .* t .* isig) .* (cospi(2*fc*t)                      );
         end         
    end
end

% helper functions
function x = vec(x), x = x(:); end
