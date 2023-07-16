% TRANSDUCER - Abstract Transducer class
%
% Superclass for a medical ultrasound transducer. Any TRANSDUCER has
% methods for the positions and orientations of it's elements as well as
% the characteristics of the elements (width/height, impulse response,
% central frequency, etc.). This class offers definitions for common
% transducers and conversion functions between real transducers (mainly
% from Verasonics) and simulation programs (k-Wave, Fullwave, FieldII).
%
% See also TRANSDUCERARRAY TRANSDUCERCONVEX TRANSDUCERMATRIX

classdef (Abstract) Transducer < matlab.mixin.Copyable
    properties
        fc (1,1) double = 5e6        % center frequency
        bw (1,2) double = [3.5e6 6.5e6] % bandwidth
        width (1,1) double = 1.5e-4 % width of an element
        height (1,1) double = 6e-3   % height of an element
        numel (1,1) double = 128     % number of elements
        offset (3,1) double = [0;0;0]% the offset from the origin
        impulse Waveform {mustBeScalarOrEmpty} = Waveform.empty() % the impulse response function of the element
    end

    properties
        el_focus (1,1) double = Inf % elevation focal depth
    end

    properties(GetAccess=public, SetAccess=protected, Dependent)
        area            % area of an element
    end

    properties(Dependent)
        bw_frac         % fractional bandwidth
    end
    properties(Dependent, Hidden)
        origin
    end

    % constructor
    methods
        function xdc = Transducer(kwargs)
            % TRANSDUCER - Transducer constructor
            %
            % xdc = TRANSDUCER(Name, Value, ...) constructs a
            % Transducer using name/value pair arguments.
            %
            % A Transducer is an abstract class and cannot be instantiated
            % directly.
            % 
            % See also TRANSDUCERARRAY TRANSDUCERCONVEX

            arguments
                kwargs.?Transducer
            end
            % name-value pair initialization
            % if width set but not height, default to 20x
            if isfield(kwargs,'width') && ~isfield(kwargs, 'height')
                kwargs.height = 20*kwargs.width;
            end

            % set all properties except bw_frac
            for f = setdiff(string(fieldnames(kwargs)), {'bw_frac'})', xdc.(f) = kwargs.(f); end

            % if frequency but not bandwidth or fractional bandwidth is
            % set, default to 60%
            if ~isfield(kwargs, 'bw') && ~isfield(kwargs, 'bw_frac')
                kwargs.bw_frac = 0.6;
            end

            % set fractional bandwidth last
            if isfield(kwargs, 'bw_frac'), xdc.bw_frac = kwargs.bw_frac; end


            % regardless of input, if impulse is empty, initialize it
            if isempty(xdc.impulse)
                xdc.impulse = xdc.ultrasoundTransducerImpulse();
            end
        end
    
        function s = obj2struct(xdc)
            % OBJ2STRUCT - Convert a QUPS object into a native MATLAB struct
            %
            % xdc = OBJ2STRUCT(xdc) converts the Transducer xdc and all of 
            % it's properties into native MATLAB structs.
            %
            % Example:
            %
            % % Create a Transducer
            % xdc = TransducerArray()
            %
            % % convert to a MATLAB struct
            % xdc = obj2struct(xdc)
            %
            arguments, xdc Transducer {mustBeScalarOrEmpty}; end
            
            W = warning('off', "MATLAB:structOnObject"); % squash warnings
            s = struct(xdc); % convert self
            if ~isempty(s), s.impulse = obj2struct(s.impulse); end % convert impulse
            s.class = class(xdc); % append class info
            warning(W); % restore warnings
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
            % % Create a Transducer
            % xdc = TransducerArray() % m, s, Hz
            %
            % % convert from meters to millimeters, hertz to megahertz
            % xdc = scale(xdc, 'dist', 1e3, 'time', 1e6) % mm, us, MHz
            %

            arguments
                xdc Transducer
                kwargs.dist (1,1) double
                kwargs.time (1,1) double
            end
            xdc = copy(xdc);
            if isfield(kwargs, 'dist')
                w = kwargs.dist;
                % scale distance (e.g. m -> mm)
                [xdc.width, xdc.height, xdc.offset, xdc.el_focus] = deal(w*xdc.width, w*xdc.height, w*xdc.offset, w*xdc.el_focus);
            end
            if isfield(kwargs, 'time')
                w = kwargs.time;
                % scale time (e.g. s -> us / Hz -> MHz)
                [xdc.fc, xdc.bw] = deal(xdc.fc/w, xdc.bw/w);
                xdc.impulse = scale(xdc.impulse, 'time', w);
            end
        end
    end

    % transducer specific methods
    methods (Abstract)
        % POSITIONS - Positions of the elements
        %
        % p = POSITIONS(xdc) returns a 3 x N array representing the 
        % positions of the center of each element.
        % 
        % See also ORIENTATIONS TRANSDUCER/PLOT
        p = positions(xdc); 

        % ORIENTATIONS - Orientation of the elements
        %
        % [theta, phi] = ORIENTATIONS(xdc) returns the azimuth and
        % elevation angles theta and phi of the elements of the transducer.
        % The angles are in degrees.
        %
        % [theta, phi, normal, width, height] = ORIENTATIONS(xdc) also
        % returns 3D vectors representing the normal vector, a vector with
        % the direction and magnitude of the width of the element, and a 
        % vector with the direction and mangitude of the height of the
        % element.
        %
        % See also POSITIONS TRANSDUCER/PLOT
        [theta, phi, normal, width, height] = orientations(xdc); % compute the orientations
    end        

    % define methods derived from position and orientation properties
    methods
        function pch = patches(xdc, sub_div)
            % PATCHES - Compute a matrix of sub-elements
            %
            % pch = PATCHES(xdc,sub_div) computes a cell matrix where each
            % element contains {X,Y,Z,C} tuples of 2x2 matrices specifying four
            % corners of each patch, where each patch is a subdivision of the
            % elements into sub-elements specified by the 1 x 2 array sub_div.
            % pch is a [Nel x Ndiv] matrix where Nel is the number of elements
            % and Ndiv is the number of sub-elements given by the
            % prod(sub_div).
            %
            % See also GETBARYCENTERS POSITIONS ORIENTATIONS
            arguments
                xdc Transducer
                sub_div (1,2) {mustBeInteger, mustBePositive} = [1,1]
            end
            % get the matrix of sub element patches
            % get unweighted difference in width and height from element center
            [dx, dy] = ndgrid(...
                ((0 : sub_div(1) - 1) / sub_div(1)) - 0.5, ...
                ((0 : sub_div(2) - 1) / sub_div(2)) - 0.5 ...
                ); % E1 x E2
            [dxs, dys] = deal(1 / sub_div(1), 1 / sub_div(2)); % sub difference

            E = prod(sub_div);

            % reshape to E x 4
            dx = dx(:); dy = dy(:);
            dx = [dx, dx + dxs, dx      , dx + dxs];
            dy = [dy, dy      , dy + dys, dy + dys];

            % move up to 1 x 1 x E x 4
            [dx, dy] = deal(shiftdim(dx, -2), shiftdim(dy, -2));

            % get the difference matrices (3 x N)
            [~, ~, ~, w, h] = xdc.orientations; % difference vectors

            % get the positions
            pn = xdc.positions;

            % get the element dimensions
            wel = xdc.width;
            hel = xdc.height;

            % get patch locations
            pc = pn + wel .* w .* dx + hel .* h .* dy; % 3 x N x E x 4
            pc(4,:) = 1; % x/y/z/c x N x E x 4 - set apodization to 1 
            % TODO: for soft baffle, this should be some distribution
            % instead
            
            % convert to patch structure {{2 x 2} x x/y/z/c} x 1 x N x E
            pc = permute(pc, [5,4,1,2,3]);  % 1 x 4 x x/y/z/c x N x E
            pc = reshape(pc, [2,2,size(pc,3:5)]); % 2 x 2 x x/y/z/c x N x E
            pch = cellfun(@(pc) shiftdim(num2cell(pc, [1:2]),1), num2cell(pc, 1:3), 'UniformOutput',false);
            pch = shiftdim(pch, 3);
                        
        end
    
        function pb = bounds(xdc)
            % BOUNDS - Compute the boundaries of the Transducer
            %
            % pb = BOUNDS(xdc) returns a 3 x 2 array of the minimum and maximum
            % cartesian coordinate in x/y/z for the Transducer xdc.
            %
            % See also POSITIONS ORIENTATIONS PATCHES
    
            arguments, xdc Transducer, end

            % transducer patches of {x,y,z,c} bound tuples
            pch = xdc.patches([1,1]);

            % get min/max bounds of the tx by iterating over each patch
            pb = [inf(3,1), -inf(3,1)];
            for i = 1:xdc.numel
                pchi = pch{i}(1:3);
                pb(:,1) = min(pb(:,1), cellfun(@(pch) min(pch, [], 'all'), pchi(:)));
                pb(:,2) = max(pb(:,2), cellfun(@(pch) max(pch, [], 'all'), pchi(:)));
            end
        end
    
        function p = getBaryCenters(xdc, el_sub_div)
            % GETBARYCENTERS - computes barycenters of the sub-elements
            % 
            % p = GETBARYCENTERS(xdc) computes the positions of the 
            % transducer sub-elements by computing the barycenters 
            % (geometric mean) of the corners of the sub-elements computed 
            % xdc.patches and returns them as a 3 x N matrix of positions
            % p.
            %
            % p = GETBARYCENTERS(xdc, sub_div) uses the 2-element
            % array sub_div to divide each element into 
            % sub_div(1) x sub_div(2) sub-elements and returns the 
            % positions as a 3 x N x E array of positions p, where 
            % E = prod(el_sub_div). The default is [1,1].
            %
            % See also PATCHES POSITIONS ORIENTATIONS

            arguments
                xdc Transducer
                el_sub_div (1,2) {mustBeInteger, mustBePositive} = [1,1]
            end

            el_patches = xdc.patches(el_sub_div);

            p = cellfun(@(patch) cellfun(@(comp) mean(comp(:)),...
                patch(1:3)'), el_patches, 'UniformOutput', false);
            p = reshape(cat(1, p{:}), [3, size(p)]);
        end
    end

    % toolbox conversion functions
    methods (Abstract)
        % GETFIELDIIAPERTURE - Create a FieldII aperture object
        %
        % ap = GETFIELDIIAPERTURE(xdc) creates a FieldII aperture ap
        % from the Transducer xdc. If FieldII is not initialized before
        % calling this function, it will be initiated.
        %
        % ap = GETFIELDIIAPERTURE(xdc, sub_div) divides the width and 
        % height of each element by into a sub_div(1) x sub_div(2)
        % sub-elements. The default is [1 1].
        %
        % ap = GETFIELDIIAPERTURE(xdc, sub_div, focus) creates an aperture
        % with an elevational focus at the point focus. A focal depth of
        % Inf represents no focus. The default is [0 0 Inf].
        % 
        % When creating the aperture, an infinite focus in instead set to
        % realmax('single').
        %
        % See also ULTRASOUNDSYSTEM.CALC_SCAT_ALL, FIELD_INIT
        aperture = getFieldIIAperture(xdc, sub_div, focus);

        % QUPS2USTB - Create a USTB compatible uff.probe object
        %
        % probe = QUPS2USTB(xdc) creates a uff.probe probe from the
        % Transducer xdc. USTB must be on the path.
        %
        % Example:
        %
        % probe = QUPS2USTB(TransducerArray.L12_3v());
        % 
        % See also UFF.PROBE
        probe = QUPS2USTB(xdc); % get the USTB probe object
    end

    % Verasonics conversion functions
    methods (Static)
        function xdc = Verasonics(Trans, c0)
            % VERASONICS - Construct a Transducer from a Verasonics struct
            %
            % xdc = Transducer.VERASONICS(Trans) constructs a Transducer 
            % from the properties defined in the Verasonics 'Trans' struct.
            %
            % xdc = Transducer.VERASONICS(Trans, c0) uses c0 in m/s as the 
            % sound speed when converting from wavelengths to meters. This
            % is typicaly set by the Verasonics property
            % 'Resource.Parameters.speedOfSound'. Be sure to explicitly set
            % this if other than 1540. 
            %
            % Example:
            %
            % Trans = struct('name', 'L12-3v', 'units', 'mm');
            % Trans = computeTrans(Trans);
            % xdc = Transducer.Verasonics(Trans);
            %
            % See also SCAN.VERASONICS TRANSDUCER.UFF TRANSDUCER.QUPS2USTB
            arguments
                Trans (1,1) struct
                c0 (1,1) double = 1540
            end

            % get distance property scaling
            switch Trans.units
                case 'wavelengths', scale = c0 / Trans.frequency * 1e-6; % lambda -> m
                case 'mm', scale = 1e-3; % mm -> m
            end

            % set the element length: search for an indication of the element height
            if isfield(Trans, 'elementLength') % best - leave as is
            elseif isfield(Trans, 'elevationApertureMm') % backup
                Trans.elementLength = Trans.elevationApertureMm * (1e-3 / scale); 
            else % fallback - make it square
                Trans.elementLength = Trans.elementWidth; 
            end

            % set required arguments
            switch Trans.type % 0=Lin(y=z=0),1=CrvdLn(y=0),2=2D,3=ann,4=R/C
                case 0,     xdc = TransducerArray.Verasonics(Trans  , c0);
                case 1,     xdc = TransducerConvex.Verasonics(Trans , c0);
                case {2,4}, xdc = TransducerMatrix.Verasonics(Trans , c0);
                otherwise,  xdc = TransducerGeneric.Verasonics(Trans, c0);
            end
            
            % set optional arguments
            % parse the impulse response
            if isfield(Trans, 'IR1wy')
                h = Trans.IR1wy; % impulse response
                t0 = - (argmax(hilbert(h))-1) / 250e6; % offset to peak time
                xdc.impulse = Waveform('t', t0 + (0:numel(h)-1) / 250e6, 'samples',h); % impulse response
            else
                xdc.impulse = Waveform.Delta();
            end

            % set the elevation focus
            if isfield(Trans, 'elevationFocusMm'), xdc.el_focus = 1e-3 * Trans.elevationFocusMm; end

            % set the bandwidth
            if isfield(Trans, 'Bandwidth'), xdc.bw = 1e6*Trans.Bandwidth([1 end]); end
        end
    end

    % SIMUS functions
    methods (Abstract)
        % GETSIMUSPARAM - Create a MUST compatible parameter struct
        %
        % p = GETSIMUSPARAM(xdc) returns a structure with properties
        % for a call to simus().
        %
        % See also ULTRASOUNDSYSTEM/SIMUS
        p = getSIMUSParam(xdc)
    end

    % UFF constructor
    methods(Static)
        function xdc = UFF(probe)
            % UFF - Construct a Transducer from a UFF probe
            %
            % xdc = TRANSDUCER.UFF(probe) converts the uff.probe probe to a
            % Transducer xdc.
            %
            % See also TRANSDUCER.QUPS2USTB
            arguments, probe uff.probe; end
            switch class(probe)
                case 'uff.linear_array',     xdc = TransducerArray.UFF(probe);
                case 'uff.curvilinear_array',xdc = TransducerConvex.UFF(probe);
                case 'uff.matrix_array',     xdc = TransducerMatrix.UFF(probe);
                case 'uff.probe',            xdc = TransducerGeneric.UFF(probe);
                otherwise,                   xdc = TransducerGeneric.UFF(probe);
            end
        end
    end

    % kWave functions (old)
    methods(Hidden)
        function [ksensor, ksensor_ind, sens_map] = getKWaveSensor(xdc, kgrid, kgrid_origin, el_sub_div)
            % GETKWAVESENSOR - Get a kWave compatible sensor struct
            % 
            % ksensor = GETKWAVESENSOR(xdc, kgrid, kgrid_origin)
            % 
            % TODO: doc this function: it gives you a k-Wave ready sensor
            % structure
            %             arguments
            %                 xdc (1,1)
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
            [width_, height_] = deal(xdc.width, xdc.height);

            % get regions for each receive transducer element
            el_cen = xdc.positions(); % center of each element
            [theta, phi, el_dir, el_wid, el_ht] = xdc.orientations(); % orientations vectors for each element

            % reorder to map between kwave and conventional ultrasound coordinates
            [el_dir, el_cen, el_wid, el_ht] = dealfun(@(v)v([3 1 2],:), el_dir, el_cen, el_wid, el_ht);

            % get edges in x/y/z (i.e. patches) for each sub element
            patches = xdc.patches(el_sub_div)'; % (E x N)
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
            parfor el = 1:xdc.numel
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
            parfor (el = 1:xdc.numel, 0)
                [~, ksensor_ind{el}] = ismember(find(ksensor{el}.mask), ind_full_sensor);
            end

            % ensure results on the host
            ksensor_ind  = cellfun(@gather, ksensor_ind, 'UniformOutput', false);
            ksensor = cellfun(@(k)struct('mask',gather(k.mask)), ksensor, 'UniformOutput', false);
        end
    end

    % FDTD functions
    methods
        function [mask, el_weight, el_dist, el_ind] = elem2grid(xdc, scan, el_sub_div)
            % ELEM2GRID - Transducer element to grid mapping
            %
            % [mask, el_weight, el_dist, el_ind] = ELEM2GRID(xdc, scan)
            % returns a binary mask and a set of indices, weights, and
            % distances defined for each element give a Transducer self and
            % a ScanCartesian scan. The element indicies are defined on the
            % vectorized non-zero indices of the mask i.e. on the indices
            % of `find(mask)`.
            %
            % [...] = ELEM2GRID(xdc, scan, el_sub_div) subdivides the
            % elements in width and height by el_sub_div. The default is
            % [1,1].
            % 
            % See also ULTRASOUNDSYSTEM/KSPACEFIRSTORDER

            arguments
                xdc (1,1) Transducer
                scan (1,1) ScanCartesian
                el_sub_div (1,2) {mustBeInteger, mustBePositive} = [1,1]
            end

            % get the sensor and source masks. This is the hard part: how 
            % do I do this for a convex probe on a grid surface?

            vec = @(x)x(:); % define locally for compatibility

            % cast grid points to single type for efficiency and base grid on the origin
            [gxv, gyv, gzv] = deal(single(scan.x), single(scan.y), single(scan.z)); % grid {dim} vector
            [gx0, gy0, gz0] = deal(scan.x(1), scan.y(1), scan.z(1)); % grid {dim} first point
            [gxd, gyd, gzd] = deal(single(scan.dx), single(scan.dy), single(scan.dz)); % grid {dim} delta
            pg = scan.positions(); % 3 x Z x X x Y
            pdims = arrayfun(@(c){find(c == scan.order)}, 'XYZ'); % dimensions mapping
            [xdim, ydim, zdim] = deal(pdims{:});
            iord = arrayfun(@(d)find([pdims{:}] == d), [1,2,3]); % inverse mapping

            % get array sizing - size '1' and step 'inf' for sliced dimensions
            [Nx, Ny, Nz] = deal(scan.nx, scan.ny, scan.nz);

            % get local element size reference
            [width_, height_, sz_] = deal(xdc.width, xdc.height, scan.size);

            % get regions for each receive transducer element
            el_cen = xdc.positions(); % center of each element
            [theta, phi, el_dir, el_wid, el_ht] = xdc.orientations(); % orientations vectors for each element

            % reorder to map between kwave and conventional ultrasound coordinates
            % [el_dir, el_cen, el_wid, el_ht] = dealfun(@(v)v([3 1 2],:), el_dir, el_cen, el_wid, el_ht);

            % get edges in x/y/z (i.e. patches) for each sub element
            patches = xdc.patches(el_sub_div)'; % (E x N)
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
            parfor el = 1:xdc.numel
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
            find_mask = find(mask);
            parfor m = 1:numel(mask_ind) %#ok<CPROPLC> 
                ind_el{m} = gather(argn(2, @ismember, mask_ind{m}, find_mask));
            end
            ind_el = reshape(ind_el, size(mask_ind));
            
            % make all into arrays, defined by mask indices
            parfor el = 1:xdc.numel
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

    % kWaveArray constructor
    methods
        function karray = kWaveArray(xdc, dim, og, karray_args)
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
            arguments
                xdc Transducer
                dim {mustBeInteger, mustBeInRange(dim, 2,3)}
                og (:,1) {mustBeNumeric}
                karray_args.Axisymmetric (1,1) logical
                karray_args.UpsamplingRate (1,1) double
                karray_args.BLITolerance (1,1) double {mustBeInRange(karray_args.BLITolerance, 0, 1)}
                karray_args.BLIType (1,1) string {mustBeMember(karray_args.BLIType, ["sinc", "exact"])}
            end

            % TODO: set axisymmetric if it makes sense
            % TODO: allow for elevation

            % initialize
            if isfield(karray_args,'BLIType')
                karray_args.BLIType = char(karray_args.BLIType);
            end
            karray_args = struct2nvpair(karray_args);
            karray = kWaveArray(karray_args{:});

            % get positions and orientations
            p = xdc.positions; % TODO: translate the entire array?
            [th, phi] = xdc.orientations;
            n = xdc.numel;
            if any(phi), warning('kWaveArray initialization with elevation angle is untested.'); end

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
            end

            % add each element to the array
            for i = 1:n
                karray.addRectElement(p(:,i), xdc.width, xdc.height, rot(:,i));
            end

            karray.setArrayPosition(arr_off, arr_rot);
        end
    end

    % fullwave functions
    methods (Abstract)
        % GETFULLWAVETRANSDUCER - define a Fullwave transducer structure
        %
        % xdc_fw = GETFULLWAVETRANSDUCER(xdc, scan) creates a fullwave
        % compatible structure xdc_fw defined on the ScanCartesian scan.
        %
        % xdc_fw is a Fullwave compatible struct with (at least) the
        % following fields:
        %       - npx           number of elements
        %       - inmap         mask of the input pixels
        %       - nInPx         number of input pixels
        %       - nOutPx        number of output pixels
        %       - incoords      (x,y,1,el) coordinate pairs of the input pixels
        %       - outcoords     (x,y,1,el) coordinate pairs of the output pixels
        %
        % See also ULTRASOUNDSYSTEM/FULLWAVESIM
        xdc_fw = getFullwaveTransducer(xdc, scan)
    end

    % Field II calls - rely on being able to get a FieldII aperture
    methods
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
        end

        function el_patches = getFieldIIPatches(xdc, el_sub_div)
            % GETFIELDIIPATCHES - get cell array of patches
            %
            % el_patches = GETFIELDIIPATCHES(xdc) returns the
            % corners of the mathematical elements represented in FieldII.
            % Each element is represented by a {x,y,z,c} tuples (a.k.a. a
            % "patch") which holds [x1, x2; x3 x4], which are the corners
            % of the elements.
            %
            % el_patches = GETFIELDIIPATCHES(xdc, el_sub_div) uses the 
            % 2-element array el_sub_div to determine the number of 
            % mathematical elements. The result is a (N x E) array where N
            % is the number elements and E is the number of subelements.
            %
            % See also GETFIELDIIBARYCENTERS
            arguments
                xdc Transducer
                el_sub_div (1,2) double {mustBeInteger, mustBePositive} = [1,1]
            end

            % get the aperture
            ap = xdc.getFieldIIAperture([0 0 0], el_sub_div);

            % extract corners
            data = xdc_get(ap, 'rect');
            ind = 1:size(data,2);

            % release the aperture
            xdc_free(ap);

            % extract the corners into {x,y,z,c} tuples (a.k.a. a patch)
            % where x/y/z/c are each 2 x 2 arrays of points
            corner_map = @(i){...
                xdc.offset(1) + ...
                [   data(11,i), data(20,i); ...
                data(14,i), data(17,i)], ...
                xdc.offset(2) + ...
                [   data(12,i), data(21,i); ...
                data(15,i), data(18,i)], ...
                xdc.offset(3) + ...
                [   data(13,i), data(22,i);...
                data(16,i), data(19,i)], ...
                data(5, i) * ones(2,2),...
                };
            el_patches = arrayfun(corner_map, ind, 'UniformOutput', false); % array of {x,y,z,c} tuples

            % reshape by element (N x E), E is subelements
            el_patches = reshape(el_patches, prod(el_sub_div), [])';
        end
    end

    % internal subroutines
    methods
        function impulse = ultrasoundTransducerImpulse(xdc)
            % ULTRASOUNDTRANSDUCERIMPULSE - create an impulse response Waveform
            % 
            % impulse = ULTRASOUNDTRANSDUCERIMPULSE(xdc) creates a gaussian
            % pulse Waveform with the bandwidth and fractional bandwidth of
            % the Transducer xdc.
            % 
            % See also WAVEFORM.DELTA()

            % defaults
            arguments
                xdc Transducer
            end

            % get energy strength cutoff time
            bwr = -6;
            tpr = -80;

            tc = gauspuls('cutoff',xdc.fc, xdc.bw_frac, bwr, tpr);
            xfc = xdc.fc;
            bwf = xdc.bw_frac;

            % get the function
            impulse_fun = Transducer.cgauspulsfun(xfc, bwf, bwr);

            % make a Waveform object
            impulse = Waveform('fun', impulse_fun, 't0', -tc, 'tend', tc);
        end
    end

    % dependent methods
    methods
        function set.bw(self, b)
            if(length(b) ~= 2)
                error("bw must be a length 2 vector of the passband cutoff frequencies.");
            end
            self.bw = b;
        end

        function a = get.area(self), a = self.width * self.height; end

        function b = get.bw_frac(self), b = range(self.bw) ./ self.fc; end

        function set.bw_frac(self, bf), self.bw = mean(self.fc) + mean(self.fc) * bf * ([-1 1] ./ 2); end

        function o = get.origin(self), o = - self.offset; end

        function set.origin(self, o), self.offset = -o; end
    end

    % plot functions
    methods
        function hp = patch(self, varargin, patch_args, xdc_args)
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
            % Inputs
            %   - el_sub_div:    2 x 1 vector of element subdivisions
            %
            % See also PATCH

            arguments
                self (1,1) Transducer
            end
            arguments(Repeating)
                varargin
            end
            arguments
                patch_args.?matlab.graphics.primitive.Patch
                patch_args.DisplayName = 'Elements'
                xdc_args.el_sub_div (1,2) double {mustBeInteger, mustBePositive} = [1 1];
            end

            % extract axes
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes') %#ok<CPROPLC> use built-in numel
                axs = varargin{1}; varargin(1) = [];
            else, axs = gca;
            end

            % get the aperture as a cell array of {x,y,z,c} pairs
            el_patches = self.patches(xdc_args.el_sub_div);

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
            % Cp = el_patches(:,:,4); % apodization - not used
            % TODO: add option for showing apodization weighting

            % plot a patch: use depth as the color
            args = struct2nvpair(patch_args);
            hp = patch(axs, 'XData', Xp, 'YData', Yp, 'ZData', Zp, 'CData', Zp, varargin{:}, args{:});

            % set default axis arguments
            xlabel(axs, 'x');
            ylabel(axs, 'y');
            zlabel(axs, 'z');
            grid  (axs, 'on');
            % zlim  (axs, [-1e-3, 1e-3] + [min(Zp(:)), max(Zp(:))])
            axis  (axs, 'equal')
            view  (axs, 3);
            
        end

        function varargout = plot(self, varargin, plot_args)
            % PLOT - overload the plot function
            %
            % PLOT(xdc) plots the locations of the Transducer xdc.
            %
            % PLOT(xdc, ax) uses the axes ax instead of the current axes.
            %
            % PLOT(..., Name, Value, ...) passes name-value pair arguments
            % to the built-in plot function so that name value pairs that
            % are valid for plot are valid here.
            %
            % h = PLOT(...) returns the handle to the plot.
            %
            % Plots only the x-z slice.
            %
            % See also TRANSDUCER/PATCH 

            arguments
                self (1,1) Transducer
            end
            arguments(Repeating)
                varargin
            end
            arguments
                plot_args.?matlab.graphics.chart.primitive.Line
                plot_args.DisplayName = 'Elements'
                plot_args.LineStyle = 'none'
                plot_args.Marker = 'square'
            end
            
            % extract axes
            if numel(varargin) >= 1 && isa(varargin{1},'matlab.graphics.axis.Axes') %#ok<CPROPLC> use built-in numel 
                axs = varargin{1}; varargin(1) = []; 
            else, axs = gca;
            end

            % get positions
            p = self.positions();

            % plot
            plot_args = struct2nvpair(plot_args);
            hp = plot(axs, p(1,:), p(3,:), varargin{:}, plot_args{:});

            % return
            if nargout > 0, varargout{1} = hp; end
        end

        function varargout = surf(self, varargin)
            % SURF - Alias for patch
            % 
            % See also TRANSDUCER/PATCH
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
