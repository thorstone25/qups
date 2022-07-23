classdef Target < Medium & Scatterers
    
    methods
        function self = Target(varargin)
            
            % initialize Medium & Scatterers
            scatterer_args = {};
            medium_args = {};
            for i = 1:2:nargin
                switch varargin{i}
                    case {'pos','amp','c_scat', 'rho_scat', 'BoA_scat', 'alpha_scat', 'scat_mode','scat'}
                        scatterer_args{end+1} = varargin{i}; %#ok<AGROW>
                        scatterer_args{end+1} = varargin{i+1}; %#ok<AGROW>
            
                    case {'soundspeed','density','non-linearity','absorption','absorption-power','perturbation-regions', ...
                            'c0', 'rho0', 'boa0', 'alpha0', 'alphap0', 'pertreg'}
                        medium_args{end+1} = varargin{i}; %#ok<AGROW>
                        medium_args{end+1} = varargin{i+1}; %#ok<AGROW>
                end
            end
            
            self@Scatterers(varargin{:});
            self@Medium(varargin{:});
            
        end
        
        function target = copy(self)
            
            args = { ...
                'pos', self.pos, ...
                'amp', self.amp, ...
                'c_scat', self.c_scat, ...
                'rho_scat', self.rho_scat, ...
                'scat_mode', self.scat_mode, ...
                'soundspeed', self.c0, ...
                'density', self.rho0, ...
                'non-linearity', self.BoA0, ...
                'absorption', self.alpha0, ...
                'perturbation-regions', self.pertreg, ...
                };
            
            target = Target(args{:});
            
        end
        
        % get SIMUS medium params
        function p = getSIMUSParam(self)
            p = struct('c', self.c0); 
            if ~isnan(self.alpha0), p.attenuation = self.alpha0; end 
        end
        
        % get kWave compatible medium struct
        function kmedium = getKWaveMedium(self, kgrid, korigin)
            % kmedium = GETKWAVEMEDIUM(self, kgrid, korigin)
            %
            %   function to create a kWave compatible medium structure for
            %   a given Target.
            %
            %   Inputs:
            %     - kgrid:      a kWaveGrid object
            %
            %     - korigin:    a length kgrid.dim array of the origin of 
            %                   the kgrid in each axis e.g.  
            %                   kgrid.x_vec + korigin(1) == target_x_axis 
            %
            %   Outputs:
            %     - kmedium:    a kWave compatible medium object.
                        
            
            % get grid point sizing
            kgrid_size = max([kgrid.Nx, kgrid.Ny, kgrid.Nz], [1,1,1]);
            
            % allocate on gpu if possible
            try
                points = gpuArray.zeros([3, kgrid_size], 'single');
            catch ME
                disp(ME);
                points = zeros([3, kgrid_size], 'single');
            end
            
            % get points in UltrasoundSystem coordinates
            points(3,:) = korigin(1) + kgrid.x(:);
            points(1,:) = korigin(2) + kgrid.y(:);
            points(2,:) = korigin(3) + kgrid.z(:);
            
            % get property map from the Medium object
            [kmedium.sound_speed, kmedium.density, kmedium.BonA, kmedium.alpha_coeff, kmedium.alpha_power] = self.getPropertyMap(points);

            % create scatterers by perturbing the region at the nearest
            % pixel/voxel
            [ind_prop, c, rho, BonA, alpha_coeff] = self.modMap(points);
            
            % set the properties - there are much better ways of doing this
            kmedium.density(ind_prop)       = rho;
            kmedium.sound_speed(ind_prop)   = c;
            kmedium.BonA(ind_prop)          = BonA;
            kmedium.alpha_coeff(ind_prop)   = alpha_coeff;

            % remove higher order terms if the coefficients are all 0s
            if all(isnan(kmedium.alpha_coeff)), kmedium = rmfield(kmedium, ["alpha_coeff", "alpha_power"]); end
            if all(isnan(kmedium.BonA)), kmedium = rmfield(kmedium, "BonA"); end

            % set alpha power if alpha coefficient is set
            if isfield(kmedium, 'alpha_coeff'), kmedium.alpha_power = self.alphap0; end
        end
    
        % get Fullwave compatible map struct
        function maps = getFullwaveMap(self, scan, varargin)
            % GETFULLWAVEMAP - Get Fullwave compatible map structure
            %
            % maps = getFullwaveMap(self, grid) returns a map sampled on
            % the grid "grid". Some assembly required. Grid is a cell array
            % of vectors {x, y, z}.

            kwargs.method = "";
            for i = 1:2:numel(varargin), kwargs.(varargin{i}) = varargin{i+1}; end

            % create samples of the grid points
            pg = scan.getImagingGrid();
            pg = cellfun(@(x) {shiftdim(x, -1)}, pg); % 1 x X x Y x Z
            pg = cat(1, pg{:}); % points 3 x X x Y x Z

            % sample all maps on the grid points
            [c, rho, BoA, alpha] = getPropertyMap(self, pg); % 1 x X x Y x Z

            % create scatterers by perturbing the region at the nearest
            % pixel/voxel
            [ind_prop, cs, rhos, BonAs, alpha_coeffs] = self.modMap(pg);
            c(ind_prop)         = cs;
            rho(ind_prop)       = rhos;
            BoA(ind_prop)       = BonAs;
            alpha(ind_prop)     = alpha_coeffs;

            % set the map properties
            maps = struct('cmap', c, 'rmap', rho, 'amap', alpha, 'nmap', 1 + BoA./2);

            % Use 0 for invalid properties in fullwave(?)
            for f = string(fieldnames(maps))', maps.(f) = nan2zero(maps.(f)); end


        end
    
    end
end