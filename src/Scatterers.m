classdef Scatterers < handle
   
    properties
        pos = zeros([3, 0])     % positions for scatterers (m) as a 3 x S vector
        amp = ones([1, 0])      % amplitude for the scatterers (linear) as a 1 x S vector
        c_scat = ones([1, 0])   % scatterring speed of sound as a 1 x S vector
        rho_scat = ones([1, 0]) % scatterring density as a 1 x S vector
        BoA_scat = ones([1, 0]) % scatterring non-linearity as a 1 x S vector
        alpha_scat = ones([1, 0]) % scatterring attenuating factor power as a 1 x S vector
        alphap_scat = ones([1, 0])  % scatterring attenuating power as a 1 x S vector
        % mode of scattering, either absolute, or
        % proportional to medium reference properties
        % {'ratio*' | 'abs'}
        scat_mode = 'ratio'
    end

    properties(Dependent, Hidden)
        xb  % x bounds as a 1 x 2 vector
        yb  % y bounds as a 1 x 2 vector
        zb  % z bounds as a 1 x 2 vector
    end
    
    % constructor and get/set methods
    methods
        function self = Scatterers(varargin)
            % name value pairs
            for i = 1:2:nargin
                switch varargin{i}
                    case {'pos', 'amp', 'c_scat', 'rho_scat', ...
                            'BoA_scat', 'alpha_scat', 'alphap_scat', 'scat_mode'}
                        self.(varargin{i}) = varargin{i+1};
                end
            end

            % infer data sizing
            Ss = ([numel(self.amp), cellfun(@(x) size(x,2), {self.pos, self.c_scat, self.rho_scat, self.BoA_scat, self.alpha_scat, self.alphap_scat})]);
            if any(Ss) % we have some input - infer sizing and reinitialize
                S = unique(Ss(Ss ~= 0));
                assert(isscalar(S), 'Cannot infer number of scatterers - check input sizes!');

                % initialize non-initialized inputs
                switch self.scat_mode % choose identity value for the operation
                    case 'ratio', id = @ones;
                    case 'abs',   id = @zeros; % this doesn't always make sense, but allows init to go through
                end

                if isempty(self.amp),           self.amp         = id([1,S]); end
                if isempty(self.c_scat),        self.c_scat      = id([1,S]); end
                if isempty(self.rho_scat),      self.rho_scat    = id([1,S]); end
                if isempty(self.BoA_scat),      self.BoA_scat    = id([1,S]); end
                if isempty(self.alpha_scat),    self.alpha_scat  = id([1,S]); end
                if isempty(self.alphap_scat),   self.alphap_scat = id([1,S]); end
            end
            

        end
        function n = numScat(self), n = numel(self.amp); end
        function bounds = getBounds(self), bounds = cat(1, self.xb, self.yb, self.zb); end
    end

    %
    methods
        function [ind_prop, c, rho, BonA, alpha_coeff] = modMap(self, points)
            % MODMAP - Modification map 
            % 
            % [ind_prop, c, rho, BonA, alpha_coeff] = MODMAP(self, points)
            % returns the indices and properties to modify a Medium at the 
            % given points to include the scatterers and their scatterering
            % properties.
            %
            % See also MEDIUM/GETPROPERTYMAP

            % create scatterers by perturbing the region at the nearest
            % pixel/voxel
            [pos] = deal(self.pos); %#ok<PROPLC> % splice
            ind_prop = nan([1, self.numScat()], 'double');
            vec = @(x) x(:);
            parfor (j = 1:self.numScat(), 0)
                ind_prop(j) = gather(argmin(vec(sum((points - pos(:,j)) .^ 2, 1)))); %#ok<PROPLC> 
            end
            
            % get the density/sound speed for the index
            switch self.scat_mode
                case 'ratio'
                    new_prop(4,:) = self.alpha_scat .* self.alpha0;
                    new_prop(3,:) = self.BoA_scat .* self.BoA0;
                    new_prop(2,:) = self.c_scat .* self.c0;
                    new_prop(1,:) = self.rho_scat .* self.rho0;
                case 'abs'
                    new_prop(4,:) = self.alpha_scat;
                    new_prop(3,:) = self.BoA_scat;
                    new_prop(2,:) = self.c_scat;
                    new_prop(1,:) = self.rho_scat;
                otherwise
                    error('Unrecognized scattering mode.');
            end
            
            % set output properties
            [c, rho, BonA, alpha_coeff] = deal(...
                new_prop(2,:), ...
                new_prop(1,:), ...
                new_prop(3,:), ...
                new_prop(4,:) ...
                );
        end
    end
    
    % plot methods
    methods 
        function h = plot(self, axs, plot_args)
            % PLOT - overload the plot function
            % 
            % PLOT(self) plots the locations of the Scatterers self.
            %
            % PLOT(self, ax) uses the axes ax instead of the current axes.
            % 
            % PLOT(..., Name, Value, ...) passes name-value pair arguments
            % to the built-in plot function so that name value pairs that 
            % are valid for plot are valid here. 
            % 
            % h = PLOT(...) returns the handle to the plot.
            % 
            % Plots only the x-z slice.
            % 
            % See also MEDIUM/IMAGESC
            arguments
                self (1,1) Scatterers
                axs (1,1) matlab.graphics.axis.Axes = gca
                plot_args.?matlab.graphics.chart.primitive.Line
            end

            % plot
            plot_args = struct2nvpair(plot_args);
            h = plot(axs, self.pos(1,:), self.pos(3,:), plot_args{:});
        end        
    end
    
    % dependent variables
    methods
        function b = get.xb(self), b = [min(self.pos(1,:)), max(self.pos(1,:))]; end
        function b = get.yb(self), b = [min(self.pos(2,:)), max(self.pos(2,:))]; end
        function b = get.zb(self), b = [min(self.pos(3,:)), max(self.pos(3,:))]; end
    end    
end