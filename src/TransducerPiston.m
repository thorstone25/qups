classdef TransducerPiston < Transducer

    properties
        radius               % the interelement distance (m)
    end
        
    % constructor and get/set methods
    methods(Access=public)
        % constructor: accept name/value pairs
        function self = TransducerPiston(varargin)
            
            % setup the transducer args
            if nargin == 1 && isa(varargin{1}, 'struct'), varargin = struct2nvpair(varargin{1}); end
            
            % initialize the (inherited) Transducer
            self@Transducer(varargin{:}) % assume we don't error on bad inputs
            
            % initialize the TransducerArray 
            for i = 1:2:numel(varargin)
                switch varargin{i}
                    case 'radius'
                        self.radius = varargin{i+1};                        
                end
            end            
        end
    end
    
    % define abstract methods
    methods
        function p = positions(self), p = self.offset; end
        function [theta, phi, normal, width, height] = orientations(self)
            [theta, phi, normal, width, height] = deal(0, 0, [0;0;1], [1;0;0], [0;1;0]);
        end
        function pb = bounds(self), pb = getBounds(self); end
        function pch = patches(self,sub_div), pch = getFieldIIPatches(self,sub_div); end
    end

    
    % define position methods
    methods    
        % get methods
        function pb = getBounds(self)
            % returns a 3 x 2 matrix of min / max values in x/y/z
            
            % transducer patches of {x,y,z,c} bound tuples
            pch = self.patches([1,1]);
            
            % get min/max bounds of the tx by iterating over each patch
            pb = [inf(3,1), -inf(3,1)];
            for i = 1:numel(pch)
                pchi = pch{i}(1:3);
                pb(:,1) = min(pb(:,1), cellfun(@(pch) min(pch, [], 'all'), pchi(:)));
                pb(:,2) = max(pb(:,2), cellfun(@(pch) max(pch, [], 'all'), pchi(:)));
            end
        end
    end
    
    % Field II conversion function
    methods(Access=public)
        function aperture = getFieldIIAperture(self, ~, el_sub_div)
            % creates a transmit / recieve aperture in Field II from the
            % transducer array
            
            % set defaults
            if nargin < 3
                % 1/10 lambda @ 1540
                element_size = self.radius * (1540 / self.fc / 10 / self.radius); 
            else
                % convert sub_divisions to element size
                element_size = self.radius / sqrt(prod(el_sub_div));
            end 
                        
            % Field II parameters
            xdc_piston_params = { ...
                self.radius, ...
                element_size, ...
                };
            
            % ensure double type
            xdc_piston_params = cellfun(@double, xdc_piston_params, 'UniformOutput', false);
            
            % Generate aperture for emission
            try evalc('field_info'); catch, field_init(-1); end
            aperture = xdc_piston(xdc_piston_params{:});
        end        
    end
    
    % USTB conversion function
    methods
        function probe = getUSTBProbe(self)
            % sub-aperture method
            %{
            pch = self.getFieldIIPatches();
            [p] = cellfun(@(pch) {cellfun(@(p) mean(p,'all'), pch(1:3))}, pch);
            [x, y, z] = cellfun(@(p) deal(p(1), p(2), p(3)), vec(p));
            [az, el, w, h] = deal(0, 0, self.radius, self.radius);
            geo = [x, y, z, repmat([az, el, w, h], [size(x, 1), 1])];
            %}

            % full aperture method
            geo = [0,0,0, 0,0, 2*self.radius,2*self.radius];
            probe = uff.probe(... 
                'geometry', geo,...
                'origin', uff.point('xyz', self.offset(:)')...
                );
        end
    end
            
    % dependents
    methods
        function set.radius(self, r)
            self.radius = r; 
            self.height = 2*r; 
            self.width  = 2*r; 
        end
    end
    
    methods(Static)
        function xdc = Panametrics()
            xdc = TransducerPiston(...
                'fc', mean([4e6 11e6]), ...
                'bandwidth', diff([4e6 11e6]), ...
                'numel', 1, ...
                'pitch', 0.2e-3, ...
                'focus', 20e-3 ...
                );
        end
    end
end
