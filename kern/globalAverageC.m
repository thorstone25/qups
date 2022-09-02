function c_avg = globalAverageC(med, scan, pa, pb)
% GLOBALAVERAGEC - Global average sound speed between a set of points.
% 
% c_avg = GLOBALAVERAGEC(med, scan, pa, pb) returns the global average
% sound speed from the point pa to all points pb in the Medium med defined 
% on the ScanCartesian scan. The points pb and pa must all lie on the 2D
% X-Z plane. 
%
% Points outside of the scan are NaN.
%
% Example:
% 
% % Create a scan
% scan = ScanCartesian('x', linspace(-5, 5, 101), 'z', linspace(0, 40, 401));
% 
% % Create a simple layered sound speed distribution
% C = 1.5*ones(scan.size); % first layer
% C(scan.z > 20, :) = 1.6; % second layer
% C(scan.z > 30, :) = 1.4; % third layer
% 
% % Construct a Medium
% med = Medium.Sampled(scan, C); % create a medium
% med.c0 = 3.0; % ambient sound speed
% 
% % Compute the averege sound speed
% p_foc = [0;1] .* scan.z(:)'; % focal points
% c_avg = globalAverageC(med, scan, [0;0], p_foc);
%
% % Plot the average sound speed as a function of depth
% figure;
% plot(p_foc(2,:), c_avg, '.-');
% grid on;
% grid minor;
% xlabel('Depth (mm)')
% ylabel('Sound Speed (mm/us)');
% title('Layered Model Average Sound Speed')
% 
% See also RAYPATHS

arguments
    med (1,1) Medium
    scan (1,1) ScanCartesian 
    pa (2,1) {mustBeNumeric} % aperture positions in 2D
    pb (2,:) {mustBeNumeric} % endpoint positions in 2D
end

% invalidate OOB points
fval = [scan.xb(1); scan.zb(1)] <= pb & pb <= [scan.xb(2); scan.zb(2)]; % OOB?
pb(:, ~all(fval, 1)) = nan; % invalidate if either dim is OOB

% sample the entire sound speed map on the grid points
cg = props(med, scan); % size(cg) <-> scan.size

% compute the raypaths from the point pa to the points pb on
% the grid points given by grid
[~, xzind] = ismember('ZX', scan.order); % get the correct order for Z and X 
[w, l] = rayPaths(scan.x, scan.z, pa, pb, 'ord', scan.order(xzind)); % [X x Z] x J

% compute the average sound speed from A to B using the
% integral of the inverse sound speed times each grid cell
% length. This can be subject to numerical precision issues
c_avg = 1 ./ (((1 ./ cg(:))' * w) ./ l); % 1 x J integral of slowness dl
c_avg(isinf(c_avg)) = nan; % disallow infinite speeds



