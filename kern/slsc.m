function z = slsc(x, dim, L, method, kdim)
% SLSC - Short-lag Spatial Coherence 
%
% z = SLSC(x) computes the short-lag spatial coherence (SLSC) across the
% data ND-array x.
%
% z = SLSC(x, dim) operates in dimension dim. The default is the last 
% non-singleton dimension.
% 
% z = SLSC(x, dim, L) uses up to L lags. If L is an array of indices, the
% indices matching L are used. Select a single lag by making an array
% including 0 e.g. L = [0 5]. The default is 
% max(1, floor(1/4 * size(x,dim)));
% 
% z = SLSC(x, dim, L, method) uses the specified method. Must be one
% of {"ensemble" | "average"*}.
% 
% z = SLSC(x, dim, L, method, kdim) interprets kdim as the time kernel
% dimension. When dimension kdim is non-singular, samples in time are
% treated as part of the same correlation.
% 
% About: 
% The short-lag spatial coherence measures how similar signals are across a
% short distance. In medical ultrasound, speckle will exhibit low coherence
% over a short region whereas tissue walls may exhibit high spatial
% coherence.
% 
% References:
% D. Hyun, A. C. Crowley and J. J. Dahl, 
% "Efficient Strategies for Estimating the Spatial Coherence of Backscatter,"
% in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
% vol. 64, no. 3, pp. 500-513, March 2017.  
% doi: <a href="matlab:web('doi.org/10.1109/TUFFC.2016.2634004')">10.1109/TUFFC.2016.2634004</a>
%
% M. A. Lediju, G. E. Trahey, B. C. Byram and J. J. Dahl, 
% "Short-lag spatial coherence of backscattered echoes: imaging characteristics," 
% in IEEE Transactions on Ultrasonics, Ferroelectrics, and Frequency Control, 
% vol. 58, no. 7, pp. 1377-1388, July 2011.  
% doi: <a href="matlab:web('doi.org/10.1109/TUFFC.2011.1957')">10.1109/TUFFC.2011.1957</a>
%
% Example 1:
% % This example requires kWave
% if ~exist('kWaveGrid', 'class')
%   warning('kWave must be on the path to run this example.');
% end
%
% % Setup a system
% sscan = ScanCartesian(...
%   'x', 1e-3*linspace(-20, 20, 1+40*2^3), ...
%   'z', 1e-3*linspace(-02, 38, 1+40*2^3) ...
%   );
% xdc = TransducerArray();
% seq = SequenceRadial('angles', -21 : 1.5 : 21, 'ranges', 1); % plane-waves
% us = UltrasoundSystem('scan', sscan, 'xdc', xdc, 'sequence', seq);
% 
% % Create a Medium to simulate
% [c, rho] = deal(1500*ones(sscan.size), 1000*ones(sscan.size));
% rho = rho + 15 * (rand(sscan.size)-0.5); % 1% variation in the density
% [Xg, ~, Zg] = sscan.getImagingGrid();
% for zj = 5:5:30
%   rho(Xg == 0 & Zg == 1e-3*zj) = 1000*2; % double the density at these depths
% end
% med = Medium.Sampled(sscan, c, rho, 'c0', 1500);
% 
% % Simulate the ChannelData
% chd = kspaceFirstOrder(us, med, sscan);
% 
% % Pre-process the data
% chd = hilbert(chd);
%
% % generate an image, preserving the receive dimension
% us.scan.dz = us.lambda / 4; % axial resolution
% us.scan.dx = us.lambda / 4; % lateral resolution
% b = DAS(us, chd, 'keep_rx', true); % preserve the receive aperture dimension
% mdim = ndims(b); % the last dimension is the receive aperture dimension
% 
% % compute and show the SLSC image
% z = slsc(b, mdim);
% figure;
% imagesc(us.scan, real(z));
% colorbar;
%
% % Example 2:
% % image with a small time kernel
% K = 3; % one-sided time kernel size
% t0 = chd.t0; % original start time
% dt = 1/chd.fs; % time shift
% for k = -K:K
%     chd.t0 = t0 + k*dt; % shift time axis
%     bk{k+K+1} = DAS(us, chd, 'keep_rx', true); % preserve the receive aperture dimension
% end
% chd.t0 = t0; % restore time axis
% mdim = ndims(bk{1}); % receive aperture dimension
% kdim = ndims(bk{1}) + 1; % find a free dimension
% bk = cat(kdim, bk{:});
%
% % compute SLSC with a time kernel
% L = floor(us.rx.numel/4); % maximum lag
% zk = slsc(bk, mdim, L, "average", kdim);
%
% % show the SLSC image
% figure;
% imagesc(us.scan, real(zk));
% colorbar;
% 
% See also DMAS COHFAC 

% defaults
arguments
    x
    dim (1,1) {mustBePositive, mustBeInteger} = prod(find(size(x,1:ndims(x)) ~= 1, 1, 'last')); % lag dimension
    L (1,:) {mustBeNonnegative, mustBeInteger} = max(1, floor(size(x,dim)/4)); % lags or maximum lag
    method (1,1) string {mustBeMember(method, ["average", "ensemble"])} = "average"; % average or ensemble estimator
    kdim (1,1) {mustBePositive, mustBeInteger} = ndims(x)+1; % time samples dimension
end

% default parpool: this is data heavy, so only use a threadpool
% don't use any pool for gpu data
if isa(gcp('nocreate'), 'parallel.threadPool') ...
        && ~isa(x, 'gpuArray') && ~isa(x, 'tall') 
    clu = gcp(); else, clu = 0; 
end

% get weighting filter across receiver pairs
A = gather(size(x,dim)); % the full size of the aperture
[M, N] = ndgrid(1:A, 1:A);
H = abs(M - N); % lags for each receiver/cross-receiver pair (M x M')
if isscalar(L), lags = 1:L; else, lags = L; end % chosen lags for adding (i.e. the short lags)
S = ismember(H,lags); % selection mask for the lags (M x M') - true if (rx, xrx) is a short lag
L = numel(lags); % number of (active) lags - this normalizes the sum of the average estimates

% choose average or ensemble
switch method
    case "average"
        % normalize magnitude per time sample (averaging) with 0/0 -> 0
        x = nan2zero(x ./ vecnorm(x,2,kdim));

        % get weighting / filter across receiver pairs
        W = S ./ (A - H) / 2 / L; % weights per pair (debiased, pos & neg averaged, multi-correlation-averaged)

        % place weights as cross receiver over cells, receiver in dim
        W = num2cell(shiftdim(W, -(dim-2)), dim); % ... x (M') x {M}

        % compute product and sum per cross-receiver
        z = 0;
        parfor (i = 1:A, clu)
            xc = sub(x,i,dim); % cross-receiver
            z = z + sum(W{i} .* conj(xc) .* x, [dim, kdim], 'omitnan');
        end

    case "ensemble"
        % scale the problem to maintain numerical accuracy
        % we do this as precisely as possible by shifting by a power of 2.
        x = x .* 2.^nextpow2(1/mean(vecnorm(x,2,dim),"all",'omitnan'));

        % place weights as cross receiver over cells, receiver in "dim"
        W = num2cell(shiftdim(S, -(dim-2)), dim);

        % compute and sum for each cross-receiver
        [z, a, b] = deal(0);
        parfor (i = 1:A, clu)
            w = W{i}; % weight for this set of indices
            xc = sub(x,i,dim); % cross-receiver
            z = z + sum(w .* conj(xc) .* x , [dim, kdim],'omitnan'); % inner-product
            a = a + sum(w .* conj(x ) .* x , [dim, kdim],'omitnan'); % a-norm
            b = b + sum(w .* conj(xc) .* xc, [dim, kdim],'omitnan'); % b-norm
        end

        % normalize by the norm of the vector with 0/0 -> 0
        z = z .* nan2zero(rsqrt(a) .* rsqrt(b));
end
end
