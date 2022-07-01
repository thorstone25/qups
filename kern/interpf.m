function y = interpf(x, t, dim, w, sdim)
% INTERPF Frequency-domain interpolation in one dimension
%
% y = INTERPF(x, t) interpolates the data x at the indices t by expressing
% the data in the frequency domain and summing over all frequencies at the
% phase corresponding to t.
%
% x is 0-base indexed implicitly i.e. x(1) in MATLAB <-> t == 0
%
% y = INTERPF(x, t, dim) samples the data in dimension dim. The default 
% is 1.
% 
% y = INTERPF(x, t, dim, w) applies the weighting array w via point-wise 
% multiplication after sampling the data. The dimensions must be 
% compatible with the data size. The default is 1.
%  
% y = INTERPF(x, t, dim, w, sdim) sums the data in the dimension(s) sdim
% after the weighting matrix has been applied. The default is [].
%
% See also INTERPN INTERP1 INTERPD WSINTERPD

% extend data if necessary
% ntau = (tau - chd.t0) * chd.fs; % sample delays (I x [1|N] x [1|M] x [1|F] x ...)
% nwrap = min(0,floor(min(ntau,[],'all'))); % amount the signal is wrapped around from negative time
% next  = max(0, ceil(max(ntau,[],'all')) - chd.T-1); % amount the signal is extended in positive time
% chd = zeropad(chd, -nwrap, next); % extend signal to ensure we sample zeros
% x = chd.data; % reference data (T x N x M x [1|F] x ...)

% clear sdim of singleton dimensions
if ~isempty(sdim), sdim = sdim(size(x, sdim) ~= 1 | size(t, sdim) ~= 1); end

L = gather(size(x, dim)); % data/fft length
% l = (0:L-1)'; % make new time vector in sample dimension
d = max(ndims(x), ndims(t)); % find max dimension
t = swapdim(t, d+1, dim); % move sampling to a free dimension
w = swapdim(w, d+1, dim); % move sampling to a free dimension

% apply phase shifts and sum (code written serially to
% request in-place operation from MATLAB)
x = fft(x, L, dim); % put data in freq domain (L x N x M x [1|F] x ... x I)
kL = complex(cospi(2*t./L),sinpi(2*t./L)); % sampling steering vector (1 x [1|N] x [1|M] x [1|F] x ... x I)
y = 0; % initialize accumulator
if isa(x, 'gpuArray') || isa(kL, 'gpuArray'), clu = 0; % avoid parfor on gpuArray
else, clu = gcp('nocreate'); if isempty(clu), clu = 0; end % otherwise, use current parpool
end

% apply phase shift and sum over freq (1 x N x M x [1|F] x ... x I)
if istall(x) || istall(kL)
    rdim = [dim, sdim]; % reduction dimensions
    l = shiftdim((1:L)', 1-dim); % each frequency index, in dim tdim
    y = matlab.tall.reduce(@(x,k,l) w .* sum(k.^(l-1) .* x, rdim), @(x)x, x, kL, l); % apply via map-reduce
else
    rdim = [d+2, sdim]; % add a singleton dimension in case sdim is empty
    xl = num2cell(x, setdiff(1:ndims(x), dim)); % splice data in freq domain (L x {N x M x [1|F] x ... x I})
    parfor (l = (1:L), clu), y = y + sum(w .* kL.^(l-1) .* xl{l}, rdim); end % apply, 1 freq at a time
end
y = swapdim(y, dim, d+1); % move samples back to first dim (I x N x 1 x M' x F x ...)
