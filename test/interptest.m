
I = 64;
T = 32;
N = 4;
M = 3;
F = 2;
t = (0:T-1)'/T;
n = (0:N-1);
m = shiftdim((0:M-1),-2);
f = shiftdim((0:F-1),-1);

wm = 1/2 + (m/(M-1));
wf = 1/2 + (f/(F-1));
x = exp(2j*pi*(1/2+f/2.*n/4).*t);
tau = 2 + (T-4) * randn([T,N,F,1]);
tau = 2 + (T-4) * ((1:I)'./I .* (1+n)./N .* (1+m)./M);

w = 1;
% w = w .* wm;
w = w .* wf;

% typefun = @double;
% typefun = @single;
% typefun = @halfT;
for wvec = {1, wm, wf, wm .* wf} % different outer-product weights
for dsum = [0 1 0 1; 0 0 1 1] % different outer-product sums
for type = string({'single', 'halfT', 'double'}) % different precisions
for i = 1:2 % gpu/cpu
    switch i
        case 1, x = gather(x);
        case 2, x = gpuArray(x);
    end
    [x, tau, w] = dealfun(str2func(type), x, tau, wvec{:});
    y{i,1} = wsinterpd(x, tau, 1, w, 'cubic', 0, 'fsum', dsum(1), 'msum', dsum(2));

%     y_all = cat(y, y_all, y);
%     figure;
%     subplot(2,1,1); plot(real(gather(y(:,:))), '.-');
%     subplot(2,1,2); plot(imag(gather(y(:,:))), '.-');
end
end
end
end