function bank = init_oustaloup_bank(lambda,mu,Ts,opts)
%INIT_OUSTALOUP_BANK Fast discrete Oustaloup realization.
% Uses bilinear (Tustin) mapping directly on the first-order factors.
% This keeps the fractional controller initialization lightweight.

nJ = numel(lambda);
bank.I = cell(nJ,1);
bank.D = cell(nJ,1);

for j = 1:nJ
    bank.I{j} = make_filter(-lambda(j),Ts,opts);
    bank.D{j} = make_filter( mu(j),Ts,opts);
end
end

function f = make_filter(alpha,Ts,opts)
N = opts.order;
wL = opts.wLow;
wH = opts.wHigh;
k = -N:N;

wz = wL*(wH/wL).^((k + N + 0.5*(1-alpha))/(2*N+1));
wp = wL*(wH/wL).^((k + N + 0.5*(1+alpha))/(2*N+1));
K = wH^alpha;

% Tustin: s = (2/Ts)*(z-1)/(z+1).
% Each factor (s+wz)/(s+wp) becomes
% (b0 + b1 z^-1)/(1 + a1 z^-1).
c = 2/Ts;
f.b0 = zeros(1,numel(wz));
f.b1 = zeros(1,numel(wz));
f.a1 = zeros(1,numel(wz));
f.u1 = zeros(1,numel(wz));
f.y1 = zeros(1,numel(wz));
f.K = K;

for r = 1:numel(wz)
    f.b0(r) = (c+wz(r))/(c+wp(r));
    f.b1(r) = (wz(r)-c)/(c+wp(r));
    f.a1(r) = (wp(r)-c)/(c+wp(r));
end
end
