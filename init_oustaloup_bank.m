function bank = init_oustaloup_bank(lambda,mu,Ts,opts)
%INIT_OUSTALOUP_BANK Fast discrete Oustaloup realization.
% Uses bilinear (Tustin) mapping directly on the first-order factors.
% This keeps the fractional controller initialization lightweight.

nJ = numel(lambda);
N = opts.order;
nSec = 2*N + 1;

bank.I_b0 = zeros(nJ, nSec);
bank.I_b1 = zeros(nJ, nSec);
bank.I_a1 = zeros(nJ, nSec);
bank.I_u1 = zeros(nJ, nSec);
bank.I_y1 = zeros(nJ, nSec);
bank.I_K  = zeros(nJ, 1);

bank.D_b0 = zeros(nJ, nSec);
bank.D_b1 = zeros(nJ, nSec);
bank.D_a1 = zeros(nJ, nSec);
bank.D_u1 = zeros(nJ, nSec);
bank.D_y1 = zeros(nJ, nSec);
bank.D_K  = zeros(nJ, 1);

for j = 1:nJ
    fI = make_filter(-lambda(j),Ts,opts);
    bank.I_b0(j,:) = fI.b0;
    bank.I_b1(j,:) = fI.b1;
    bank.I_a1(j,:) = fI.a1;
    bank.I_K(j)    = fI.K;

    fD = make_filter(mu(j),Ts,opts);
    bank.D_b0(j,:) = fD.b0;
    bank.D_b1(j,:) = fD.b1;
    bank.D_a1(j,:) = fD.a1;
    bank.D_K(j)    = fD.K;
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
