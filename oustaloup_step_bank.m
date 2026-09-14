function [If,Df,bank] = oustaloup_step_bank(bank,e)
%OUSTALOUP_STEP_BANK Advance six fractional filters by one sample.
% Cascaded first-order Tustin sections are used for speed and numerical
% stability during repeated controller simulations.

nJ = numel(e);
If = zeros(nJ,1);
Df = zeros(nJ,1);

for j = 1:nJ
    [If(j),bank.I{j}] = filter_step(bank.I{j},e(j));
    [Df(j),bank.D{j}] = filter_step(bank.D{j},e(j));
end
end

function [y,f] = filter_step(f,u)
% Cascade the 11 first-order sections.
v = u;
for r = 1:numel(f.b0)
    y = f.b0(r)*v + f.b1(r)*f.u1(r) - f.a1(r)*f.y1(r);
    f.u1(r) = v;
    f.y1(r) = y;
    v = y;
end
y = f.K*v;
end
