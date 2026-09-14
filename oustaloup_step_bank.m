function [If,Df,bank] = oustaloup_step_bank(bank,e)
%OUSTALOUP_STEP_BANK Advance six fractional filters by one sample.
% Cascaded first-order Tustin sections are used for speed and numerical
% stability during repeated controller simulations.

vI = e(:);
nSec = size(bank.I_b0, 2);
for r = 1:nSec
    yI = bank.I_b0(:,r).*vI + bank.I_b1(:,r).*bank.I_u1(:,r) - bank.I_a1(:,r).*bank.I_y1(:,r);
    bank.I_u1(:,r) = vI;
    bank.I_y1(:,r) = yI;
    vI = yI;
end
If = bank.I_K .* vI;

vD = e(:);
for r = 1:nSec
    yD = bank.D_b0(:,r).*vD + bank.D_b1(:,r).*bank.D_u1(:,r) - bank.D_a1(:,r).*bank.D_y1(:,r);
    bank.D_u1(:,r) = vD;
    bank.D_y1(:,r) = yD;
    vD = yD;
end
Df = bank.D_K .* vD;
end
