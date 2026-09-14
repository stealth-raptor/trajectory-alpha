function J = compute_itae(t,qRef,q)
%COMPUTE_ITAE Paper ITAE fitness for six-joint tracking.
%
%   Jiang et al., Int. J. Dynam. Control (2025), Eq. (29):
%
%       ITAE = \int_0^\infty t * sum_{j=1}^{6} |e_j(t)| dt
%
%   Discrete trapezoid-equivalent form used here:
%
%       ITAE \approx sum_k t(k) * sum_j |e_j(k)| * Ts

t = t(:);
e = qRef - q;
Ts = t(2) - t(1);
J = sum(t .* sum(abs(e),2)) * Ts;
end
