function [M, Cqd, G, fr] = paper_lagrange_ur5(q, qd, friction)
%PAPER_LAGRANGE_UR5  Paper-faithful coupled 6-DOF Lagrange dynamics for UR5.
%
%   [M, Cqd, G, fr] = paper_lagrange_ur5(q, qd, friction)
%
%   Implements the dynamic equation required by the reference paper:
%
%       M(q) * qdd + C(q,qd)*qd + G(q) + tau_f(qd) = tau
%
%   where tau_f,i = fc_i * sign(qd_i) + b_i * qd_i
%
%   The kinetic and potential energies follow the Lagrange derivation
%   outlined in the paper (homogeneous transforms, link inertias, masses).
%   Because the paper does not publish a complete DH table or COM locations,
%   a documented standard UR5-like kinematic convention is used for the
%   missing geometric data. All numerical mass / length / inertia values
%   are taken EXACTLY from Table 2 of the paper.
%
%   INPUTS
%     q         6x1 joint positions [rad]
%     qd        6x1 joint velocities [rad/s]
%     friction  struct with fields
%                 .viscous  6x1  (b_i)
%                 .coulomb  6x1  (fc_i)
%
%   OUTPUTS
%     M         6x6 configuration-dependent inertia matrix
%     Cqd       6x1 Coriolis/centrifugal vector  (= C(q,qd)*qd)
%     G         6x1 gravity vector
%     fr        6x1 friction torque vector
%
%   Geometric / modelling assumptions (NOT stated in the paper)
%   ----------------------------------------------------------
%   - Classic DH parameters with paper lengths mapped to a/d:
%       a     = [0, 0.612, 0.571, 0, 0, 0]
%       d     = [0.128, 0, 0, 0.164, 0.115, 0.092]
%       alpha = [pi/2, 0, 0, pi/2, -pi/2, 0]
%   - Centre-of-mass of each link expressed in its own link frame,
%     placed at the geometric mid-point of the corresponding length
%     parameter (standard approximation when only a scalar length is given):
%       Link 1: [0; 0; d1/2]
%       Link 2: [a2/2; 0; 0]
%       Link 3: [a3/2; 0; 0]
%       Link 4: [0; 0; d4/2]
%       Link 5: [0; 0; d5/2]
%       Link 6: [0; 0; d6/2]
%   - The diagonal inertia matrices of Table 2 are interpreted as the
%     inertia tensors about the COM, expressed in principal axes that
%     coincide with the link frame axes.
%   - Gravity vector in the base frame: g = [0; 0; -9.81] m/s^2
%
%   These geometric choices are isolated here so that any future audit
%   can replace them without touching the paper-supplied numerical data.

q  = q(:);
qd = qd(:);
assert(numel(q)==6 && numel(qd)==6, 'q and qd must be 6x1');

%% ---------------------------------------------------------------------
%  Paper-supplied numerical parameters (Table 2) - DO NOT CHANGE
%  ---------------------------------------------------------------------
m = [2.0; 2.5; 5.7; 3.9; 2.5; 2.5];                 % kg

% Diagonal inertia tensors about COM (kg·m²) exactly as published
I_com = cell(6,1);
I_com{1} = diag([1.0, 1.0, 1.0]) * 1e-3;
I_com{2} = diag([4.0, 4.0, 4.0]) * 1e-2;
I_com{3} = diag([6.0, 6.0, 4.0]) * 1e-2;
I_com{4} = diag([5.5, 5.5, 4.0]) * 1e-2;
I_com{5} = diag([4.0, 4.0, 4.0]) * 1e-2;
I_com{6} = diag([4.0, 4.0, 4.0]) * 1e-2;

%% ---------------------------------------------------------------------
%  Documented geometric assumptions (paper does not provide these)
%  ---------------------------------------------------------------------
a     = [0;     0.612; 0.571; 0;     0;     0    ];
d     = [0.128; 0;     0;     0.164; 0.115; 0.092];
alpha = [pi/2;  0;     0;     pi/2; -pi/2;  0    ];

% COM locations in each link frame (mid-point of the length parameter)
com_rel = [ 0,        a(2)/2,  a(3)/2,  0,        0,        0;
            0,        0,       0,       0,        0,        0;
            d(1)/2,   0,       0,       d(4)/2,   d(5)/2,   d(6)/2 ];

g = [0; 0; -9.81];   % base-frame gravity

%% ---------------------------------------------------------------------
%  Forward kinematics + geometric Jacobians of the COMs
%  ---------------------------------------------------------------------
[T, z_axes, origins] = classic_dh_fk(q, a, d, alpha);

% Pre-allocate
Jv = cell(6,1);
Jw = cell(6,1);
R  = cell(6,1);
p_com = zeros(3,6);

for i = 1:6
    Ri = T{i}(1:3,1:3);
    R{i} = Ri;
    p_com(:,i) = Ri * com_rel(:,i) + T{i}(1:3,4);

    Jv{i} = zeros(3,6);
    Jw{i} = zeros(3,6);
    for j = 1:i
        zj = z_axes(:,j);                 % z_{j-1} in base
        oj = origins(:,j);                % origin of frame j-1
        Jw{i}(:,j) = zj;
        Jv{i}(:,j) = cross(zj, p_com(:,i) - oj);
    end
end

%% ---------------------------------------------------------------------
%  Inertia matrix M(q)  (fully coupled)
%  ---------------------------------------------------------------------
M = zeros(6,6);
for i = 1:6
    Jvi = Jv{i};
    Jwi = Jw{i};
    Ii_world = R{i} * I_com{i} * R{i}.';
    M = M + m(i) * (Jvi.' * Jvi) + (Jwi.' * Ii_world * Jwi);
end
% Enforce symmetry (numerical round-off)
M = 0.5 * (M + M.');

%% ---------------------------------------------------------------------
%  Gravity vector G(q)
%  ---------------------------------------------------------------------
G = zeros(6,1);
for i = 1:6
    G = G + m(i) * (Jv{i}.' * g);
end

%% ---------------------------------------------------------------------
%  Coriolis / centrifugal vector C(q,qd)*qd
%  Computed via Christoffel symbols with central finite differences
%  of M (n=6, inexpensive and energy-consistent to first order).
%  ---------------------------------------------------------------------
Cqd = christoffel_Cqd(q, qd, M, a, d, alpha, m, I_com, com_rel);

%% ---------------------------------------------------------------------
%  Friction (exact structural form required by the paper)
%  ---------------------------------------------------------------------
fc = friction.coulomb(:);
b  = friction.viscous(:);
fr = fc .* sign(qd) + b .* qd;

end

% =====================================================================
%  Local helpers
% =====================================================================

function [Tcell, z_axes, origins] = classic_dh_fk(q, a, d, alpha)
% Classic DH forward kinematics.
% Returns:
%   Tcell{i}  = 0T_i   (4x4)
%   z_axes    = [z0 ... z5]  (3x6)   z_{j-1} for joint j
%   origins   = [o0 ... o5]  (3x6)   origin of frame j-1

Tcell   = cell(6,1);
z_axes  = zeros(3,6);
origins = zeros(3,6);

T = eye(4);
z_axes(:,1)  = [0;0;1];          % z0
origins(:,1) = [0;0;0];          % o0

for i = 1:6
    ct = cos(q(i));  st = sin(q(i));
    ca = cos(alpha(i)); sa = sin(alpha(i));

    A = [ ct, -st*ca,  st*sa, a(i)*ct;
          st,  ct*ca, -ct*sa, a(i)*st;
           0,     sa,     ca,    d(i);
           0,      0,      0,      1 ];

    T = T * A;
    Tcell{i} = T;

    if i < 6
        z_axes(:,i+1)  = T(1:3,3);
        origins(:,i+1) = T(1:3,4);
    end
end
end

function Cqd = christoffel_Cqd(q, qd, M0, a, d, alpha, m, I_com, com_rel)
% Assemble C*qd from Christoffel symbols using vector tensor formulation:
%   (C*qd)_i = (Mdot*qd)_i - 0.5 * (qd.' * (dM/dq_i) * qd)
% where Mdot = sum_k (dM/dq_k * qd_k).

n   = 6;
eps = 1e-7;
dM = cell(n,1);
Mdot = zeros(n,n);

for k = 1:n
    qp = q; qp(k) = qp(k) + eps;
    qm = q; qm(k) = qm(k) - eps;
    Mp = inertia_only(qp, a, d, alpha, m, I_com, com_rel);
    Mm = inertia_only(qm, a, d, alpha, m, I_com, com_rel);
    dM{k} = (Mp - Mm) / (2*eps);
    Mdot = Mdot + dM{k} * qd(k);
end

Cqd = Mdot * qd;
for i = 1:n
    Cqd(i) = Cqd(i) - 0.5 * (qd.' * dM{i} * qd);
end
end

function M = inertia_only(q, a, d, alpha, m, I_com, com_rel)
% Lightweight evaluation of M(q) only (used by finite-difference).
[T, z_axes, origins] = classic_dh_fk(q, a, d, alpha);
M = zeros(6,6);
for i = 1:6
    Ri = T{i}(1:3,1:3);
    p  = Ri * com_rel(:,i) + T{i}(1:3,4);
    Jv = zeros(3,6);
    Jw = zeros(3,6);
    for j = 1:i
        zj = z_axes(:,j);
        oj = origins(:,j);
        Jw(:,j) = zj;
        Jv(:,j) = cross(zj, p - oj);
    end
    Ii_w = Ri * I_com{i} * Ri.';
    M = M + m(i)*(Jv.'*Jv) + (Jw.'*Ii_w*Jw);
end
M = 0.5*(M + M.');
end
