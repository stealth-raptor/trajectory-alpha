function plant = make_ur5_plant(cfg)
%MAKE_UR5_PLANT  Single definitive paper-faithful Lagrange plant for UR5.
%
%   This is the FINAL plant architecture.  There are NO alternative models,
%   NO rigidBodyTree / universalUR5 path, and NO decoupled approximation.
%
%   The dynamics follow the Lagrange formulation of the reference paper:
%
%       M(q) qdd + C(q,qd) qd + G(q) + tau_f(qd) = tau
%
%   with tau_f,i = fc_i * sign(qd_i) + b_i * qd_i
%
%   All mass, length and inertia values are taken exactly from Table 2.
%   Missing geometric data (DH parameters, COM locations) are supplied by
%   clearly documented assumptions inside paper_lagrange_ur5.m.
%
%   Interface preserved for simulate_ur5_controller.m:
%       plant.accel(q, dq, tau)  returns 6x1 joint accelerations.

% Friction defaults (paper does not publish numerical fc, b values).
% Callers may override via cfg.friction.
if ~isfield(cfg, 'friction')
    cfg.friction.viscous = 0.05 * ones(6,1);
    cfg.friction.coulomb = 0.20 * ones(6,1);
else
    if ~isfield(cfg.friction, 'viscous')
        cfg.friction.viscous = 0.05 * ones(6,1);
    end
    if ~isfield(cfg.friction, 'coulomb')
        cfg.friction.coulomb = 0.20 * ones(6,1);
    end
end

plant.type        = 'paper_lagrange';
plant.description = ['paper-faithful Lagrange coupled 6-DOF UR5 ', ...
                     '(Table 2 params + documented DH/COM assumptions)'];
plant.friction    = cfg.friction;

% Acceleration map required by the existing simulator
plant.accel = @(q, dq, tau) paper_accel(q, dq, tau, cfg.friction);

% Gravity feedforward map used by the controller
plant.gravity = @(q) paper_gravity(q, cfg.friction);

end

% -------------------------------------------------------------------------
function qdd = paper_accel(q, dq, tau, friction)
% Solve  M(q) qdd = tau - C(q,dq)*dq - G(q) - fr(dq)

q   = q(:);
dq  = dq(:);
tau = tau(:);

[M, Cqd, G, fr] = paper_lagrange_ur5(q, dq, friction);

rhs = tau - Cqd - G - fr;

% Guard against numerical singularity (should never happen for a physical robot)
if rcond(M) < 1e-12
    warning('paper_accel:IllConditionedM', ...
        'Inertia matrix is nearly singular (rcond = %.3e).', rcond(M));
end

qdd = M \ rhs;
qdd = qdd(:);

if any(~isfinite(qdd))
    error('paper_accel:NonFinite', ...
        'Non-finite acceleration produced at q = [%s].', sprintf(' %.4g', q));
end
end

% -------------------------------------------------------------------------
function G = paper_gravity(q, friction)
% Return the model gravity torque without recomputing controller dynamics.
[~, ~, G, ~] = paper_lagrange_ur5(q, zeros(6,1), friction);
G = G(:);
end
