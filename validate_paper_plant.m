%% VALIDATE_PAPER_PLANT
% Exhaustive numerical checks that the single Lagrange plant satisfies
% every requirement of the paper-faithful 6-DOF model.
%
% Run from the project root:
%   validate_paper_plant

clear; clc;

fprintf('=== Paper-faithful Lagrange UR5 plant validation ===\n\n');

% Minimal cfg just for friction
cfg.friction.viscous = 0.05*ones(6,1);
cfg.friction.coulomb = 0.20*ones(6,1);

plant = make_ur5_plant(cfg);
fprintf('Plant description: %s\n\n', plant.description);

% ------------------------------------------------------------------
% 1-3. M(q) properties at several configurations
% ------------------------------------------------------------------
test_configs = {
    zeros(6,1), ...
    [0.3; -0.5; 0.8; -0.2; 0.4; -0.1], ...
    [pi/4; pi/6; -pi/3; pi/8; -pi/5; pi/7], ...
    ones(6,1)*0.1
};

all_ok = true;

for ic = 1:numel(test_configs)
    q = test_configs{ic};
    qd = zeros(6,1);                 % pure configuration test
    [M, Cqd, G, fr] = paper_lagrange_ur5(q, qd, cfg.friction);

    fprintf('Config %d: q = [%s]\n', ic, sprintf(' %.3f', q));

    % size
    if ~isequal(size(M), [6 6])
        fprintf('  FAIL: M is not 6x6\n'); all_ok = false;
    else
        fprintf('  PASS: M is 6x6\n');
    end

    % symmetry
    sym_err = norm(M - M.', 'fro') / max(norm(M,'fro'), eps);
    if sym_err > 1e-10
        fprintf('  FAIL: M not symmetric (rel err = %.3e)\n', sym_err); all_ok = false;
    else
        fprintf('  PASS: M approximately symmetric (rel err = %.3e)\n', sym_err);
    end

    % positive definite
    eigM = eig(M);
    if any(eigM <= 0)
        fprintf('  FAIL: M not positive definite (min eig = %.3e)\n', min(eigM)); all_ok = false;
    else
        fprintf('  PASS: M positive definite (min eig = %.3e)\n', min(eigM));
    end

    % coupling: off-diagonal elements
    off_diag = M - diag(diag(M));
    max_off = max(abs(off_diag(:)));
    if max_off < 1e-8
        fprintf('  WARN: M appears diagonal (max off-diag = %.3e) - weak coupling?\n', max_off);
    else
        fprintf('  PASS: M is coupled (max |off-diag| = %.3e)\n', max_off);
    end

    % sizes of other terms
    if ~isequal(size(Cqd), [6 1]) || ~isequal(size(G), [6 1]) || ~isequal(size(fr), [6 1])
        fprintf('  FAIL: Cqd/G/fr size incorrect\n'); all_ok = false;
    else
        fprintf('  PASS: Cqd, G, fr are all 6x1\n');
    end
    fprintf('\n');
end

% ------------------------------------------------------------------
% 4. Explicit G(q) at the zero configuration
% ------------------------------------------------------------------
q0  = zeros(6,1);
qd0 = zeros(6,1);
[~, ~, G0, ~] = paper_lagrange_ur5(q0, qd0, cfg.friction);
fprintf('G(q=0) = \n');
disp(G0);

% At the stretched configuration the gravity torques on the distal
% links should be non-zero for a realistic arm; we only print the value.
if any(abs(G0) > 1e-6)
    fprintf('PASS: Gravity torque is non-zero at the zero pose (as expected for a serial arm).\n\n');
else
    fprintf('NOTE: G(q=0) is numerically zero.  This can occur with the chosen COM mid-point\n');
    fprintf('      assumption when all links are aligned with the base axes; not necessarily an error.\n\n');
end

% Gravity compensation must cancel the model gravity term at rest.
q_gravity = [0.3; -0.5; 0.8; -0.2; 0.4; -0.1];
qdd_gravity = plant.accel(q_gravity, zeros(6,1), plant.gravity(q_gravity));
if norm(qdd_gravity, inf) < 1e-10
    fprintf('PASS: Gravity compensation cancels static acceleration (||qdd||inf = %.3e).\n\n', ...
        norm(qdd_gravity, inf));
else
    fprintf('FAIL: Gravity compensation residual is %.3e.\n\n', ...
        norm(qdd_gravity, inf));
    all_ok = false;
end

% ------------------------------------------------------------------
% 5. Coupling demonstration: change of one joint affects another joint's accel
% ------------------------------------------------------------------
q_base = [0.2; 0.3; -0.4; 0.1; -0.2; 0.05];
qd_base = [0.1; -0.15; 0.05; 0; 0.02; -0.01];
tau_base = zeros(6,1);

qdd1 = plant.accel(q_base, qd_base, tau_base);

% Perturb only joint 2 position
q_pert = q_base;
q_pert(2) = q_pert(2) + 0.05;
qdd2 = plant.accel(q_pert, qd_base, tau_base);

delta_qdd = qdd2 - qdd1;
fprintf('Coupling test (perturb q2 by +0.05 rad):\n');
fprintf('  Delta qdd = [%s]\n', sprintf(' %.4e', delta_qdd));
if any(abs(delta_qdd([1 3 4 5 6])) > 1e-6)
    fprintf('  PASS: Acceleration of other joints changed -> dynamics are coupled.\n\n');
else
    fprintf('  FAIL: No visible coupling.\n\n'); all_ok = false;
end

% ------------------------------------------------------------------
% 6. Confirm no universalUR5 / rigidBodyTree dependency
% ------------------------------------------------------------------
% Strip comments before searching so documentation that mentions the
% *absence* of these symbols does not trigger a false positive.
src  = regexprep(fileread('make_ur5_plant.m'), '%.*$', '', 'lineanchors');
src2 = regexprep(fileread('paper_lagrange_ur5.m'), '%.*$', '', 'lineanchors');
if contains(src, 'loadrobot') || contains(src, 'universalUR5') || ...
   contains(src, 'rigidBodyTree') || contains(src2, 'loadrobot') || ...
   contains(src2, 'universalUR5') || contains(src2, 'rigidBodyTree')
    fprintf('FAIL: Forbidden toolbox symbols still present in plant files.\n');
    all_ok = false;
else
    fprintf('PASS: No universalUR5 / rigidBodyTree / loadrobot references remain.\n\n');
end

% ------------------------------------------------------------------
% 7. Friction form
% ------------------------------------------------------------------
qd_test = [0.5; -0.3; 0; 0.1; -0.2; 0.4];
[~,~,~,fr_test] = paper_lagrange_ur5(zeros(6,1), qd_test, cfg.friction);
expected_fr = cfg.friction.coulomb .* sign(qd_test) + cfg.friction.viscous .* qd_test;
if norm(fr_test - expected_fr) < 1e-12
    fprintf('PASS: Friction exactly matches fc*sign(qd) + b*qd.\n\n');
else
    fprintf('FAIL: Friction form incorrect.\n\n'); all_ok = false;
end

% ------------------------------------------------------------------
% Final verdict
% ------------------------------------------------------------------
if all_ok
    fprintf('=== ALL VALIDATION CHECKS PASSED ===\n');
    fprintf('The plant is ready to be frozen for subsequent PID/FOPID work.\n');
else
    fprintf('=== SOME CHECKS FAILED - INSPECT OUTPUT ABOVE ===\n');
end