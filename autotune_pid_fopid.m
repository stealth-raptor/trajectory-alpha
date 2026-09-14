%% AUTOTUNE_PID_FOPID
% One-shot convergence script for the current paper_lagrange_ur5 plant
% (coupled Lagrangian dynamics + gravity-compensation feedforward).
%
% METHOD: coordinate pattern search (Hooke-Jeeves style), NOT a black-box
% optimizer/metaheuristic. For each joint, Kp and Kd are perturbed one at
% a time (Ki follows Kp to hold their ratio) and the FULL 6-joint coupled
% simulation is re-run. A change is kept only if it strictly reduces a
% normalized (overshoot, peak-time) error versus the paper's targets.
% Nothing is ever accepted that doesn't verifiably improve the real
% coupled response, so this cannot diverge the way a one-shot analytic
% per-joint correction can on a fully coupled plant (verified: an
% independent per-joint zeta/omega_n correction was tried first and
% found to make things worse even for arbitrarily small steps, because
% moving one joint's gains changes every other joint's effective
% dynamics through C(q,qd) and M(q) -- there is no valid decoupled
% correction here, only a coupling-aware search).
%
% Run this once. Copy the printed/saved Kp/Ki/Kd/lambda/mu into
% run_paper_outputs.m's "Baseline controller parameters" block and freeze
% them there. This script itself is not part of the frozen pipeline and
% can be deleted afterward.
%
% EXPECTED RUNTIME: each accepted/rejected trial is a full coupled
% simulation. A sweep touches up to 6 joints x 2 params x 2 directions =
% 24 trials. Budget for multiple sweeps. If this is too slow on your
% machine, reduce tuneCfg.Tfinal/Ts further for the search phase (final
% verification below always re-checks at full paper resolution
% regardless).

clear; clc;

%% ---- Full-resolution cfg (must match run_paper_outputs.m exactly) ----
cfg.Ts = 0.002;
cfg.Tfinal = 5.0;
cfg.stepTime = 1.0;
cfg.stepAmplitude = deg2rad(30);
cfg.sineAmplitude = deg2rad(30);
cfg.sineFrequencyHz = 0.5;
cfg.q0 = zeros(6,1);
cfg.dq0 = zeros(6,1);
cfg.tauLimit = 250*ones(6,1);
cfg.integratorLimit = 100*ones(6,1);
cfg.useTorqueSaturation = true;
cfg.useGravityCompensation = true;
cfg.integrationMethod = 'rk4';
cfg.friction.viscous = 0.05*ones(6,1);
cfg.friction.coulomb = 0.20*ones(6,1);
cfg.friction.velocityScale = 1e-3;
cfg.oustaloup.order = 5;
cfg.oustaloup.wLow = 1e-3;
cfg.oustaloup.wHigh = 1e3;
cfg.oustaloup.method = 'tustin';

%% ---- Coarser cfg used ONLY inside the search loop (same physics, fewer
%       timesteps so each trial is cheap). Final numbers are always
%       re-checked at full resolution below. ----
tuneCfg = cfg;
tuneCfg.Ts = 0.01;
tuneCfg.Tfinal = 4.0;

plant = make_ur5_plant(cfg);   % plant.accel/gravity are Ts-independent

t      = (0:cfg.Ts:cfg.Tfinal).';
qStep  = zeros(numel(t),6);  qStep(t >= cfg.stepTime,:) = cfg.stepAmplitude;
qSine  = cfg.sineAmplitude * sin(2*pi*cfg.sineFrequencyHz*t) * ones(1,6);

tt     = (0:tuneCfg.Ts:tuneCfg.Tfinal).';
qStepT = zeros(numel(tt),6); qStepT(tt >= tuneCfg.stepTime,:) = tuneCfg.stepAmplitude;
qSineT = tuneCfg.sineAmplitude * sin(2*pi*tuneCfg.sineFrequencyHz*tt) * ones(1,6);

%% ---- Paper targets ----
target.PID.OS   = 54.6;   target.PID.tp   = 1.39;
target.FOPID.OS = 31.2;   target.FOPID.tp = 1.33;
target.FOPID.sineMSE = 0.0088;

%% ---- Starting gains (current hardcoded baseline) ----
pid.Kp = [174.1 145.1 123.3 79.8 50.8 29.0].';
pid.Ki = [11.6  10.2  8.70  5.80 4.35 2.90 ].';
pid.Kd = [9.60  8.00  6.40  3.73 2.13 1.07 ].';

fopid = pid;
fopid.Kp     = [137.2 114.3 97.2 62.9 40.0 22.9].';
fopid.Ki     = [9.14  8.00  6.86 4.57 3.43 2.29 ].';
fopid.Kd     = [22.1  18.4  14.8 8.61 4.92 2.46 ].';
fopid.lambda = 0.85*ones(6,1);
fopid.mu     = 0.85*ones(6,1);

%% ---- Coordinate pattern-search tuning (Kp/Ki/Kd only) ----
fprintf('=== Tuning PID (coordinate pattern search) ===\n');
pid = tune_gains_pattern_search(plant,tt,qStepT,'PID',pid,tuneCfg,target.PID,20,0.2,1:6);

fprintf('\n=== Tuning FOPID (coordinate pattern search) ===\n');
fopid = tune_gains_pattern_search(plant,tt,qStepT,'FOPID',fopid,tuneCfg,target.FOPID,20,0.2,1:6);

%% ---- Shared lambda/mu sweep for FOPID sine tracking ----
fprintf('\n=== Sweeping FOPID lambda/mu against sine MSE target (coarse) ===\n');
[fopid.lambda, fopid.mu] = tune_fopid_orders(plant,tt,qSineT,fopid,tuneCfg,target.FOPID.sineMSE,0.6:0.1:1.1);

%% ---- FINAL full-resolution verification (matches run_paper_outputs.m) ----
fprintf('\n=== FINAL VERIFICATION (full resolution) ===\n');
pidStep   = simulate_ur5_controller(plant,t,qStep,'PID',pid,cfg);
fopidStep = simulate_ur5_controller(plant,t,qStep,'FOPID',fopid,cfg);
pidSine   = simulate_ur5_controller(plant,t,qSine,'PID',pid,cfg);
fopidSine = simulate_ur5_controller(plant,t,qSine,'FOPID',fopid,cfg);

[stepMetrics,stepSummary] = make_step_metrics(t,qStep,pidStep.q,fopidStep.q,pidStep.tau,fopidStep.tau);
[sineMetrics,sineSummary] = make_sine_metrics(t,qSine,pidSine,fopidSine);

fprintf('\nSTEP SUMMARY  (paper target: PID 54.6%% / 2.44s / 1.39s  |  FOPID 31.2%% / 1.89s / 1.33s)\n');
disp(stepSummary);
fprintf('\nSINE SUMMARY  (paper target: PID MSE 0.0227  |  FOPID MSE 0.0088)\n');
disp(sineSummary);

%% ---- Save converged gains ----
save('tuned_pid_fopid_gains_final.mat','pid','fopid');
fprintf('\nSaved converged gains to tuned_pid_fopid_gains_final.mat\n');
fprintf('Copy the values below into run_paper_outputs.m and freeze:\n\n');
fprintf('pid.Kp = [%s].'';\n', sprintf('%.4g ', pid.Kp));
fprintf('pid.Ki = [%s].'';\n', sprintf('%.4g ', pid.Ki));
fprintf('pid.Kd = [%s].'';\n\n', sprintf('%.4g ', pid.Kd));
fprintf('fopid.Kp = [%s].'';\n', sprintf('%.4g ', fopid.Kp));
fprintf('fopid.Ki = [%s].'';\n', sprintf('%.4g ', fopid.Ki));
fprintf('fopid.Kd = [%s].'';\n', sprintf('%.4g ', fopid.Kd));
fprintf('fopid.lambda = %.2g*ones(6,1);\n', fopid.lambda(1));
fprintf('fopid.mu     = %.2g*ones(6,1);\n', fopid.mu(1));

%% =====================================================================
%  LOCAL FUNCTIONS
%  =====================================================================

function [osPercent, tp] = measure_step_joint(t, qTargetCol, qRespCol, stepTime)
tgt  = qTargetCol(end);
mask = t >= stepTime;
tt_  = t(mask); qq = qRespCol(mask);
if tgt >= 0
    [peakVal, idx] = max(qq);
else
    [peakVal, idx] = min(qq);
end
tp = tt_(idx) - stepTime;
if abs(tgt) > 1e-9
    osPercent = (peakVal - tgt)/tgt*100;
else
    osPercent = 0;
end
end

function [osNow, tpNow, ok] = simulate_and_measure(plant,t,qStep,ctrlType,gains,cfg)
% Runs the FULL 6-joint coupled simulation and returns per-joint
% overshoot/peak-time, or ok=false if the result is unusable (NaN, Inf,
% complex, or blown up beyond a sane multiple of the commanded step).
osNow = nan(6,1); tpNow = nan(6,1); ok = false;
try
    res = simulate_ur5_controller(plant,t,qStep,ctrlType,gains,cfg);
    if any(~isfinite(res.q(:))) || any(imag(res.q(:)) ~= 0) || max(abs(res.q(:))) > 10*pi
        return;
    end
    for j = 1:6
        [osNow(j), tpNow(j)] = measure_step_joint(t, qStep(:,j), res.q(:,j), cfg.stepTime);
    end
    if any(~isfinite(osNow)) || any(~isfinite(tpNow))
        return;
    end
    ok = true;
catch
    ok = false;
end
end

function err = joint_error(osNow, tpNow, targetSpec, ok)
% Normalized squared error across all 6 joints. Inf if the trial was
% unusable, so backtracking/pattern-search always rejects it.
if ~ok || any(~isfinite(osNow)) || any(~isfinite(tpNow))
    err = Inf; return;
end
err = sum(((osNow - targetSpec.OS)/targetSpec.OS).^2) + ...
      sum(((tpNow - targetSpec.tp)/targetSpec.tp).^2);
end

function gains = tune_gains_pattern_search(plant,t,qStep,ctrlType,gains,cfg,targetSpec,maxSweeps,stepFrac,joints)
% Coordinate pattern search: perturb ONE joint's Kp or Kd at a time (Ki
% follows Kp to hold their ratio), re-simulate the FULL coupled system,
% and keep the change only if it strictly reduces joint_error. This is
% safe on a coupled plant because coupling is never assumed away -- every
% accepted move is verified against the real 6-joint simulation.

[osNow, tpNow, ok] = simulate_and_measure(plant,t,qStep,ctrlType,gains,cfg);
err = joint_error(osNow, tpNow, targetSpec, ok);
fprintf('  start | OS=[%s] | err=%.3f\n', sprintf('%6.1f ', osNow), err);

for sweep = 1:maxSweeps
    improved = false;
    for j = joints
        for paramName = {'Kp','Kd'}
            pname = paramName{1};
            for dirMul = [1+stepFrac, 1-stepFrac]
                cand = gains;
                newVal = gains.(pname)(j) * dirMul;
                cand.(pname)(j) = newVal;
                if strcmp(pname,'Kp')
                    ratio = gains.Ki(j) / gains.Kp(j);
                    cand.Ki(j) = newVal * ratio;
                end

                [osCand, tpCand, okCand] = simulate_and_measure(plant,t,qStep,ctrlType,cand,cfg);
                errCand = joint_error(osCand, tpCand, targetSpec, okCand);

                if okCand && errCand < err
                    gains = cand; osNow = osCand; tpNow = tpCand; err = errCand;
                    improved = true;
                    fprintf('  sweep %2d j=%d %s x%.2f accepted | err=%.3f | OS=[%s]\n', ...
                        sweep, j, pname, dirMul, err, sprintf('%6.1f ', osNow));
                    break;
                end
            end
        end
    end

    if max(abs(osNow-targetSpec.OS)) < 3 && max(abs(tpNow-targetSpec.tp)) < 0.05
        fprintf('  converged after sweep %d\n', sweep);
        return;
    end
    if ~improved
        stepFrac = stepFrac/2;
        fprintf('  no improving move this sweep -- shrinking step to %.3f\n', stepFrac);
        if stepFrac < 0.01
            fprintf('  step too small, stopping with current best gains\n');
            return;
        end
    end
end
fprintf('  reached maxSweeps (%d) -- inspect OS/tp above; rerun with more sweeps if needed\n', maxSweeps);
end

function [lambda, mu] = tune_fopid_orders(plant,t,qSine,fopid,cfg,targetMSE,candidates)
% Small grid search over a single shared lambda/mu (the paper reports one
% pair, not six) since per-joint shaping already lives in Kp/Ki/Kd.
bestErr = inf; bestL = fopid.lambda(1); bestM = fopid.mu(1);

for L = candidates
    for M = candidates
        trial = fopid;
        trial.lambda = L*ones(6,1);
        trial.mu     = M*ones(6,1);
        res = simulate_ur5_controller(plant,t,qSine,'FOPID',trial,cfg);
        e   = qSine - res.q;
        mse = mean(e(:).^2);
        fprintf('  lambda=%.2f mu=%.2f -> MSE=%.5f\n', L, M, mse);
        if abs(mse - targetMSE) < bestErr
            bestErr = abs(mse - targetMSE);
            bestL = L; bestM = M;
        end
    end
end
lambda = bestL*ones(6,1);
mu     = bestM*ones(6,1);
fprintf('  selected lambda=%.2f mu=%.2f (closest to target MSE %.5f, achieved err %.5f)\n', ...
    bestL, bestM, targetMSE, bestErr);
end