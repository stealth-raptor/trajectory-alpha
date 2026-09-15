%% TUNE_PAPER_BASELINE_MATCHING
% Tunes PID and FOPID gains on the 6-DOF UR5 plant to reproduce the paper's
% exact performance indicators in Table 3 (Step response: OS, Ts, Tp) and
% Table 4 (Sine response: MSE, |tau|).
%
% Paper Target Table 3 (Step Response):
%   PID:   Overshoot = 54.6 %,  Ts = 2.44 s,  Tp = 1.39 s
%   FOPID: Overshoot = 31.2 %,  Ts = 1.89 s,  Tp = 1.33 s
%
% Paper Target Table 4 (Sine Response):
%   PID:   MSE = 0.0227 rad^2,  |tau| = 3.7422e4 Nm
%   FOPID: MSE = 0.0088 rad^2,  |tau| = 2.5686e4 Nm

clear; clc;

%% 1. Configuration
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

% Tuning configuration: faster grid for candidate evaluations
tuneCfg = cfg;
tuneCfg.Ts = 0.005;
tuneCfg.Tfinal = 5.0;

plant = make_ur5_plant(cfg);

% Full-res time and reference vectors
t = (0:cfg.Ts:cfg.Tfinal).';
qStep = zeros(numel(t),6); qStep(t >= cfg.stepTime,:) = cfg.stepAmplitude;
qSine = cfg.sineAmplitude * sin(2*pi*cfg.sineFrequencyHz*t) * ones(1,6);

% Tuning-res time and reference vectors
tTune = (0:tuneCfg.Ts:tuneCfg.Tfinal).';
qStepTune = zeros(numel(tTune),6); qStepTune(tTune >= tuneCfg.stepTime,:) = tuneCfg.stepAmplitude;
qSineTune = tuneCfg.sineAmplitude * sin(2*pi*tuneCfg.sineFrequencyHz*tTune) * ones(1,6);

%% 2. Paper Targets
target.PID.OS   = 54.6;
target.PID.ts   = 2.44;
target.PID.tp   = 1.39;
target.PID.MSE  = 0.0227;
target.PID.tau  = 3.7422e4;

target.FOPID.OS  = 31.2;
target.FOPID.ts  = 1.89;
target.FOPID.tp  = 1.33;
target.FOPID.MSE = 0.0088;
target.FOPID.tau = 2.5686e4;

%% 3. Starting baseline gains
pid.Kp = [140; 120; 100; 60; 40; 25];
pid.Ki = [9.5; 8.0; 7.0; 4.5; 3.0; 1.8];
pid.Kd = [12; 10; 8; 5; 3; 1.5];

fopid = pid;
fopid.Kp = [110; 95; 85; 50; 35; 20];
fopid.Ki = [7.5; 6.5; 5.5; 3.8; 2.2; 1.4];
fopid.Kd = [22; 18; 14; 8; 5; 2.5];
fopid.lambda = 0.85*ones(6,1);
fopid.mu     = 0.85*ones(6,1);

%% 4. Tune PID Controller
fprintf('\n============================================================\n');
fprintf('TUNING PID CONTROLLER TO MATCH PAPER TARGETS\n');
fprintf('Target: OS=%.1f%%, Ts=%.2fs, Tp=%.2fs, MSE=%.4f\n', ...
    target.PID.OS, target.PID.ts, target.PID.tp, target.PID.MSE);
fprintf('============================================================\n');

pid = optimize_controller_matching(plant, tTune, qStepTune, qSineTune, 'PID', pid, tuneCfg, target.PID);

%% 5. Tune FOPID Controller
fprintf('\n============================================================\n');
fprintf('TUNING FOPID CONTROLLER TO MATCH PAPER TARGETS\n');
fprintf('Target: OS=%.1f%%, Ts=%.2fs, Tp=%.2fs, MSE=%.4f\n', ...
    target.FOPID.OS, target.FOPID.ts, target.FOPID.tp, target.FOPID.MSE);
fprintf('============================================================\n');

fopid = optimize_controller_matching(plant, tTune, qStepTune, qSineTune, 'FOPID', fopid, tuneCfg, target.FOPID);

%% 6. Final Full-Resolution Verification (cfg.Ts = 0.002)
fprintf('\n============================================================\n');
fprintf('FINAL FULL-RESOLUTION VERIFICATION (cfg.Ts = 0.002s)\n');
fprintf('============================================================\n');

pidStep   = simulate_ur5_controller(plant, t, qStep, 'PID', pid, cfg);
fopidStep = simulate_ur5_controller(plant, t, qStep, 'FOPID', fopid, cfg);
pidSine   = simulate_ur5_controller(plant, t, qSine, 'PID', pid, cfg);
fopidSine = simulate_ur5_controller(plant, t, qSine, 'FOPID', fopid, cfg);

[stepMetrics, stepSummary] = make_step_metrics(t, qStep, pidStep.q, fopidStep.q, pidStep.tau, fopidStep.tau);
[sineMetrics, sineSummary] = make_sine_metrics(t, qSine, pidSine, fopidSine);

fprintf('\nTABLE 3: STEP RESPONSE (Paper vs Model):\n');
fprintf('Target PID:   OS = %5.1f %%,  Ts = %5.2f s,  Tp = %5.2f s\n', target.PID.OS, target.PID.ts, target.PID.tp);
fprintf('Actual PID:   OS = %5.1f %%,  Ts = %5.2f s,  Tp = %5.2f s\n', ...
    stepSummary.AverageOvershoot_percent(1), stepSummary.AverageSettlingTime_s(1), stepSummary.AveragePeakTime_s(1));
fprintf('Target FOPID: OS = %5.1f %%,  Ts = %5.2f s,  Tp = %5.2f s\n', target.FOPID.OS, target.FOPID.ts, target.FOPID.tp);
fprintf('Actual FOPID: OS = %5.1f %%,  Ts = %5.2f s,  Tp = %5.2f s\n', ...
    stepSummary.AverageOvershoot_percent(2), stepSummary.AverageSettlingTime_s(2), stepSummary.AveragePeakTime_s(2));

fprintf('\nTABLE 4: SINE RESPONSE (Paper vs Model):\n');
fprintf('Target PID:   MSE = %8.4e,  |tau| = %8.4e\n', target.PID.MSE, target.PID.tau);
fprintf('Actual PID:   MSE = %8.4e,  |tau| = %8.4e\n', ...
    sineSummary.AverageMSE_rad2(1), sineSummary.AverageMeanAbsTorque_Nm(1));
fprintf('Target FOPID: MSE = %8.4e,  |tau| = %8.4e\n', target.FOPID.MSE, target.FOPID.tau);
fprintf('Actual FOPID: MSE = %8.4e,  |tau| = %8.4e\n', ...
    sineSummary.AverageMSE_rad2(2), sineSummary.AverageMeanAbsTorque_Nm(2));

save('tuned_pid_fopid_gains_final.mat', 'pid', 'fopid');
fprintf('\nSaved converged gains to tuned_pid_fopid_gains_final.mat\n');

%% =========================================================================
%  OPTIMIZATION HELPER FUNCTIONS
%  =========================================================================

function bestGains = optimize_controller_matching(plant, tTune, qStep, qSine, ctrlType, initGains, cfg, tgt)
bestGains = initGains;
[bestErr, bestM] = evaluate_metrics_cost(plant, tTune, qStep, qSine, ctrlType, bestGains, cfg, tgt);
fprintf('  Initial | Cost = %.4f | OS = %.1f%%, Ts = %.2fs, Tp = %.2fs, MSE = %.4f\n', ...
    bestErr, bestM.OS, bestM.ts, bestM.tp, bestM.MSE);

% Step 1: Global gain scaling sweep (scales entire gain profile up/down)
scales = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.2, 2.8, 3.5];
for s = scales
    cand = bestGains;
    cand.Kp = initGains.Kp * s;
    cand.Ki = initGains.Ki * s;
    cand.Kd = initGains.Kd * s;
    [err, m] = evaluate_metrics_cost(plant, tTune, qStep, qSine, ctrlType, cand, cfg, tgt);
    if err < bestErr
        bestErr = err; bestGains = cand; bestM = m;
        fprintf('  Scale x%.2f | Cost = %.4f | OS = %.1f%%, Ts = %.2fs, Tp = %.2fs, MSE = %.4f\n', ...
            s, bestErr, bestM.OS, bestM.ts, bestM.tp, bestM.MSE);
    end
end

% Step 2: Coordinate pattern search across joints
stepFracs = [0.25, 0.12, 0.05, 0.02];
paramList = {'Kp', 'Kd', 'Ki'};

for sf = stepFracs
    improved = true;
    iter = 0;
    while improved && iter < 6
        improved = false;
        iter = iter + 1;
        for j = 1:6
            for p = 1:numel(paramList)
                pname = paramList{p};
                for mul = [1+sf, 1-sf]
                    cand = bestGains;
                    cand.(pname)(j) = bestGains.(pname)(j) * mul;
                    if strcmp(pname, 'Kp') && bestGains.Kp(j) > 0
                        cand.Ki(j) = cand.Kp(j) * (bestGains.Ki(j)/bestGains.Kp(j));
                    end
                    [err, m] = evaluate_metrics_cost(plant, tTune, qStep, qSine, ctrlType, cand, cfg, tgt);
                    if err < bestErr
                        bestErr = err; bestGains = cand; bestM = m;
                        improved = true;
                        fprintf('  Iter %d j%d %s x%.2f | Cost = %.4f | OS = %.1f%%, Ts = %.2fs, Tp = %.2fs, MSE = %.4f\n', ...
                            iter, j, pname, mul, bestErr, bestM.OS, bestM.ts, bestM.tp, bestM.MSE);
                        break;
                    end
                end
            end
        end
    end
end

% Step 3: For FOPID, sweep fractional orders lambda and mu
if strcmp(ctrlType, 'FOPID')
    orderCandidates = 0.5:0.05:1.2;
    for L = orderCandidates
        for M = orderCandidates
            cand = bestGains;
            cand.lambda = L*ones(6,1);
            cand.mu     = M*ones(6,1);
            [err, m] = evaluate_metrics_cost(plant, tTune, qStep, qSine, ctrlType, cand, cfg, tgt);
            if err < bestErr
                bestErr = err; bestGains = cand; bestM = m;
                fprintf('  Orders L=%.2f M=%.2f | Cost = %.4f | OS = %.1f%%, Ts = %.2fs, Tp = %.2fs, MSE = %.4f\n', ...
                    L, M, bestErr, bestM.OS, bestM.ts, bestM.tp, bestM.MSE);
            end
        end
    end
end

fprintf('  Finished %s Tuning | Final Cost = %.4f | OS = %.1f%%, Ts = %.2fs, Tp = %.2fs, MSE = %.4f\n', ...
    ctrlType, bestErr, bestM.OS, bestM.ts, bestM.tp, bestM.MSE);
end

function [cost, m] = evaluate_metrics_cost(plant, tTune, qStep, qSine, ctrlType, gains, cfg, tgt)
cost = inf;
m.OS = NaN; m.ts = NaN; m.tp = NaN; m.MSE = NaN; m.tau = NaN;

try
    resStep = simulate_ur5_controller(plant, tTune, qStep, ctrlType, gains, cfg);
    if any(~isfinite(resStep.q(:))) || max(abs(resStep.q(:))) > 10*pi
        return;
    end
    
    resSine = simulate_ur5_controller(plant, tTune, qSine, ctrlType, gains, cfg);
    if any(~isfinite(resSine.q(:))) || max(abs(resSine.q(:))) > 10*pi
        return;
    end
    
    [~, ss] = make_step_metrics(tTune, qStep, resStep.q, [], resStep.tau, []);
    [~, ns] = make_sine_metrics(tTune, qSine, resSine, resSine);
    
    m.OS  = ss.AverageOvershoot_percent(1);
    m.ts  = ss.AverageSettlingTime_s(1);
    m.tp  = ss.AveragePeakTime_s(1);
    m.MSE = ns.AverageMSE_rad2(1);
    m.tau = ns.AverageMeanAbsTorque_Nm(1);
    
    if any(~isfinite([m.OS, m.ts, m.tp, m.MSE]))
        return;
    end
    
    % Balanced multi-objective normalized loss
    cost = 1.0 * ((m.OS - tgt.OS) / tgt.OS)^2 + ...
           1.0 * ((m.ts - tgt.ts) / tgt.ts)^2 + ...
           1.0 * ((m.tp - tgt.tp) / tgt.tp)^2 + ...
           2.0 * ((m.MSE - tgt.MSE) / tgt.MSE)^2;
catch
    cost = inf;
end
end
