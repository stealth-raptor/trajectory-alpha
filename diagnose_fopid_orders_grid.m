%% DIAGNOSE_FOPID_ORDERS_GRID
% Diagnostic-only grid search over FOPID lambda and mu.
% Kp, Ki and Kd remain fixed at the current project baseline values.
% This script does not modify run_paper_outputs.m or any baseline file.

clear; clc; close all;

%% Current simulation configuration
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
cfg.integrationMethod = 'rk4';

cfg.friction.viscous = 0.05*ones(6,1);
cfg.friction.coulomb = 0.20*ones(6,1);
cfg.friction.velocityScale = 1e-3;

cfg.oustaloup.order = 5;
cfg.oustaloup.wLow = 1e-3;
cfg.oustaloup.wHigh = 1e3;
cfg.oustaloup.method = 'tustin';

%% Fixed FOPID gains
% Do not change these values in this diagnostic.
fopid.Kp = [137.2 114.3 97.2 62.9 40.0 22.9].';
fopid.Ki = [9.14 8.00 6.86 4.57 3.43 2.29].';
fopid.Kd = [22.1 18.4 14.8 8.61 4.92 2.46].';

%% Grid and step-response reference targets
lambdaValues = [0.6 0.7 0.8 0.9 1.0 1.1];
muValues = [0.6 0.7 0.8 0.9 1.0 1.1];

baselineStep.Overshoot_percent = 33.03;
baselineStep.SettlingTime_s = 1.828;
baselineStep.PeakTime_s = 1.368;

% The sine MSE is the primary objective. Step deviations are normalized by
% practical tolerances and contribute a smaller secondary penalty.
stepTolerance.Overshoot_percent = 10.0;
stepTolerance.SettlingTime_s = 0.50;
stepTolerance.PeakTime_s = 0.30;
stepPenaltyWeight = 0.25;

%% Plant and references
plant = make_ur5_plant(cfg);
t = (0:cfg.Ts:cfg.Tfinal).';
qStep = zeros(numel(t),6);
qStep(t >= cfg.stepTime,:) = cfg.stepAmplitude;
qSine = cfg.sineAmplitude * ...
    sin(2*pi*cfg.sineFrequencyHz*t) * ones(1,6);

fprintf('\nFOPID lambda/mu diagnostic grid search\n');
fprintf('Plant: %s\n',plant.description);
fprintf('Kp, Ki and Kd are fixed for all 36 trials.\n');

%% Allocate result table
nTrials = numel(lambdaValues)*numel(muValues);
rows = cell(nTrials,10);
trial = 0;

for lambda = lambdaValues
    for mu = muValues
        trial = trial + 1;
        gains = fopid;
        gains.lambda = lambda*ones(6,1);
        gains.mu = mu*ones(6,1);

        fprintf('Trial %2d/%2d: lambda = %.1f, mu = %.1f\n', ...
            trial,nTrials,lambda,mu);

        sineResult = simulate_ur5_controller( ...
            plant,t,qSine,'FOPID',gains,cfg);
        stepResult = simulate_ur5_controller( ...
            plant,t,qStep,'FOPID',gains,cfg);

        eSine = qSine-sineResult.q;
        sineMSE = mean(eSine.^2,'all');
        sineRMSE = sqrt(sineMSE);
        sineTorque = mean(sum(abs(sineResult.tau),1));

        [~,stepSummary] = make_step_metrics( ...
            t,qStep,[],stepResult.q,[],stepResult.tau);
        stepOvershoot = stepSummary.AverageOvershoot_percent(1);
        stepSettling = stepSummary.AverageSettlingTime_s(1);
        stepPeak = stepSummary.AveragePeakTime_s(1);

        stepPenalty = ...
            ((stepOvershoot-baselineStep.Overshoot_percent)/ ...
                stepTolerance.Overshoot_percent)^2 + ...
            ((stepSettling-baselineStep.SettlingTime_s)/ ...
                stepTolerance.SettlingTime_s)^2 + ...
            ((stepPeak-baselineStep.PeakTime_s)/ ...
                stepTolerance.PeakTime_s)^2;

        % MSE is primary. The penalty is scaled relative to the best MSE
        % after all trials are complete, below.
        rows(trial,:) = {lambda,mu,sineMSE,sineRMSE,sineTorque, ...
            stepOvershoot,stepSettling,stepPeak,stepPenalty,NaN};
    end
end

results = cell2table(rows,'VariableNames', ...
    {'Lambda','Mu','SineMSE_rad2','SineRMSE_rad', ...
     'SineMeanAbsTorque_Nm','StepOvershoot_percent', ...
     'StepSettlingTime_s','StepPeakTime_s','StepPenalty','OverallScore'});

% Normalize MSE so the score remains interpretable while preserving MSE as
% the dominant ranking quantity. The secondary step term cannot dominate.
mseScale = min(results.SineMSE_rad2);
results.OverallScore = results.SineMSE_rad2/mseScale + ...
    stepPenaltyWeight*results.StepPenalty;
results = sortrows(results,{'OverallScore','SineMSE_rad2'});
results.Rank = (1:height(results)).';
results = movevars(results,'Rank','Before','Lambda');

fprintf('\nGRID RESULTS (ranked by overall score)\n');
disp(results);

best = results(1,:);
fprintf('\nBEST LAMBDA/MU COMBINATION\n');
fprintf('lambda = %.1f, mu = %.1f\n',best.Lambda,best.Mu);
fprintf('Sine MSE = %.9g rad^2\n',best.SineMSE_rad2);
fprintf('Sine RMSE = %.9g rad\n',best.SineRMSE_rad);
fprintf('Sine accumulated absolute torque = %.9g Nm\n', ...
    best.SineMeanAbsTorque_Nm);
fprintf('Step overshoot = %.6g %%\n',best.StepOvershoot_percent);
fprintf('Step settling time = %.6g s\n',best.StepSettlingTime_s);
fprintf('Step peak time = %.6g s\n',best.StepPeakTime_s);
fprintf('Overall score = %.9g\n',best.OverallScore);

% Save only diagnostic results, never controller baseline parameters.
outFile = fullfile(pwd,'fopid_order_grid_results.csv');
writetable(results,outFile);
fprintf('\nDiagnostic table saved to %s\n',outFile);
