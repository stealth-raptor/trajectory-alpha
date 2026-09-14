%% RUN_PAPER_OUTPUTS
% Generate the PID-versus-FOPID results used for the paper-style comparison.
%
% The paper reports six-joint step and sine comparisons. Some numerical
% experiment details are not published, so the assumptions are kept visible.

clear; clc; close all;

%% Experiment configuration
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

cfg.optimizeGains = true;
cfg.forceRetune = false;
cfg.tunedGainsFile = fullfile(pwd,'tuned_pid_fopid_gains.mat');
cfg.tuneTs = 0.02;
cfg.tuneTfinal = 2.0;
cfg.tuneMaxIter = 5;
cfg.tuneMaxFunEvals = 12;

%% Baseline controller parameters
% The paper does not publish the final gains. These are explicit baseline
% values used for a reproducible PID-versus-FOPID comparison.
pid.Kp = [120 100 85 55 35 20].';
pid.Ki = [8 7 6 4 3 2].';
pid.Kd = [18 15 12 7 4 2].';

fopid = pid;
fopid.lambda = 0.85*ones(6,1);
fopid.mu = 0.85*ones(6,1);

%% Plant and references
plant = make_ur5_plant(cfg);
t = (0:cfg.Ts:cfg.Tfinal).';
qStep = zeros(numel(t),6);
qStep(t >= cfg.stepTime,:) = cfg.stepAmplitude;
qSine = cfg.sineAmplitude * ...
    sin(2*pi*cfg.sineFrequencyHz*t) * ones(1,6);

fprintf('\nPaper-style PID versus FOPID output generation\n');
fprintf('Plant: %s\n',plant.description);
fprintf('Step: %.3f rad at %.2f s\n', ...
    cfg.stepAmplitude,cfg.stepTime);
fprintf('Sine: %.3f rad amplitude at %.3f Hz\n', ...
    cfg.sineAmplitude,cfg.sineFrequencyHz);

%% Load tuned gains or optimize them once
if cfg.optimizeGains && exist(cfg.tunedGainsFile,'file') == 2 && ...
        ~cfg.forceRetune
    tuned = load(cfg.tunedGainsFile,'pid','fopid');
    pid = tuned.pid;
    fopid = tuned.fopid;
    fprintf('\nLoaded previously tuned gains from %s\n', ...
        cfg.tunedGainsFile);
elseif cfg.optimizeGains
    fprintf('\nTuning gains once using a short bounded search...\n');
    [pid,fopid] = tunePaperGains(pid,fopid,cfg,plant);
    save(cfg.tunedGainsFile,'pid','fopid');
    fprintf('Gain tuning complete and saved to %s\n', ...
        cfg.tunedGainsFile);
else
    fprintf('\nUsing baseline gains; gain tuning is disabled.\n');
end

%% Run the two controllers
fprintf('\nRunning PID step response...\n');
pidStep = simulate_ur5_controller(plant,t,qStep,'PID',pid,cfg);
fprintf('Running FOPID step response...\n');
fopidStep = simulate_ur5_controller(plant,t,qStep,'FOPID',fopid,cfg);

fprintf('Running PID sine response...\n');
pidSine = simulate_ur5_controller(plant,t,qSine,'PID',pid,cfg);
fprintf('Running FOPID sine response...\n');
fopidSine = simulate_ur5_controller(plant,t,qSine,'FOPID',fopid,cfg);

%% Metrics
[stepMetrics,stepSummary] = make_step_metrics( ...
    t,qStep,pidStep.q,fopidStep.q,pidStep.tau,fopidStep.tau);
[sineMetrics,sineSummary] = make_sine_metrics( ...
    t,qSine,pidSine,fopidSine);

fprintf('\nSTEP SUMMARY\n');
disp(stepSummary);
fprintf('\nSINE SUMMARY\n');
disp(sineSummary);

%% Output folders
outRoot = fullfile(pwd,'outputs','paper');
stepFigRoot = fullfile(outRoot,'step','figures');
stepMetricRoot = fullfile(outRoot,'step','metrics');
sineFigRoot = fullfile(outRoot,'sine','figures');
sineMetricRoot = fullfile(outRoot,'sine','metrics');
resultRoot = fullfile(outRoot,'results');

folders = {stepFigRoot,stepMetricRoot,sineFigRoot,sineMetricRoot,resultRoot};
for k = 1:numel(folders)
    if ~exist(folders{k},'dir')
        mkdir(folders{k});
    end
end

colors.pid = [0.10 0.55 0.20];
colors.fopid = [0.00 0.30 0.75];

%% Step position and error figures
for j = 1:6
    fig = figure('Visible','off','Color','w');
    plot(t,qStep(:,j),'k--','LineWidth',1.2); hold on;
    plot(t,pidStep.q(:,j),'Color',colors.pid,'LineWidth',1.1);
    plot(t,fopidStep.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    grid on; box on;
    xlabel('Time (s)'); ylabel('Joint angle (rad)');
    title(sprintf('Joint %d step position tracking',j));
    legend('Target','PID','FOPID','Location','best');
    saveFigure(fig,stepFigRoot,sprintf('joint_%02d_position_tracking',j));

    fig = figure('Visible','off','Color','w');
    plot(t,qStep(:,j)-pidStep.q(:,j),'Color',colors.pid,'LineWidth',1.1); hold on;
    plot(t,qStep(:,j)-fopidStep.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    yline(0,'k:'); grid on; box on;
    xlabel('Time (s)'); ylabel('Tracking error (rad)');
    title(sprintf('Joint %d step tracking error',j));
    legend('PID error','FOPID error','Location','best');
    saveFigure(fig,stepFigRoot,sprintf('joint_%02d_tracking_error',j));
end

fig = figure('Visible','off','Color','w');
tiledlayout(3,2,'TileSpacing','compact');
for j = 1:6
    nexttile;
    plot(t,pidStep.tau(:,j),'Color',colors.pid,'LineWidth',1.0); hold on;
    plot(t,fopidStep.tau(:,j),'Color',colors.fopid,'LineWidth',1.0);
    grid on; box on; title(sprintf('Joint %d torque',j));
    xlabel('Time (s)'); ylabel('\tau (N m)');
    if j == 1
        legend('PID','FOPID','Location','best');
    end
end
saveFigure(fig,stepFigRoot,'all_joints_torque_comparison');

%% Sine position and error figures
for j = 1:6
    fig = figure('Visible','off','Color','w');
    plot(t,qSine(:,j),'k--','LineWidth',1.2); hold on;
    plot(t,pidSine.q(:,j),'Color',colors.pid,'LineWidth',1.1);
    plot(t,fopidSine.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    grid on; box on;
    xlabel('Time (s)'); ylabel('Joint angle (rad)');
    title(sprintf('Joint %d sine position tracking',j));
    legend('Target','PID','FOPID','Location','best');
    saveFigure(fig,sineFigRoot,sprintf('joint_%02d_position_tracking',j));

    fig = figure('Visible','off','Color','w');
    plot(t,qSine(:,j)-pidSine.q(:,j),'Color',colors.pid,'LineWidth',1.1); hold on;
    plot(t,qSine(:,j)-fopidSine.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    yline(0,'k:'); grid on; box on;
    xlabel('Time (s)'); ylabel('Tracking error (rad)');
    title(sprintf('Joint %d sine tracking error',j));
    legend('PID error','FOPID error','Location','best');
    saveFigure(fig,sineFigRoot,sprintf('joint_%02d_tracking_error',j));
end

fig = figure('Visible','off','Color','w');
tiledlayout(3,2,'TileSpacing','compact');
for j = 1:6
    nexttile;
    plot(t,pidSine.tau(:,j),'Color',colors.pid,'LineWidth',1.0); hold on;
    plot(t,fopidSine.tau(:,j),'Color',colors.fopid,'LineWidth',1.0);
    grid on; box on; title(sprintf('Joint %d torque',j));
    xlabel('Time (s)'); ylabel('\tau (N m)');
    if j == 1
        legend('PID','FOPID','Location','best');
    end
end
saveFigure(fig,sineFigRoot,'all_joints_torque_comparison');

%% Save paper-style tables and complete results
writetable(stepMetrics,fullfile(stepMetricRoot,'step_joint_metrics.csv'));
writetable(stepSummary,fullfile(stepMetricRoot,'step_summary_metrics.csv'));
writetable(sineMetrics,fullfile(sineMetricRoot,'sine_joint_metrics.csv'));
writetable(sineSummary,fullfile(sineMetricRoot,'sine_summary_metrics.csv'));

save(fullfile(resultRoot,'paper_pid_fopid_results.mat'), ...
    'cfg','plant','pid','fopid','t','qStep','qSine', ...
    'pidStep','fopidStep','pidSine','fopidSine', ...
    'stepMetrics','stepSummary','sineMetrics','sineSummary');

fprintf('\nSaved paper-style PID/FOPID results under outputs/paper/\n');

%% Local helper
function [pid,fopid] = tunePaperGains(pid,fopid,cfg,plant)
cfgTune = cfg;
cfgTune.Ts = cfg.tuneTs;
cfgTune.Tfinal = cfg.tuneTfinal;
cfgTune.integrationMethod = 'euler';
tTune = (0:cfgTune.Ts:cfgTune.Tfinal).';
qStepTune = zeros(numel(tTune),6);
qStepTune(tTune >= cfgTune.stepTime,:) = cfgTune.stepAmplitude;
qSineTune = cfgTune.sineAmplitude * ...
    sin(2*pi*cfgTune.sineFrequencyHz*tTune) * ones(1,6);

% The paper gives target averages but not final gains. Search positive
% shared gain scales so the six-joint baseline remains interpretable.
pidBase = pid;
fopidBase = fopid;
opts = optimset('Display','off','MaxIter',cfg.tuneMaxIter, ...
    'MaxFunEvals',cfg.tuneMaxFunEvals,'TolX',1e-2,'TolFun',1e-3);

xPid = fminsearch(@(x) paperObjective(x,'PID',pidBase,fopidBase, ...
    plant,tTune,qStepTune,qSineTune,cfgTune),zeros(3,1),opts);
pid.Kp = pidBase.Kp*exp(xPid(1));
pid.Ki = pidBase.Ki*exp(xPid(2));
pid.Kd = pidBase.Kd*exp(xPid(3));

xFopid = fminsearch(@(x) paperObjective(x,'FOPID',pidBase,fopidBase, ...
    plant,tTune,qStepTune,qSineTune,cfgTune),zeros(5,1),opts);
fopid.Kp = fopidBase.Kp*exp(xFopid(1));
fopid.Ki = fopidBase.Ki*exp(xFopid(2));
fopid.Kd = fopidBase.Kd*exp(xFopid(3));
fopid.lambda = min(max(fopidBase.lambda*exp(xFopid(4)),0.20),1.50);
fopid.mu = min(max(fopidBase.mu*exp(xFopid(5)),0.20),1.50);
end

function score = paperObjective(x,kind,pidBase,fopidBase,plant,t,qStep,qSine,cfg)
% Bound logarithmic gain/order changes before evaluating the robot model.
% This prevents fminsearch from probing physically meaningless unstable
% controllers that can drive the rigid-body state to Inf or NaN.
if any(abs(x) > 0.5)
    score = 1e9 + 1e6*sum(max(abs(x)-0.5,0).^2);
    return;
end

if strcmpi(kind,'PID')
    gains = pidBase;
    gains.Kp = gains.Kp*exp(x(1));
    gains.Ki = gains.Ki*exp(x(2));
    gains.Kd = gains.Kd*exp(x(3));
else
    gains = fopidBase;
    gains.Kp = gains.Kp*exp(x(1));
    gains.Ki = gains.Ki*exp(x(2));
    gains.Kd = gains.Kd*exp(x(3));
    gains.lambda = min(max(gains.lambda*exp(x(4)),0.20),1.50);
    gains.mu = min(max(gains.mu*exp(x(5)),0.20),1.50);
end

try
    warningState = warning;
    warning('off','all');
    warningCleanup = onCleanup(@()warning(warningState)); %#ok<NASGU>
    stepResult = simulate_ur5_controller(plant,t,qStep,kind,gains,cfg);
    sineResult = simulate_ur5_controller(plant,t,qSine,kind,gains,cfg);
    if strcmpi(kind,'PID')
        stepMetrics = make_step_metrics(t,qStep,stepResult.q,[]);
    else
        stepMetrics = make_step_metrics(t,qStep,[],stepResult.q);
    end
    overshoot = mean(stepMetrics.Overshoot_percent,'omitnan');
    settling = mean(stepMetrics.SettlingTime_s,'omitnan');
    peak = mean(stepMetrics.PeakTime_s,'omitnan');
    sineMSE = mean((qSine-sineResult.q).^2,'all');
    if ~isfinite(overshoot), overshoot = 1e3; end
    if ~isfinite(settling), settling = 1e3; end
    if ~isfinite(peak), peak = 1e3; end
    score = ((overshoot-targetOvershoot(kind))/25)^2 + ...
        ((settling-targetSettling(kind))/1.5)^2 + ...
        ((peak-targetPeak(kind))/0.7)^2 + ...
        ((sineMSE-targetMSE(kind))/0.03)^2;
    if ~isfinite(score), score = 1e12; end
catch
    score = 1e12;
end
end

function value = targetOvershoot(kind)
if strcmpi(kind,'PID'), value = 54.6; else, value = 31.2; end
end

function value = targetSettling(kind)
if strcmpi(kind,'PID'), value = 2.44; else, value = 1.89; end
end

function value = targetPeak(kind)
if strcmpi(kind,'PID'), value = 1.39; else, value = 1.33; end
end

function value = targetMSE(kind)
if strcmpi(kind,'PID'), value = 2.27e-2; else, value = 0.88e-2; end
end

function saveFigure(fig,folder,baseName)
print(fig,fullfile(folder,[baseName,'.png']),'-dpng','-r300');
savefig(fig,fullfile(folder,[baseName,'.fig']));
close(fig);
end
