%% RUN_PAPER_OUTPUTS
% Generate the PID-versus-FOPID results used for the paper-style comparison,
% then tune the same controllers with a standard Grey Wolf Optimizer (GWO)
% using the paper ITAE fitness. Plant, controller, trajectories, simulation
% settings and metric definitions are unchanged.

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
cfg.useGravityCompensation = true;
cfg.integrationMethod = 'rk4';

cfg.friction.viscous = 0.05*ones(6,1);
cfg.friction.coulomb = 0.20*ones(6,1);
cfg.friction.velocityScale = 1e-3;

cfg.oustaloup.order = 5;
cfg.oustaloup.wLow = 1e-3;
cfg.oustaloup.wHigh = 1e3;
cfg.oustaloup.method = 'tustin';


%% Baseline controller parameters
% Hardcoded hand-tuned baseline values. The paper does not publish its
% final gains, so these values are selected by direct simulation trials.
pid.Kp = [174.1 145.1 123.3 79.8 50.8 29.0].';
pid.Ki = [11.6 10.2 8.70 5.80 4.35 2.90].';
pid.Kd = [9.60 8.00 6.40 3.73 2.13 1.07].';

fopid = pid;
fopid.Kp = [137.2 114.3 97.2 62.9 40.0 22.9].';
fopid.Ki = [9.14 8.00 6.86 4.57 3.43 2.29].';
fopid.Kd = [22.1 18.4 14.8 8.61 4.92 2.46].';
fopid.lambda = 0.85*ones(6,1);
fopid.mu = 0.85*ones(6,1);

%% GWO settings
% Standard GWO (Mirjalili 2014). Fitness is closed-loop ITAE on the same
% step trajectory and plant. Search uses a coarser Ts (tuneTs) only to make
% the swarm evaluations tractable; all printed metrics are full cfg.Ts.
gwo.enabled = true;
gwo.fastMode = false;  % Full deep optimization search
if gwo.fastMode
    gwo.nWolves = 4;
    gwo.maxIter = 6;
    gwo.tuneTs = 0.02;
else
    gwo.nWolves = 10;
    gwo.maxIter = 15;
    gwo.tuneTs = 0.005;  % Finer grid matching simulation timestep
end
gwo.rngSeed = 42;
gwo.useCache = false;  % Disable cache to ensure fresh GWO output and graphs are generated

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

fprintf('\nBaseline PID/FOPID use hardcoded gains; GWO then retunes copies of those gains.\n');

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

%% GWO-PID and GWO-FOPID (same plant, trajectories, cfg, metrics)
gwoPid = pid;
gwoFopid = fopid;
gwoPidTune = [];
gwoFopidTune = [];

if gwo.enabled
    gwoCacheFile = fullfile(pwd,'outputs','paper','results','gwo_pid_fopid_cache.mat');
    loadedCache = false;
    if gwo.useCache && exist(gwoCacheFile,'file')
        cache = load(gwoCacheFile);
        if isfield(cache,'gwo') && isequal(cache.gwo.nWolves,gwo.nWolves) ...
                && isequal(cache.gwo.maxIter,gwo.maxIter) ...
                && isequal(cache.gwo.rngSeed,gwo.rngSeed) ...
                && isequal(cache.gwo.tuneTs,gwo.tuneTs)
            gwoPid = cache.gwoPid;
            gwoFopid = cache.gwoFopid;
            gwoPidTune = cache.gwoPidTune;
            gwoFopidTune = cache.gwoFopidTune;
            loadedCache = true;
            fprintf('\nLoaded cached GWO gains from %s\n',gwoCacheFile);
        end
    end

    if ~loadedCache
        fprintf('\nRunning GWO-PID on step ITAE (%d wolves, %d iterations)...\n', ...
            gwo.nWolves,gwo.maxIter);
        [gwoPid,gwoPidTune] = gwo_optimize_controller( ...
            plant,t,qStep,'PID',pid,cfg,gwo);

        fprintf('\nRunning GWO-FOPID on step ITAE (%d wolves, %d iterations)...\n', ...
            gwo.nWolves,gwo.maxIter);
        [gwoFopid,gwoFopidTune] = gwo_optimize_controller( ...
            plant,t,qStep,'FOPID',fopid,cfg,gwo);
    end
else
    fprintf('\nGWO disabled; GWO-PID/GWO-FOPID reuse the hardcoded baselines.\n');
end

fprintf('\nRunning GWO-PID step response...\n');
gwoPidStep = simulate_ur5_controller(plant,t,qStep,'PID',gwoPid,cfg);
fprintf('Running GWO-FOPID step response...\n');
gwoFopidStep = simulate_ur5_controller(plant,t,qStep,'FOPID',gwoFopid,cfg);
fprintf('Running GWO-PID sine response...\n');
gwoPidSine = simulate_ur5_controller(plant,t,qSine,'PID',gwoPid,cfg);
fprintf('Running GWO-FOPID sine response...\n');
gwoFopidSine = simulate_ur5_controller(plant,t,qSine,'FOPID',gwoFopid,cfg);

[stepMetricsGwo,stepSummaryGwo] = make_step_metrics( ...
    t,qStep,gwoPidStep.q,gwoFopidStep.q,gwoPidStep.tau,gwoFopidStep.tau);
[sineMetricsGwo,sineSummaryGwo] = make_sine_metrics( ...
    t,qSine,gwoPidSine,gwoFopidSine);

stepSummaryGwo.Controller = ["GWO-PID";"GWO-FOPID"];
sineSummaryGwo.Controller = ["GWO-PID";"GWO-FOPID"];
stepMetricsGwo.Controller = renameGwoControllers(stepMetricsGwo.Controller);
sineMetricsGwo.Controller = renameGwoControllers(sineMetricsGwo.Controller);

printFourWayComparison( ...
    t,qStep,qSine, ...
    pidStep,fopidStep,gwoPidStep,gwoFopidStep, ...
    pidSine,fopidSine,gwoPidSine,gwoFopidSine, ...
    stepSummary,stepSummaryGwo,sineSummary,sineSummaryGwo);

fprintf('\nGWO-PID best gains (Kp, Ki, Kd):\n');
disp([gwoPid.Kp, gwoPid.Ki, gwoPid.Kd]);
fprintf('GWO-FOPID best gains (Kp, Ki, Kd, lambda, mu):\n');
disp([gwoFopid.Kp, gwoFopid.Ki, gwoFopid.Kd, gwoFopid.lambda, gwoFopid.mu]);
if ~isempty(gwoPidTune) && isfield(gwoPidTune,'bestFitness')
    fprintf('GWO-PID best step ITAE = %.6g\n',gwoPidTune.bestFitness);
end
if ~isempty(gwoFopidTune) && isfield(gwoFopidTune,'bestFitness')
    fprintf('GWO-FOPID best step ITAE = %.6g\n',gwoFopidTune.bestFitness);
end

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
colors.gwoPid = [0.85 0.35 0.00];
colors.gwoFopid = [0.60 0.00 0.60];

%% Step position and error figures (4-way comparison)
for j = 1:6
    fig = figure('Visible','off','Color','w');
    plot(t,qStep(:,j),'k--','LineWidth',1.2); hold on;
    plot(t,pidStep.q(:,j),'Color',colors.pid,'LineWidth',1.1);
    plot(t,fopidStep.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    plot(t,gwoPidStep.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.');
    plot(t,gwoFopidStep.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1,'LineStyle','-');
    grid on; box on;
    xlabel('Time (s)'); ylabel('Joint angle (rad)');
    title(sprintf('Joint %d step position tracking',j));
    legend('Target','PID','FOPID','GWO-PID','GWO-FOPID','Location','best');
    saveFigure(fig,stepFigRoot,sprintf('joint_%02d_position_tracking',j));

    fig = figure('Visible','off','Color','w');
    plot(t,qStep(:,j)-pidStep.q(:,j),'Color',colors.pid,'LineWidth',1.1); hold on;
    plot(t,qStep(:,j)-fopidStep.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    plot(t,qStep(:,j)-gwoPidStep.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.');
    plot(t,qStep(:,j)-gwoFopidStep.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1,'LineStyle','-');
    yline(0,'k:'); grid on; box on;
    xlabel('Time (s)'); ylabel('Tracking error (rad)');
    title(sprintf('Joint %d step tracking error',j));
    legend('PID error','FOPID error','GWO-PID error','GWO-FOPID error','Location','best');
    saveFigure(fig,stepFigRoot,sprintf('joint_%02d_tracking_error',j));
end

fig = figure('Visible','off','Color','w');
tiledlayout(3,2,'TileSpacing','compact');
for j = 1:6
    nexttile;
    plot(t,pidStep.tau(:,j),'Color',colors.pid,'LineWidth',1.0); hold on;
    plot(t,fopidStep.tau(:,j),'Color',colors.fopid,'LineWidth',1.0);
    plot(t,gwoPidStep.tau(:,j),'Color',colors.gwoPid,'LineWidth',1.0,'LineStyle','-.');
    plot(t,gwoFopidStep.tau(:,j),'Color',colors.gwoFopid,'LineWidth',1.0);
    grid on; box on; title(sprintf('Joint %d torque',j));
    xlabel('Time (s)'); ylabel('\tau (N m)');
    if j == 1
        legend('PID','FOPID','GWO-PID','GWO-FOPID','Location','best');
    end
end
saveFigure(fig,stepFigRoot,'all_joints_torque_comparison');

%% Sine position and error figures (4-way comparison)
for j = 1:6
    fig = figure('Visible','off','Color','w');
    plot(t,qSine(:,j),'k--','LineWidth',1.2); hold on;
    plot(t,pidSine.q(:,j),'Color',colors.pid,'LineWidth',1.1);
    plot(t,fopidSine.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    plot(t,gwoPidSine.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.');
    plot(t,gwoFopidSine.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1,'LineStyle','-');
    grid on; box on;
    xlabel('Time (s)'); ylabel('Joint angle (rad)');
    title(sprintf('Joint %d sine position tracking',j));
    legend('Target','PID','FOPID','GWO-PID','GWO-FOPID','Location','best');
    saveFigure(fig,sineFigRoot,sprintf('joint_%02d_position_tracking',j));

    fig = figure('Visible','off','Color','w');
    plot(t,qSine(:,j)-pidSine.q(:,j),'Color',colors.pid,'LineWidth',1.1); hold on;
    plot(t,qSine(:,j)-fopidSine.q(:,j),'Color',colors.fopid,'LineWidth',1.1);
    plot(t,qSine(:,j)-gwoPidSine.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.');
    plot(t,qSine(:,j)-gwoFopidSine.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1,'LineStyle','-');
    yline(0,'k:'); grid on; box on;
    xlabel('Time (s)'); ylabel('Tracking error (rad)');
    title(sprintf('Joint %d sine tracking error',j));
    legend('PID error','FOPID error','GWO-PID error','GWO-FOPID error','Location','best');
    saveFigure(fig,sineFigRoot,sprintf('joint_%02d_tracking_error',j));
end

fig = figure('Visible','off','Color','w');
tiledlayout(3,2,'TileSpacing','compact');
for j = 1:6
    nexttile;
    plot(t,pidSine.tau(:,j),'Color',colors.pid,'LineWidth',1.0); hold on;
    plot(t,fopidSine.tau(:,j),'Color',colors.fopid,'LineWidth',1.0);
    plot(t,gwoPidSine.tau(:,j),'Color',colors.gwoPid,'LineWidth',1.0,'LineStyle','-.');
    plot(t,gwoFopidSine.tau(:,j),'Color',colors.gwoFopid,'LineWidth',1.0);
    grid on; box on; title(sprintf('Joint %d torque',j));
    xlabel('Time (s)'); ylabel('\tau (N m)');
    if j == 1
        legend('PID','FOPID','GWO-PID','GWO-FOPID','Location','best');
    end
end
saveFigure(fig,sineFigRoot,'all_joints_torque_comparison');

%% Save paper-style tables and complete results
writetable(stepMetrics,fullfile(stepMetricRoot,'step_joint_metrics.csv'));
writetable(stepSummary,fullfile(stepMetricRoot,'step_summary_metrics.csv'));
writetable(sineMetrics,fullfile(sineMetricRoot,'sine_joint_metrics.csv'));
writetable(sineSummary,fullfile(sineMetricRoot,'sine_summary_metrics.csv'));

save(fullfile(resultRoot,'paper_pid_fopid_results.mat'), ...
    'cfg','plant','pid','fopid','gwo','gwoPid','gwoFopid', ...
    'gwoPidTune','gwoFopidTune','t','qStep','qSine', ...
    'pidStep','fopidStep','gwoPidStep','gwoFopidStep', ...
    'pidSine','fopidSine','gwoPidSine','gwoFopidSine', ...
    'stepMetrics','stepSummary','sineMetrics','sineSummary', ...
    'stepMetricsGwo','stepSummaryGwo','sineMetricsGwo','sineSummaryGwo');

if gwo.enabled && ~isempty(gwoPidTune) && ~isempty(gwoFopidTune)
    save(fullfile(resultRoot,'gwo_pid_fopid_cache.mat'), ...
        'gwo','gwoPid','gwoFopid','gwoPidTune','gwoFopidTune');

    gwoFigRoot = fullfile(outRoot,'gwo','figures');
    if ~exist(gwoFigRoot,'dir')
        mkdir(gwoFigRoot);
    end
    
    % Convergence plot
    fig = figure('Visible','off','Color','w');
    plot(1:numel(gwoPidTune.fitnessHistory),gwoPidTune.fitnessHistory, ...
        'LineWidth',1.2); hold on;
    plot(1:numel(gwoFopidTune.fitnessHistory),gwoFopidTune.fitnessHistory, ...
        'LineWidth',1.2);
    grid on; box on;
    xlabel('Iteration'); ylabel('Best ITAE');
    title('GWO convergence (step ITAE)');
    legend('GWO-PID','GWO-FOPID','Location','best');
    saveFigure(fig,gwoFigRoot,'gwo_itae_convergence');

    % Dedicated GWO Step position tracking (3x2 grid)
    fig = figure('Visible','off','Color','w');
    tiledlayout(3,2,'TileSpacing','compact');
    for j = 1:6
        nexttile;
        plot(t,qStep(:,j),'k--','LineWidth',1.0); hold on;
        plot(t,gwoPidStep.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.');
        plot(t,gwoFopidStep.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1);
        grid on; box on; title(sprintf('Joint %d',j));
        xlabel('Time (s)'); ylabel('q (rad)');
        if j == 1, legend('Target','GWO-PID','GWO-FOPID','Location','best'); end
    end
    saveFigure(fig,gwoFigRoot,'gwo_step_position_tracking');

    % Dedicated GWO Step tracking error (3x2 grid)
    fig = figure('Visible','off','Color','w');
    tiledlayout(3,2,'TileSpacing','compact');
    for j = 1:6
        nexttile;
        plot(t,qStep(:,j)-gwoPidStep.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.'); hold on;
        plot(t,qStep(:,j)-gwoFopidStep.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1);
        yline(0,'k:'); grid on; box on; title(sprintf('Joint %d error',j));
        xlabel('Time (s)'); ylabel('Error (rad)');
        if j == 1, legend('GWO-PID','GWO-FOPID','Location','best'); end
    end
    saveFigure(fig,gwoFigRoot,'gwo_step_tracking_error');

    % Dedicated GWO Sine position tracking (3x2 grid)
    fig = figure('Visible','off','Color','w');
    tiledlayout(3,2,'TileSpacing','compact');
    for j = 1:6
        nexttile;
        plot(t,qSine(:,j),'k--','LineWidth',1.0); hold on;
        plot(t,gwoPidSine.q(:,j),'Color',colors.gwoPid,'LineWidth',1.1,'LineStyle','-.');
        plot(t,gwoFopidSine.q(:,j),'Color',colors.gwoFopid,'LineWidth',1.1);
        grid on; box on; title(sprintf('Joint %d',j));
        xlabel('Time (s)'); ylabel('q (rad)');
        if j == 1, legend('Target','GWO-PID','GWO-FOPID','Location','best'); end
    end
    saveFigure(fig,gwoFigRoot,'gwo_sine_position_tracking');

    % Dedicated GWO Torque comparison (3x2 grid)
    fig = figure('Visible','off','Color','w');
    tiledlayout(3,2,'TileSpacing','compact');
    for j = 1:6
        nexttile;
        plot(t,gwoPidStep.tau(:,j),'Color',colors.gwoPid,'LineWidth',1.0,'LineStyle','-.'); hold on;
        plot(t,gwoFopidStep.tau(:,j),'Color',colors.gwoFopid,'LineWidth',1.0);
        grid on; box on; title(sprintf('Joint %d torque',j));
        xlabel('Time (s)'); ylabel('\tau (N m)');
        if j == 1, legend('GWO-PID','GWO-FOPID','Location','best'); end
    end
    saveFigure(fig,gwoFigRoot,'gwo_torque_comparison');
end

writetable(stepMetricsGwo,fullfile(stepMetricRoot,'step_joint_metrics_gwo.csv'));
writetable(stepSummaryGwo,fullfile(stepMetricRoot,'step_summary_metrics_gwo.csv'));
writetable(sineMetricsGwo,fullfile(sineMetricRoot,'sine_joint_metrics_gwo.csv'));
writetable(sineSummaryGwo,fullfile(sineMetricRoot,'sine_summary_metrics_gwo.csv'));

fprintf('\nSaved paper-style PID/FOPID and GWO results under outputs/paper/\n');

%% Local helper
function saveFigure(fig,folder,baseName)
print(fig,fullfile(folder,[baseName,'.png']),'-dpng','-r300');
savefig(fig,fullfile(folder,[baseName,'.fig']));
close(fig);
end

function names = renameGwoControllers(names)
names = string(names);
names(names=="PID") = "GWO-PID";
names(names=="FOPID") = "GWO-FOPID";
end

function printFourWayComparison( ...
    t,qStep,qSine, ...
    pidStep,fopidStep,gwoPidStep,gwoFopidStep, ...
    pidSine,fopidSine,gwoPidSine,gwoFopidSine, ...
    stepSummary,stepSummaryGwo,sineSummary,sineSummaryGwo)

labels = ["PID","FOPID","GWO-PID","GWO-FOPID"];
qStepAll = {pidStep.q,fopidStep.q,gwoPidStep.q,gwoFopidStep.q};
qSineAll = {pidSine.q,fopidSine.q,gwoPidSine.q,gwoFopidSine.q};
stepTab = [stepSummary; stepSummaryGwo];
sineTab = [sineSummary; sineSummaryGwo];
stepTab.Controller = string(stepTab.Controller);
sineTab.Controller = string(sineTab.Controller);

itaeStep = zeros(4,1);
itaeSine = zeros(4,1);
for k = 1:4
    itaeStep(k) = compute_itae(t,qStep,qStepAll{k});
    itaeSine(k) = compute_itae(t,qSine,qSineAll{k});
end

fprintf('\n============================================================\n');
fprintf('FOUR-WAY COMPARISON  (PID, FOPID, GWO-PID, GWO-FOPID)\n');
fprintf('ITAE is the paper fitness on the same trajectories.\n');
fprintf('Step metrics: overshoot, settling time, peak time, |tau|.\n');
fprintf('Sine metrics: MSE and |tau|.\n');
fprintf('============================================================\n');
fprintf('%-10s %12s %12s %8s %8s %8s %12s %14s\n', ...
    'Ctrl','ITAE_step','ITAE_sine','OS_%','ts_s','tp_s','MSE_rad2','|tau|_sine');

for k = 1:4
    sIdx = find(stepTab.Controller==labels(k),1);
    nIdx = find(sineTab.Controller==labels(k),1);
    fprintf('%-10s %12.6g %12.6g %8.3f %8.3f %8.3f %12.6g %14.6g\n', ...
        labels(k),itaeStep(k),itaeSine(k), ...
        stepTab.AverageOvershoot_percent(sIdx), ...
        stepTab.AverageSettlingTime_s(sIdx), ...
        stepTab.AveragePeakTime_s(sIdx), ...
        sineTab.AverageMSE_rad2(nIdx), ...
        sineTab.AverageMeanAbsTorque_Nm(nIdx));
end

fprintf('\nStep average |tau| (N m):\n');
for k = 1:4
    sIdx = find(stepTab.Controller==labels(k),1);
    fprintf('  %-10s  %.6g\n',labels(k),stepTab.AverageAbsTorque_Nm(sIdx));
end

if isfield(pidStep,'tau')
    % Keep the original per-controller summary tables visible as well.
    fprintf('\nSTEP SUMMARY (PID / FOPID)\n');
    disp(stepSummary);
    fprintf('STEP SUMMARY (GWO-PID / GWO-FOPID)\n');
    disp(stepSummaryGwo);
    fprintf('SINE SUMMARY (PID / FOPID)\n');
    disp(sineSummary);
    fprintf('SINE SUMMARY (GWO-PID / GWO-FOPID)\n');
    disp(sineSummaryGwo);
end
end
