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

fprintf('\nUsing hardcoded PID/FOPID gains; no optimizer is used.\n');

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
function saveFigure(fig,folder,baseName)
print(fig,fullfile(folder,[baseName,'.png']),'-dpng','-r300');
savefig(fig,fullfile(folder,[baseName,'.fig']));
close(fig);
end
