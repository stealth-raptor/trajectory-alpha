function result = simulate_ur5_controller(plant,t,qRef,controllerType,gains,cfg)
%SIMULATE_UR5_CONTROLLER
% Six-joint UR5 simulation for PID/FOPID.
%
% FOPID uses the Oustaloup fractional approximation.

n = numel(t);
h = t(2)-t(1);

%% ============================================================
% Allocate states
% ============================================================

q = zeros(6,n);
dq = zeros(6,n);
tauHist = zeros(6,n);

q(:,1) = cfg.q0(:);
dq(:,1) = cfg.dq0(:);

%% ============================================================
% Controller states
% ============================================================

intState = zeros(6,1);
ePrev = zeros(6,1);

%% ============================================================
% Initialize Oustaloup filters
% ============================================================

if strcmpi(controllerType,'FOPID')

    opts.order = cfg.oustaloup.order;
    opts.wLow = cfg.oustaloup.wLow;
    opts.wHigh = cfg.oustaloup.wHigh;

    fracBank = init_oustaloup_bank( ...
        gains.lambda(:), ...
        gains.mu(:), ...
        h, ...
        opts);

end

%% ============================================================
% Main simulation loop
% ============================================================

for k = 1:n-1

    %% --------------------------------------------------------
    % Current state
    % --------------------------------------------------------

    qk = q(:,k);
    dqk = dq(:,k);

    r = qRef(k,:).';
    r = r(:);

    %% --------------------------------------------------------
    % Tracking error
    % --------------------------------------------------------

    e = r-qk;
    e = e(:);

    %% ========================================================
    % CONTROLLER
    % ========================================================

    if strcmpi(controllerType,'PID')

        %% Standard PID

        intState = intState+h*e;

        intState = min( ...
            max(intState,-cfg.integratorLimit(:)), ...
            cfg.integratorLimit(:));

        de = (e-ePrev)/h;

        tau = ...
            gains.Kp(:).*e + ...
            gains.Ki(:).*intState + ...
            gains.Kd(:).*de;

    elseif strcmpi(controllerType,'FOPID')

        %% ====================================================
        % Oustaloup fractional operators
        % ====================================================

        [If,Df,fracBank] = ...
            oustaloup_step_bank(fracBank,e);

        If = If(:);
        Df = Df(:);

        %% ====================================================
        % FOPID control law
        % ====================================================

        tau = ...
            gains.Kp(:).*e + ...
            gains.Ki(:).*If + ...
            gains.Kd(:).*Df;

    else

        error('Unknown controller type: %s',controllerType);

    end

    %% --------------------------------------------------------
    % Force torque to column vector
    % --------------------------------------------------------

    tau = tau(:);

    % Gravity compensation is applied as model-based feedforward torque.
    if isfield(cfg,'useGravityCompensation') && cfg.useGravityCompensation
        if ~isfield(plant,'gravity')
            error('simulate_ur5_controller:MissingGravityMap', ...
                'Gravity compensation requires plant.gravity(q).');
        end
        tau = tau + plant.gravity(qk);
    end

    %% --------------------------------------------------------
    % Torque saturation
    % --------------------------------------------------------

    if cfg.useTorqueSaturation

        tau = min( ...
            max(tau,-cfg.tauLimit(:)), ...
            cfg.tauLimit(:));

    end

    tauHist(:,k) = tau;

    %% ========================================================
    % State vector
    % ========================================================

    x = [qk;dqk];
    x = x(:);

    %% ========================================================
    % Dynamics
    % ========================================================

    f = @(xx) stateDerivative(xx,tau,plant);

    %% ========================================================
    % Numerical integration
    % ========================================================

    switch lower(cfg.integrationMethod)

        case 'euler'

            dx = f(x);
            dx = dx(:);

            xnext = x+h*dx;

        case 'midpoint'

            k1 = f(x);
            k1 = k1(:);

            k2 = f(x+0.5*h*k1);
            k2 = k2(:);

            xnext = x+h*k2;

        case 'rk4'

            k1 = f(x);
            k1 = k1(:);

            k2 = f(x+0.5*h*k1);
            k2 = k2(:);

            k3 = f(x+0.5*h*k2);
            k3 = k3(:);

            k4 = f(x+h*k3);
            k4 = k4(:);

            xnext = x + ...
                (h/6)*(k1+2*k2+2*k3+k4);

        otherwise

            error( ...
                'Unknown integration method: %s', ...
                cfg.integrationMethod);

    end

    xnext = xnext(:);

    %% --------------------------------------------------------
    % Safety checks
    % --------------------------------------------------------

    if numel(xnext) ~= 12

        error( ...
            'State vector has %d elements instead of 12.', ...
            numel(xnext));

    end

    if any(~isfinite(xnext))

        error( ...
            'Simulation produced NaN/Inf at t = %.4f s.', ...
            t(k));

    end

    %% --------------------------------------------------------
    % Store state
    % --------------------------------------------------------

    q(:,k+1) = xnext(1:6);
    dq(:,k+1) = xnext(7:12);

    ePrev = e;

end

%% ============================================================
% Final torque
% ============================================================

tauHist(:,n) = tauHist(:,n-1);

%% ============================================================
% Output
% ============================================================

result.t = t;
result.q = q.';
result.dq = dq.';
result.tau = tauHist.';

end


%% ============================================================
% STATE DERIVATIVE
% ============================================================

function dx = stateDerivative(x,tau,plant)

x = x(:);
tau = tau(:);

if numel(x) ~= 12

    error('State vector must contain 12 elements.');

end

q = x(1:6);
dq = x(7:12);

ddq = plant.accel(q,dq,tau);
ddq = ddq(:);

if numel(ddq) ~= 6

    error( ...
        'Plant acceleration must contain 6 elements.');

end

dx = [dq;ddq];
dx = dx(:);

end