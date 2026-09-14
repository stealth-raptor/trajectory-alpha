function plant = make_ur5_plant(cfg)
%MAKE_UR5_PLANT Prefer a coupled rigidBodyTree UR5 plant.
% If Robotics System Toolbox/loadrobot is unavailable, use a clearly marked
% decoupled approximation based on the masses, lengths and diagonal inertia
% values listed in Table 2 of the paper.

% Older experiment scripts do not define friction explicitly. Use the
% project's documented defaults while allowing callers to override them.
if ~isfield(cfg,'friction')
    cfg.friction.viscous = 0.05*ones(6,1);
    cfg.friction.coulomb = 0.20*ones(6,1);
    cfg.friction.velocityScale = 1e-3;
else
    if ~isfield(cfg.friction,'viscous')
        cfg.friction.viscous = 0.05*ones(6,1);
    end
    if ~isfield(cfg.friction,'coulomb')
        cfg.friction.coulomb = 0.20*ones(6,1);
    end
    if ~isfield(cfg.friction,'velocityScale')
        cfg.friction.velocityScale = 1e-3;
    end
end

plant.friction = cfg.friction;
hasRST = exist('loadrobot','file') == 2;

if hasRST
    try
        robot = loadrobot('universalUR5','DataFormat','column', ...
            'Gravity',[0 0 -9.81]);

        % Reject malformed or singular robot assets before RK4 can spread
        % invalid dynamics values through the whole state.
        qTest = zeros(6,1);
        dqTest = zeros(6,1);
        MTest = massMatrix(robot,qTest);
        hTest = velocityProduct(robot,qTest,dqTest);
        GTest = gravityTorque(robot,qTest);
        valid = all(isfinite(MTest),'all') && all(isfinite(hTest)) && ...
            all(isfinite(GTest)) && rcond(MTest) > 1e-12;

        if valid
            plant.type = 'rigidBodyTree';
            plant.robot = robot;
            plant.description = 'coupled rigidBodyTree universalUR5';
            plant.accel = @(q,dq,tau) ...
                rigidBodyTreeAccel(robot,q,dq,tau,cfg);
            return;
        end

        warning(['universalUR5 preflight returned invalid or singular ', ...
            'dynamics. Falling back to the Table-2 approximation.']);
    catch ME
        warning('UR5Plant:LoadRobotFailed', ...
            'loadrobot failed (%s). Falling back to decoupled model.', ...
            ME.message);
    end
end

% Table 2 values from the paper. The paper does not provide a complete
% executable UR5 parameter set, so this is only a fallback approximation.
m = [2.0 2.5 5.7 3.9 2.5 2.5].';
L = [0.128 0.612 0.571 0.164 0.115 0.092].';
Izz = [1.0e-3 4.0e-2 4.0e-2 4.0e-2 4.0e-2 4.0e-2].';
J = Izz + m.*(L/2).^2;
g = 9.81;

plant.type = 'decoupledFallback';
plant.description = ...
    'decoupled Table-2 approximation (not coupled UR5 dynamics)';
plant.J = J;
plant.m = m;
plant.L = L;
plant.g = g;
plant.accel = @(q,dq,tau) ...
    decoupledAccel(q,dq,tau,J,m,L,g,cfg);
end

function ddq = rigidBodyTreeAccel(robot,q,dq,tau,cfg)
validateattributes(q,{'double','single'}, ...
    {'real','finite','vector','numel',6});
validateattributes(dq,{'double','single'}, ...
    {'real','finite','vector','numel',6});
validateattributes(tau,{'double','single'}, ...
    {'real','finite','vector','numel',6});

M = massMatrix(robot,q);
h = velocityProduct(robot,q,dq);
G = gravityTorque(robot,q);
fr = cfg.friction.viscous.*dq + cfg.friction.coulomb .* ...
    tanh(dq/cfg.friction.velocityScale);

ddq = M \ (tau - h - G - fr);

if any(~isfinite(ddq))
    error(['UR5 rigidBodyTree dynamics produced NaN/Inf acceleration ', ...
        'at q = [%s].'],sprintf(' %.5g',q));
end
end

function ddq = decoupledAccel(q,dq,tau,J,m,L,g,cfg)
q = q(:);
dq = dq(:);
tau = tau(:);

fr = cfg.friction.viscous.*dq + cfg.friction.coulomb .* ...
    tanh(dq/cfg.friction.velocityScale);
gravity = m.*g.*(L/2).*sin(q);
ddq = (tau - fr - gravity)./J;
ddq = ddq(:);
end