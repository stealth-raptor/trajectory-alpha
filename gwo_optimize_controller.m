function [bestGains,gwoOut] = gwo_optimize_controller( ...
    plant,t,qRef,controllerType,seedGains,cfg,gwoCfg)
%GWO_OPTIMIZE_CONTROLLER Standard Grey Wolf Optimizer for PID/FOPID gains.
%
%   Canonical GWO (Mirjalili, Mirjalili & Lewis, 2014): each search agent
%   (wolf) is a complete controller-parameter vector. Fitness is the
%   paper ITAE from a full closed-loop UR5 simulation. The plant,
%   controller, trajectory and simulation settings are not modified.
%
%   PID decision vector  (18 x 1):  [Kp(1:6); Ki(1:6); Kd(1:6)]
%   FOPID decision vector (30 x 1): [Kp; Ki; Kd; lambda; mu]

controllerType = upper(strtrim(controllerType));
gwoCfg = applyGwoDefaults(gwoCfg);

if isfield(gwoCfg,'rngSeed') && ~isempty(gwoCfg.rngSeed)
    rng(gwoCfg.rngSeed);
end

% Fitness uses the frozen plant/controller. A coarser Ts may be used only
% inside the search; callers still evaluate the best gains at cfg.Ts.
[cfgFit,tFit,qFit] = makeFitnessGrid(plant,t,qRef,cfg,gwoCfg);

[lb,ub] = gwo_parameter_bounds(controllerType,seedGains,gwoCfg);
dim = numel(lb);
nWolves = gwoCfg.nWolves;
maxIter = gwoCfg.maxIter;

positions = zeros(nWolves,dim);
fitness = inf(nWolves,1);

% Seed wolf 1 with the existing (frozen) gains so GWO cannot do worse than
% the baseline unless a strictly better ITAE vector is found.
positions(1,:) = packGains(controllerType,seedGains).';
positions(1,:) = min(max(positions(1,:),lb),ub);
for i = 2:nWolves
    positions(i,:) = lb + rand(1,dim).*(ub-lb);
end

fprintf('GWO-%s: %d wolves, %d iterations, %d parameters\n', ...
    controllerType,nWolves,maxIter,dim);
fprintf('  fitness grid Ts = %.4g s, Tfinal = %.3f s (%d samples)\n', ...
    tFit(2)-tFit(1),tFit(end),numel(tFit));

for i = 1:nWolves
    fitness(i) = evaluateItae(plant,tFit,qFit,controllerType, ...
        positions(i,:),cfgFit);
end

[alphaScore,alphaPos,betaScore,betaPos,deltaScore,deltaPos] = ...
    rankLeaders(positions,fitness);

bestHistory = zeros(maxIter,1);
meanHistory = zeros(maxIter,1);

for iter = 1:maxIter
    a = 2 - iter*(2/maxIter);

    positions = gwoUpdatePositions(positions, ...
        alphaPos,betaPos,deltaPos,a,lb,ub);

    for i = 1:nWolves
        fitness(i) = evaluateItae(plant,tFit,qFit,controllerType, ...
            positions(i,:),cfgFit);
        [alphaScore,alphaPos,betaScore,betaPos,deltaScore,deltaPos] = ...
            updateLeaders(positions(i,:),fitness(i), ...
            alphaScore,alphaPos,betaScore,betaPos,deltaScore,deltaPos);
    end

    bestHistory(iter) = alphaScore;
    meanHistory(iter) = mean(fitness(isfinite(fitness)));
    if ~isfinite(meanHistory(iter))
        meanHistory(iter) = inf;
    end

    fprintf('  iter %3d/%d  best ITAE = %.6g  mean ITAE = %.6g\n', ...
        iter,maxIter,bestHistory(iter),meanHistory(iter));
end

bestGains = unpackGains(controllerType,alphaPos(:),seedGains);

gwoOut.controllerType = controllerType;
gwoOut.bestVector = alphaPos(:);
gwoOut.bestFitness = alphaScore;
gwoOut.bestGains = bestGains;
gwoOut.fitnessHistory = bestHistory;
gwoOut.meanFitnessHistory = meanHistory;
gwoOut.lb = lb(:);
gwoOut.ub = ub(:);
gwoOut.nWolves = nWolves;
gwoOut.maxIter = maxIter;
gwoOut.dim = dim;
gwoOut.fitnessTs = tFit(2)-tFit(1);
gwoOut.fitnessTfinal = tFit(end);
end

%% ------------------------------------------------------------------------
function gwoCfg = applyGwoDefaults(gwoCfg)
if nargin < 1 || isempty(gwoCfg)
    gwoCfg = struct();
end
if ~isfield(gwoCfg,'nWolves') || isempty(gwoCfg.nWolves)
    gwoCfg.nWolves = 8;
end
if ~isfield(gwoCfg,'maxIter') || isempty(gwoCfg.maxIter)
    gwoCfg.maxIter = 12;
end
if ~isfield(gwoCfg,'rngSeed')
    gwoCfg.rngSeed = 42;
end
if ~isfield(gwoCfg,'tuneTs')
    gwoCfg.tuneTs = 0.01;
end
end

function [cfgFit,tFit,qFit] = makeFitnessGrid(~,t,qRef,cfg,gwoCfg)
cfgFit = cfg;
if isempty(gwoCfg.tuneTs) || gwoCfg.tuneTs == cfg.Ts
    tFit = t;
    qFit = qRef;
    return;
end
cfgFit.Ts = gwoCfg.tuneTs;
tFit = (t(1):cfgFit.Ts:t(end)).';
qFit = interp1(t,qRef,tFit,'previous','extrap');
end

function [lb,ub] = gwo_parameter_bounds(controllerType,seedGains,gwoCfg)
% Joint-wise boxes around the frozen baseline. Orders stay in (0, 2).
if isfield(gwoCfg,'lb') && isfield(gwoCfg,'ub') ...
        && ~isempty(gwoCfg.lb) && ~isempty(gwoCfg.ub)
    lb = gwoCfg.lb(:).';
    ub = gwoCfg.ub(:).';
    return;
end

Kp = seedGains.Kp(:);
Ki = seedGains.Ki(:);
Kd = seedGains.Kd(:);

lbKp = max(1,0.25*Kp);
ubKp = max(lbKp+1,min(400,2.5*Kp));
lbKi = max(0.05,0.25*Ki);
ubKi = max(lbKi+0.05,min(80,2.5*Ki));
lbKd = max(0.05,0.25*Kd);
ubKd = max(lbKd+0.05,min(80,2.5*Kd));

switch controllerType
    case 'PID'
        lb = [lbKp; lbKi; lbKd].';
        ub = [ubKp; ubKi; ubKd].';
    case 'FOPID'
        lbLam = 0.20*ones(6,1);
        ubLam = 1.50*ones(6,1);
        lbMu = 0.20*ones(6,1);
        ubMu = 1.50*ones(6,1);
        lb = [lbKp; lbKi; lbKd; lbLam; lbMu].';
        ub = [ubKp; ubKi; ubKd; ubLam; ubMu].';
    otherwise
        error('gwo_optimize_controller:UnknownType', ...
            'Unknown controller type: %s',controllerType);
end
end

function x = packGains(controllerType,gains)
switch controllerType
    case 'PID'
        x = [gains.Kp(:); gains.Ki(:); gains.Kd(:)];
    case 'FOPID'
        x = [gains.Kp(:); gains.Ki(:); gains.Kd(:); ...
             gains.lambda(:); gains.mu(:)];
    otherwise
        error('gwo_optimize_controller:UnknownType', ...
            'Unknown controller type: %s',controllerType);
end
end

function gains = unpackGains(controllerType,x,template)
gains = template;
x = x(:);
gains.Kp = x(1:6);
gains.Ki = x(7:12);
gains.Kd = x(13:18);
switch controllerType
    case 'PID'
        if isfield(gains,'lambda'), gains = rmfield(gains,'lambda'); end
        if isfield(gains,'mu'), gains = rmfield(gains,'mu'); end
    case 'FOPID'
        gains.lambda = x(19:24);
        gains.mu = x(25:30);
    otherwise
        error('gwo_optimize_controller:UnknownType', ...
            'Unknown controller type: %s',controllerType);
end
end

function J = evaluateItae(plant,t,qRef,controllerType,x,cfg)
gains = unpackGains(controllerType,x,struct( ...
    'Kp',zeros(6,1),'Ki',zeros(6,1),'Kd',zeros(6,1), ...
    'lambda',ones(6,1),'mu',ones(6,1)));
try
    result = simulate_ur5_controller(plant,t,qRef,controllerType,gains,cfg);
    if any(~isfinite(result.q(:))) || any(imag(result.q(:)) ~= 0) ...
            || max(abs(result.q(:))) > 10*pi
        J = inf;
        return;
    end
    J = compute_itae(t,qRef,result.q);
    if ~isfinite(J)
        J = inf;
    end
catch
    J = inf;
end
end

function positions = gwoUpdatePositions(positions,alphaPos,betaPos,deltaPos,a,lb,ub)
[nWolves,dim] = size(positions);

r1 = rand(nWolves,dim);
r2 = rand(nWolves,dim);
A1 = 2*a*r1 - a;
C1 = 2*r2;
X1 = alphaPos - A1.*abs(C1.*alphaPos - positions);

r1 = rand(nWolves,dim);
r2 = rand(nWolves,dim);
A2 = 2*a*r1 - a;
C2 = 2*r2;
X2 = betaPos - A2.*abs(C2.*betaPos - positions);

r1 = rand(nWolves,dim);
r2 = rand(nWolves,dim);
A3 = 2*a*r1 - a;
C3 = 2*r2;
X3 = deltaPos - A3.*abs(C3.*deltaPos - positions);

positions = (X1 + X2 + X3)/3;
positions = min(max(positions,lb),ub);
end

function [aS,aP,bS,bP,dS,dP] = rankLeaders(positions,fitness)
[sortedFit,idx] = sort(fitness,'ascend');
n = numel(idx);
aS = sortedFit(1); aP = positions(idx(1),:);
if n >= 2
    bS = sortedFit(2); bP = positions(idx(2),:);
else
    bS = aS; bP = aP;
end
if n >= 3
    dS = sortedFit(3); dP = positions(idx(3),:);
else
    dS = bS; dP = bP;
end
end

function [aS,aP,bS,bP,dS,dP] = updateLeaders(x,J,aS,aP,bS,bP,dS,dP)
if J < aS
    dS = bS; dP = bP;
    bS = aS; bP = aP;
    aS = J;  aP = x;
elseif J < bS
    dS = bS; dP = bP;
    bS = J;  bP = x;
elseif J < dS
    dS = J;  dP = x;
end
end
