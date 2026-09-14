function [metrics,summary] = make_step_metrics(t,qRef,qPID,qFOPID,tauPID,tauFOPID)
%MAKE_STEP_METRICS Calculate step-response metrics for PID and/or FOPID.
%
% If both qPID and qFOPID are supplied, the returned table contains 12 rows
% (6 joints per controller), matching the paper's six-joint comparison.

if nargin < 3, qPID = []; end
if nargin < 4, qFOPID = []; end
if nargin < 5, tauPID = []; end
if nargin < 6, tauFOPID = []; end

rows = {};

if ~isempty(qPID)
    rows = [rows; controllerRows(t,qRef,qPID,'PID')]; %#ok<AGROW>
end
if ~isempty(qFOPID)
    rows = [rows; controllerRows(t,qRef,qFOPID,'FOPID')]; %#ok<AGROW>
end

if isempty(rows)
    error('At least one controller response must be supplied.');
end

metrics = cell2table(rows,'VariableNames', ...
    {'Joint','Controller','Overshoot_percent','SettlingTime_s', ...
     'PeakTime_s','Target_rad'});

controllers = unique(metrics.Controller,'stable');
summary = table('Size',[numel(controllers) 5], ...
    'VariableTypes',{'string','double','double','double','double'}, ...
    'VariableNames',{'Controller','AverageOvershoot_percent', ...
                     'AverageSettlingTime_s','AveragePeakTime_s', ...
                     'AverageAbsTorque_Nm'});

for c = 1:numel(controllers)
    idx = strcmp(metrics.Controller,controllers{c});
    summary.Controller(c) = string(controllers{c});
    summary.AverageOvershoot_percent(c) = mean(metrics.Overshoot_percent(idx),'omitnan');
    summary.AverageSettlingTime_s(c) = mean(metrics.SettlingTime_s(idx),'omitnan');
    summary.AveragePeakTime_s(c) = mean(metrics.PeakTime_s(idx),'omitnan');
    if strcmp(controllers{c},'PID') && ~isempty(tauPID)
        summary.AverageAbsTorque_Nm(c) = mean(sum(abs(tauPID),1));
    elseif strcmp(controllers{c},'FOPID') && ~isempty(tauFOPID)
        summary.AverageAbsTorque_Nm(c) = mean(sum(abs(tauFOPID),1));
    else
        summary.AverageAbsTorque_Nm(c) = NaN;
    end
end
end

function rows = controllerRows(t,qRef,q,name)
nJ = size(q,2);
rows = cell(nJ,6);
for j = 1:nJ
    y = q(:,j);
    target = qRef(end,j);

    try
        s = stepinfo(y,t,target,'SettlingTimeThreshold',0.02);
        overshoot = s.Overshoot;
        settling = s.SettlingTime;
        peakTime = s.PeakTime;
    catch
        peak = max(y);
        if abs(target) > 1e-12
            overshoot = max(0,(peak-target)/abs(target)*100);
        else
            overshoot = NaN;
        end
        settling = NaN;
        peakTime = NaN;
    end

    rows(j,:) = {j,name,overshoot,settling,peakTime,target};
end
end
