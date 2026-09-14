function [metrics,summary] = make_sine_metrics(~,qRef,pidResult,fopidResult)
%MAKE_SINE_METRICS Compute sine-response metrics used by the paper's Table 4.
controllers = {'PID','FOPID'};
results = {pidResult,fopidResult};
rows = cell(12,8); r = 0;
for c = 1:numel(controllers)
    for j = 1:6
        r = r + 1;
        e = qRef(:,j) - results{c}.q(:,j);
        tau = results{c}.tau(:,j);
        rows(r,:) = {j,controllers{c},mean(e.^2),sqrt(mean(e.^2)), ...
            mean(abs(e)),max(abs(e)),sum(abs(tau)),max(abs(tau))};
    end
end
metrics = cell2table(rows,'VariableNames', ...
    {'Joint','Controller','MSE_rad2','RMSE_rad','MeanAbsError_rad', ...
     'MaxAbsError_rad','MeanAbsTorque_Nm','MaxAbsTorque_Nm'});

summary = table('Size',[2 5], ...
    'VariableTypes',{'string','double','double','double','double'}, ...
    'VariableNames',{'Controller','AverageMSE_rad2','AverageRMSE_rad', ...
                     'AverageMeanAbsTorque_Nm','AverageMaxAbsTorque_Nm'});
for c = 1:numel(controllers)
    idx = strcmp(metrics.Controller,controllers{c});
    summary.Controller(c) = string(controllers{c});
    summary.AverageMSE_rad2(c) = mean(metrics.MSE_rad2(idx));
    summary.AverageRMSE_rad(c) = mean(metrics.RMSE_rad(idx));
    summary.AverageMeanAbsTorque_Nm(c) = mean(metrics.MeanAbsTorque_Nm(idx));
    summary.AverageMaxAbsTorque_Nm(c) = mean(metrics.MaxAbsTorque_Nm(idx));
end
end
