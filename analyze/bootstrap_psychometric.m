% bootstrap & fit psychometric curve
clear all;

subj = 5;
ADAPT = '22.5'; % '45', '22.5'

load(['../data/data_' ADAPT '_sub' num2str(subj) '.mat'], 'adaptor', 'test', 'dtheta', ...
    'nTot_ctrl', 'nRight_ctrl', 'pRight_ctrl', 'nTot_exp', 'nRight_exp', 'pRight_exp'); 

%% fit psychometric to data
PF = @PAL_CumulativeNormal;
vlb = [0];
vub = [inf];
options = optimset('MaxFunEvals', 10000);

param_ctrl = NaN(length(test), 4);
thresh_ctrl = NaN(1, length(test));
pRight_pred_ctrl = NaN(length(test), length(dtheta));
for t = 1:length(test)
    beta_fit = NaN(1,5);
    negLog = NaN(1,5);
    for i = 1:5
        param0 = [i/5];
        [ beta_fit(i), negLog(i), exitflag, output ] = fminsearchbnd(@(beta)negLog_psychometric( beta, dtheta, nRight_ctrl(t, :), nTot_ctrl(t, :), PF ),param0,vlb,vub, options);
    end
    [M, I] = min(negLog);
    param_ctrl(t, :) = [0, beta_fit(I), 0, 0];
    thresh_ctrl(t) = (PF(param_ctrl(t, :), 0.75, 'inverse')-PF(param_ctrl(t, :), 0.25, 'inverse'))/2;
    pRight_pred_ctrl(t,:) = PF(param_ctrl(t, :), dtheta);
end

param_exp = NaN(length(test), 4);
thresh_exp = NaN(1, length(test));
pRight_pred_exp = NaN(length(test), length(dtheta));
for t = 1:length(test)
    beta_fit = NaN(1,5);
    negLog = NaN(1,5);
    for i = 1:5
        param0 = [i/5];
        [ beta_fit(i), negLog(i), exitflag, output ] = fminsearchbnd(@(beta)negLog_psychometric( beta, dtheta, nRight_exp(t, :), nTot_exp(t, :), PF ),param0,vlb,vub, options);
    end
    [M, I] = min(negLog);
    param_exp(t, :) = [0, beta_fit(I), 0, 0];
    thresh_exp(t) = (PF(param_exp(t, :), 0.75, 'inverse')-PF(param_exp(t, :), 0.25, 'inverse'))/2;
    pRight_pred_exp(t,:) = PF(param_exp(t, :), dtheta);
end

save(['boot_psychometric_' ADAPT '_sub' num2str(subj) '.mat'], 'adaptor', 'test', 'dtheta', 'param_ctrl', 'thresh_ctrl', 'pRight_pred_ctrl', 'param_exp', 'thresh_exp', 'pRight_pred_exp')

%% bootstrap
nBoot = 1000;

nTrl = sum(nTot_ctrl(1, :)); 
thresh_ctrl_boot = NaN(nBoot, length(test));
thresh_exp_boot = NaN(nBoot, length(test));
    
for t = 1:length(test)
    nRight_nLeft = reshape([nRight_ctrl(t, :);nTot_ctrl(t, :)-nRight_ctrl(t, :)], 1, 2*length(dtheta));
    edges = 1 + [0, cumsum(nRight_nLeft)];
    for i = 1:nBoot
        Trl_num_boot = randi(nTrl,1,nTrl);
        [nRight_nLeft_boot, edge] = histcounts(Trl_num_boot,edges);
        nRight_nLeft_boot = reshape(nRight_nLeft_boot, 2, length(dtheta));
        nRight_boot = nRight_nLeft_boot(1,:);
        nTot_boot = sum(nRight_nLeft_boot, 1);
        
        beta_fit = NaN(1,5);
        negLog = NaN(1,5);
        for j = 1:5
            param0 = [j/5];
            [ beta_fit(j), negLog(j), exitflag, output ] = fminsearchbnd(@(beta)negLog_psychometric( beta, dtheta, nRight_boot, nTot_boot, PF ),param0,vlb,vub, options);
        end
        [M, I] = min(negLog);
        params_fit = [0, beta_fit(I), 0, 0];
        thresh_ctrl_boot(i,t) = (PF(params_fit, 0.75, 'inverse')-PF(params_fit, 0.25, 'inverse'))/2;
    end

    nRight_nLeft = reshape([nRight_exp(t, :);nTot_exp(t, :)-nRight_exp(t, :)], 1, 2*length(dtheta));
    edges = 1 + [0, cumsum(nRight_nLeft)];
    for i = 1:nBoot
        Trl_num_boot = randi(nTrl,1,nTrl);
        [nRight_nLeft_boot, edge] = histcounts(Trl_num_boot,edges);
        nRight_nLeft_boot = reshape(nRight_nLeft_boot, 2, length(dtheta));
        nRight_boot = nRight_nLeft_boot(1,:);
        nTot_boot = sum(nRight_nLeft_boot, 1);
        
        beta_fit = NaN(1,5);
        negLog = NaN(1,5);
        for j = 1:5
            param0 = [j/5];
            [ beta_fit(j), negLog(j), exitflag, output ] = fminsearchbnd(@(beta)negLog_psychometric( beta, dtheta, nRight_boot, nTot_boot, PF ),param0,vlb,vub, options);
        end
        [M, I] = min(negLog);
        params_fit = [0, beta_fit(I), 0, 0];
        thresh_exp_boot(i,t) = (PF(params_fit, 0.75, 'inverse')-PF(params_fit, 0.25, 'inverse'))/2;
    end
end

thresh_ctrl_confid = prctile(thresh_ctrl_boot, [2.5, 50, 97.5], 1);
thresh_exp_confid = prctile(thresh_exp_boot, [2.5, 50, 97.5], 1);

thresh_ratio_boot = thresh_exp_boot./thresh_ctrl_boot;
thresh_ratio_confid = prctile(thresh_ratio_boot, [2.5, 50, 97.5], 1);
thresh_diff_boot = thresh_exp_boot-thresh_ctrl_boot;
thresh_diff_confid = prctile(thresh_diff_boot, [2.5, 50, 97.5], 1);

save(['boot_psychometric_' ADAPT '_sub' num2str(subj) '.mat'], 'nBoot', 'thresh_ctrl_boot', 'thresh_exp_boot', ...
    'thresh_ctrl_confid', 'thresh_exp_confid', 'thresh_ratio_boot', 'thresh_ratio_confid', 'thresh_diff_boot', 'thresh_diff_confid', '-append')

%% plot
ctrl_color = [0, 113, 188]/255;
exp_color = [216, 82, 24]/255;

figure(1)
set(gcf,'Position',[100, 100, 800, 600]);
hold on
errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_ctrl, thresh_ctrl(1)], ...
    [thresh_ctrl, thresh_ctrl(1)]-[thresh_ctrl_confid(1,:), thresh_ctrl_confid(1,1)], [thresh_ctrl_confid(3,:), thresh_ctrl_confid(3,1)]-[thresh_ctrl, thresh_ctrl(1)], ...
    'o', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', ctrl_color, 'Color', ctrl_color)
plot([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_ctrl, thresh_ctrl(1)], '--', 'LineWidth', 2, 'Color', ctrl_color)
errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_exp, thresh_exp(1)], ...
    [thresh_exp, thresh_exp(1)]-[thresh_exp_confid(1,:), thresh_exp_confid(1,1)], [thresh_exp_confid(3,:), thresh_exp_confid(3,1)]-[thresh_exp, thresh_exp(1)], ...
    'o', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', exp_color, 'Color', exp_color)
plot([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_exp, thresh_exp(1)], '--', 'LineWidth', 2, 'Color', exp_color)
xlim([adaptor(2)-90 adaptor(2)+90])
ylim([0 1.1*max([thresh_ctrl_confid(3,:), thresh_exp_confid(3,:), thresh_ctrl, thresh_exp])])
set(gca,'XTick',adaptor(2)-90:45:adaptor(2)+90)
xlabel('Test orientation (deg)')
ylabel('Threshold (deg)')
set(gca, 'FontSize', 24)


%% function
function negLog = negLog_psychometric( beta, StimLevels, NumPos, OutOfNum, PF )
    epsilon = 1e-2;
    negLog = 0;
    ProbPos = PF( [0,beta,0,0], StimLevels );
    
    for j = 1:length(StimLevels)
        pPos = ProbPos(j);
        pNeg = 1-pPos;
        if pPos < epsilon
            pPos = epsilon;
        end
        if pNeg < epsilon
            pNeg = epsilon;
        end
        negLog = negLog - NumPos(j)*log(pPos) - (OutOfNum(j)-NumPos(j))*log(pNeg);
    end
end



