% 1peak: increased coding accuracy at adaptor only
% fit to 45 & 22.5 oblique conditions jointly, with parallel computing
clear all;
epsilon = 2e-2;
res = 0.5;
FIT_VER = '1peak';

for subj = 5 %1:5

load(['../data/data_22.5_sub' num2str(subj) '.mat'])
load(['../analyze/boot_psychometric_22.5_sub' num2str(subj) '.mat'])

adaptor_225 = adaptor;
test_225 = test;
dtheta_225 = dtheta;
% nTot_ctrl_225 = nTot_ctrl;
% nRight_ctrl_225 = nRight_ctrl;
% pRight_ctrl_225 = pRight_ctrl;
nTot_exp_225 = nTot_exp;
nRight_exp_225 = nRight_exp;
pRight_exp_225 = pRight_exp;
% thresh_ctrl_225 = thresh_ctrl;
% thresh_ctrl_confid_225 = thresh_ctrl_confid;
thresh_exp_225 = thresh_exp;
thresh_exp_confid_225 = thresh_exp_confid;
% pRight_pred_ctrl_225 = pRight_pred_ctrl;
pRight_pred_exp_225 = pRight_pred_exp;

load(['../data/data_45_sub' num2str(subj) '.mat'])
load(['../analyze/boot_psychometric_45_sub' num2str(subj) '.mat'])

adaptor_45 = adaptor;
test_45 = test;
dtheta_45 = dtheta;
% nTot_ctrl_45 = nTot_ctrl;
% nRight_ctrl_45 = nRight_ctrl;
% pRight_ctrl_45 = pRight_ctrl;
nTot_exp_45 = nTot_exp;
nRight_exp_45 = nRight_exp;
pRight_exp_45 = pRight_exp;
% thresh_ctrl_45 = thresh_ctrl;
% thresh_ctrl_confid_45 = thresh_ctrl_confid;
thresh_exp_45 = thresh_exp;
thresh_exp_confid_45 = thresh_exp_confid;
% pRight_pred_ctrl_45 = pRight_pred_ctrl;
pRight_pred_exp_45 = pRight_pred_exp;

load(['fit_ctrl_4522.5_sub' num2str(subj) '.mat'], 'theta_o', 'kappa_e', 'k_o', 'kappa_o', 'kappa_i', 'x');

%% adaptation
theta = [0];
params0 = [0.1, 30];
vlb = [-inf, 0.01]; 
vub = [inf, 700];

currPool = gcp('nocreate');
if isempty(currPool)
    parpool(6)
end

negLog_adapt_prior( params0, k_o, theta_o, kappa_o, theta, kappa_i, kappa_e, adaptor_45(1), test_45, dtheta_45, nTot_exp_45, nRight_exp_45, adaptor_225(1), test_225, dtheta_225, nTot_exp_225, nRight_exp_225, res, epsilon )
profile on
% [fit_params, negLog_exp] = fminsearchbnd(@(params)negLog_adapt_prior( params, k_o, theta_o, kappa_o, theta, kappa_i, kappa_e, adaptor_45(1), test_45, dtheta_45, nTot_exp_45, nRight_exp_45, adaptor_225(1), test_225, dtheta_225, nTot_exp_225, nRight_exp_225, res, epsilon ),params0,vlb,vub)
[fit_params, negLog_exp] = fmincon(@(params)negLog_adapt_prior( params, k_o, theta_o, kappa_o, theta, kappa_i, kappa_e, adaptor_45(1), test_45, dtheta_45, nTot_exp_45, nRight_exp_45, adaptor_225(1), test_225, dtheta_225, nTot_exp_225, nRight_exp_225, res, epsilon ),...
    params0,[],[],[],[],vlb,vub,@prior_con)
profile off
BIC_exp = 2*negLog_exp + length(vlb)*log(sum(nTot_exp_45,'all')+sum(nTot_exp_225,'all'));

k = fit_params(1);
kappa = fit_params(2);

fit_kernel = sum_n_vmpdf_180( x, k, theta, kappa );

[ ~, ~, ~, thresh_exp_pred_45 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor_45(2), adaptor_45(2)+x );
[ ~, ~, ~, thresh_exp_pred_225 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor_225(2), adaptor_225(2)+x );

save(['fit_' FIT_VER '_4522.5_sub' num2str(subj) '.mat'], 'theta', 'k', 'kappa', 'negLog_exp', 'BIC_exp', 'x', 'fit_kernel', 'thresh_exp_pred_45', 'thresh_exp_pred_225');



x2 = -90:res:90-res; % detheta grid

[ x1, probCmpCw ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor_45(1), adaptor_45(1)+test_45 );
probCmpCw_exp_45 = NaN(size(probCmpCw));
for i = 1:length(test_45)
    probCmpCw_exp_45(i,:) = circshift(probCmpCw(i,:), -round((adaptor_45(1)+90+test_45(i))/res));
end

[ x1, probCmpCw ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor_225(1), adaptor_225(1)+test_225 );
probCmpCw_exp_225 = NaN(size(probCmpCw));
for i = 1:length(test_225)
    probCmpCw_exp_225(i,:) = circshift(probCmpCw(i,:), -round((adaptor_225(1)+90+test_225(i))/res));
end

save(['fit_' FIT_VER '_4522.5_sub' num2str(subj) '.mat'], 'x2', 'probCmpCw_exp_45', 'probCmpCw_exp_225', '-append');

%% plot
ctrl_color = [0, 113, 188]/255;
exp_color = [216, 82, 24]/255;

for i = 1:length(test_45)
    figure(i)
    hold on
    scatter(dtheta_45, pRight_exp_45(i, :), nTot_exp_45(i, :)*10+1, 'MarkerEdgeColor', exp_color);
    plot(dtheta_45, pRight_pred_exp_45(i, :), '--', 'LineWidth', 2, 'Color', exp_color);
    plot(x2, probCmpCw_exp_45(i,:), '-', 'LineWidth', 2, 'Color', exp_color);
    hold off
    xlim([dtheta_45(1), dtheta_45(end)])
    ylim([0, 1.01])
end
for i = 1:length(test_225)
    figure(length(test_45)+i)
    hold on
    scatter(dtheta_225, pRight_exp_225(i, :), nTot_exp_225(i, :)*10+1, 'MarkerEdgeColor', exp_color);
    plot(dtheta_45, pRight_pred_exp_225(i, :), '--', 'LineWidth', 2, 'Color', exp_color);
    plot(x2, probCmpCw_exp_225(i,:), '-', 'LineWidth', 2, 'Color', exp_color);
    hold off
    xlim([dtheta_225(1), dtheta_225(end)])
    ylim([0, 1.01])
end

figure(length(test_45)+length(test_225)+1)
set(gcf,'Position',[0, 0, 1200, 450]);
ymax = ceil(max([thresh_exp_confid_45(3,:), thresh_exp_confid_225(3,:), thresh_exp_pred_45, thresh_exp_pred_225]));

subplot(1,2,1)
hold on
errorbar([-test_45+adaptor_45(2), -test_45(1)-180+adaptor_45(2)], [thresh_exp_45, thresh_exp_45(1)], ...
    [thresh_exp_45, thresh_exp_45(1)]-[thresh_exp_confid_45(1,:), thresh_exp_confid_45(1,1)], [thresh_exp_confid_45(3,:), thresh_exp_confid_45(3,1)]-[thresh_exp_45, thresh_exp_45(1)], ...
    'o', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', exp_color, 'Color', exp_color)
plot(x+adaptor_45(2), thresh_exp_pred_45, 'LineWidth', 2, 'Color', exp_color)
xlim([adaptor_45(2)-90 adaptor_45(2)+90])
ylim([0 ymax])
set(gca,'XTick',adaptor_45(2)-90:45:adaptor_45(2)+90, 'YTick',0:2:ymax)
xlabel('Test orientation (deg)')
ylabel('Threshold (deg)')
set(gca, 'FontSize', 18)

subplot(1,2,2)
hold on
errorbar([-test_225+adaptor_225(2), -test_225(1)-180+adaptor_225(2)], [thresh_exp_225, thresh_exp_225(1)], ...
    [thresh_exp_225, thresh_exp_225(1)]-[thresh_exp_confid_225(1,:), thresh_exp_confid_225(1,1)], [thresh_exp_confid_225(3,:), thresh_exp_confid_225(3,1)]-[thresh_exp_225, thresh_exp_225(1)], ...
    'o', 'MarkerSize', 12, 'LineWidth', 2, 'MarkerEdgeColor', exp_color, 'Color', exp_color)
plot(x+adaptor_225(2), thresh_exp_pred_225, 'LineWidth', 2, 'Color', exp_color)
xlim([adaptor_225(2)-90 adaptor_225(2)+90])
ylim([0 ymax])
set(gca,'XTick',adaptor_225(2)-90:45:adaptor_225(2)+90, 'YTick',0:2:ymax)
xlabel('Test orientation (deg)')
ylabel('Threshold (deg)')
set(gca, 'FontSize', 18)

end
%%
function negLog = negLog_overall( params, theta_o, kappa_e, adaptor1, test1, dtheta1, nTot1, nRight1, adaptor2, test2, dtheta2, nTot2, nRight2, res, epsilon )
    k_o = params(1);
    kappa_o = params(2);
    kappa_i = params(3);
    
    [ x1, probCmpCw01 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], 0, 0, 5, kappa_i, kappa_e, res, adaptor1, adaptor1+test1 ); % ProbRefCw: absolute orientation of ref, not dtheta_pool
    [ x1, probCmpCw02 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], 0, 0, 5, kappa_i, kappa_e, res, adaptor2, adaptor2+test2 ); % ProbRefCw: absolute orientation of ref, not dtheta_pool

    x2 = -90:res:90-res; % dtheta grid
    probCmpCw1 = NaN(size(probCmpCw01));
    probCmpCw2 = NaN(size(probCmpCw02));
    for i = 1:length(test1)
        probCmpCw1(i,:) = circshift(probCmpCw01(i,:),-round((adaptor1+90+test1(i))/res)); % re-centered at dtheta=0
    end
    for i = 1:length(test2)
        probCmpCw2(i,:) = circshift(probCmpCw02(i,:),-round((adaptor2+90+test2(i))/res));
    end
    
    negLog = 0;
    for i = 1:length(test1)
        for j = 1:length(dtheta1)
            pCW = interp1(x2, probCmpCw1(i,:), dtheta1(j));
            pCCW = 1-pCW;
            if pCW < epsilon
                pCW = epsilon;
            end
            if pCCW < epsilon
                pCCW = epsilon;
            end
            negLog = negLog - nRight1(i,j)*log(pCW) - (nTot1(i,j)-nRight1(i,j))*log(pCCW);
        end
    end
    for i = 1:length(test2)
        for j = 1:length(dtheta2)
            pCW = interp1(x2, probCmpCw2(i,:), dtheta2(j));
            pCCW = 1-pCW;
            if pCW < epsilon
                pCW = epsilon;
            end
            if pCCW < epsilon
                pCCW = epsilon;
            end
            negLog = negLog - nRight2(i,j)*log(pCW) - (nTot2(i,j)-nRight2(i,j))*log(pCCW);
        end
    end
end

function negLog = negLog_adapt_prior( params, k_o, theta_o, kappa_o, theta, kappa_i, kappa_e, adaptor1, test1, dtheta1, nTot1, nRight1, adaptor2, test2, dtheta2, nTot2, nRight2, res, epsilon )
    k = params(1:length(theta));
    kappa = params(length(theta)+1:length(theta)*2);
    
    [ x1, probCmpCw01 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor1, adaptor1+test1 ); % ProbRefCw: absolute orientation of ref, not dtheta_pool
    [ x1, probCmpCw02 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor2, adaptor2+test2 ); % ProbRefCw: absolute orientation of ref, not dtheta_pool

    x2 = -90:res:90-res;
    probCmpCw1 = NaN(size(probCmpCw01));
    probCmpCw2 = NaN(size(probCmpCw02));
    for i = 1:length(test1)
        probCmpCw1(i,:) = circshift(probCmpCw01(i,:),-round((adaptor1+90+test1(i))/res));
    end
    for i = 1:length(test2)
        probCmpCw2(i,:) = circshift(probCmpCw02(i,:),-round((adaptor2+90+test2(i))/res));
    end
    
    negLog = 0;
    for i = 1:length(test1)
        for j = 1:length(dtheta1)
            pCW = interp1(x2, probCmpCw1(i,:), dtheta1(j));
            pCCW = 1-pCW;
            if pCW < epsilon
                pCW = epsilon;
            end
            if pCCW < epsilon
                pCCW = epsilon;
            end
            negLog = negLog - nRight1(i,j)*log(pCW) - (nTot1(i,j)-nRight1(i,j))*log(pCCW);
        end
    end
    for i = 1:length(test2)
        for j = 1:length(dtheta2)
            pCW = interp1(x2, probCmpCw2(i,:), dtheta2(j));
            pCCW = 1-pCW;
            if pCW < epsilon
                pCW = epsilon;
            end
            if pCCW < epsilon
                pCCW = epsilon;
            end
            negLog = negLog - nRight2(i,j)*log(pCW) - (nTot2(i,j)-nRight2(i,j))*log(pCCW);
        end
    end
end

function negLog = negLog_overall_one( params, theta_o, kappa_e, adaptor, test, dtheta_pool, nTot, nRight, res, epsilon )
    k_o = params(1);
    kappa_o = params(2);
    kappa_i = params(3);
    
    [ x1, probCmpCw0 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], 0, 0, 5, kappa_i, kappa_e, res, adaptor, adaptor+test ); % ProbRefCw: absolute orientation of ref, not dtheta_pool

    x2 = -90:res:90-res;
    probCmpCw = NaN(size(probCmpCw0));
    for i = 1:length(test)
        probCmpCw(i,:) = circshift(probCmpCw0(i,:),-round((adaptor+90+test(i))/res));
    end
    
    negLog = 0;
    for i = 1:length(test)
        for j = 1:length(dtheta_pool)
            pCW = interp1(x2, probCmpCw(i,:), dtheta_pool(j));
            pCCW = 1-pCW;
            if pCW < epsilon
                pCW = epsilon;
            end
            if pCCW < epsilon
                pCCW = epsilon;
            end
            negLog = negLog - nRight(i,j)*log(pCW) - (nTot(i,j)-nRight(i,j))*log(pCCW);
        end
    end
end

function negLog = negLog_adapt_prior_one( params, k_o, theta_o, kappa_o, theta, kappa_i, kappa_e, adaptor, test, dtheta_pool, nTot, nRight, res, epsilon )
    k = params(1:length(theta));
    kappa = params(length(theta)+1:length(theta)*2);
    
    [ x1, probCmpCw0 ] = ECAdapt_2AFC_par( [k_o k_o], theta_o, [kappa_o kappa_o], k, theta, kappa, kappa_i, kappa_e, res, adaptor, adaptor+test ); % ProbRefCw: absolute orientation of ref, not dtheta_pool

    x2 = -90:res:90-res;
    probCmpCw = NaN(size(probCmpCw0));
    for i = 1:length(test)
        probCmpCw(i,:) = circshift(probCmpCw0(i,:),-round((adaptor+90+test(i))/res));
    end
    
    negLog = 0;
    for i = 1:length(test)
        for j = 1:length(dtheta_pool)
            pCW = interp1(x2, probCmpCw(i,:), dtheta_pool(j));
            pCCW = 1-pCW;
            if pCW < epsilon
                pCW = epsilon;
            end
            if pCCW < epsilon
                pCCW = epsilon;
            end
            negLog = negLog - nRight(i,j)*log(pCW) - (nTot(i,j)-nRight(i,j))*log(pCCW);
        end
    end
end

function [minprior, ceq] = prior_con( params )
    k = params(1);
    kappa = params(2);
    minprior = 1e-8 - min(sum_n_vmpdf_180( [0,90], k, 0, kappa ));
    ceq = 0;
end
