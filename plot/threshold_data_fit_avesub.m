% Fig.3
clear all;
ctrl_color = [0, 113, 188]/255;
exp_color = [216, 82, 24]/255;

subjs = 1:5;
nTest = 8;
nBoot = 1000;

figure(1)
set(gcf,'Position',[00, 00, 1000, 400]);

%% 45 deg
thresh_ctrl_combo = NaN(length(subjs), nTest);
thresh_exp_combo = NaN(length(subjs), nTest);
thresh_ctrl_boot_combo = NaN(nBoot, nTest, length(subjs));
thresh_exp_boot_combo = NaN(nBoot, nTest, length(subjs));
thresh_ctrl_pred_combo = NaN(length(subjs), 181);
thresh_exp_pred_combo = NaN(length(subjs), 181);

for subj = 1:length(subjs)
    load(['../analyze/boot_psychometric_45_sub' num2str(subj) '.mat'], 'adaptor', 'test', 'thresh_ctrl', 'thresh_exp', 'thresh_ctrl_boot', 'thresh_exp_boot');
    load(['../model/fit_ctrl_4522.5_sub' num2str(subj) '.mat'], 'x', 'thresh_ctrl_pred_45');
    load(['../model/fit_2peak_4522.5_sub' num2str(subj) '.mat'], 'thresh_exp_pred_45');
    
    if subj == 1
        thresh_ctrl_combo(subj,:) = thresh_ctrl([1,2,3,4,6,8,9,10]);
        thresh_exp_combo(subj,:) = thresh_exp([1,2,3,4,6,8,9,10]);
        thresh_ctrl_boot_combo(:,:,subj) = thresh_ctrl_boot(:,[1,2,3,4,6,8,9,10]);
        thresh_exp_boot_combo(:,:,subj) = thresh_exp_boot(:,[1,2,3,4,6,8,9,10]);
    else
        thresh_ctrl_combo(subj,:) = thresh_ctrl;
        thresh_exp_combo(subj,:) = thresh_exp;
        thresh_ctrl_boot_combo(:,:,subj) = thresh_ctrl_boot;
        thresh_exp_boot_combo(:,:,subj) = thresh_exp_boot;
    end
    thresh_ctrl_pred_combo(subj,:) = thresh_ctrl_pred_45;
    thresh_exp_pred_combo(subj,:) = thresh_exp_pred_45;
end
thresh_ctrl_mean = mean(thresh_ctrl_combo,1);
thresh_exp_mean = mean(thresh_exp_combo,1);
thresh_ctrl_boot_mean = mean(thresh_ctrl_boot_combo,3);
thresh_exp_boot_mean = mean(thresh_exp_boot_combo,3);
thresh_ctrl_pred_mean = mean(thresh_ctrl_pred_combo,1);
thresh_exp_pred_mean = mean(thresh_exp_pred_combo,1);

thresh_ctrl_mean_confid = prctile(thresh_ctrl_boot_mean, [2.5, 50, 97.5], 1);
thresh_exp_mean_confid = prctile(thresh_exp_boot_mean, [2.5, 50, 97.5], 1);

subplot(1,2,1)
hold on
plot(x+adaptor(2), thresh_ctrl_pred_mean, 'LineWidth', 2, 'Color', ctrl_color);
plot(x+adaptor(2), thresh_exp_pred_mean, 'LineWidth', 2, 'Color', exp_color);
errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_ctrl_mean, thresh_ctrl_mean(1)], ...
    [thresh_ctrl_mean, thresh_ctrl_mean(1)]-[thresh_ctrl_mean_confid(1,:), thresh_ctrl_mean_confid(1,1)], ...
    [thresh_ctrl_mean_confid(3,:), thresh_ctrl_mean_confid(3,1)]-[thresh_ctrl_mean, thresh_ctrl_mean(1)], ...
    'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', ctrl_color, 'Color', ctrl_color)
errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_exp_mean, thresh_exp_mean(1)], ...
    [thresh_exp_mean, thresh_exp_mean(1)]-[thresh_exp_mean_confid(1,:), thresh_exp_mean_confid(1,1)], ...
    [thresh_exp_mean_confid(3,:), thresh_exp_mean_confid(3,1)]-[thresh_exp_mean, thresh_exp_mean(1)], ...
    'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', exp_color, 'Color', exp_color)

xlim([adaptor(2)-90-5 adaptor(2)+90+5])
ylim([0 8])
set(gca,'XTick',adaptor(2)-90:45:adaptor(2)+90, 'YTick',0:2:8)
xlabel('Test orientation (deg)')
ylabel('Discrimination threshold (deg)')
set(gca, 'FontSize', 18)

%% 22.5 deg
thresh_ctrl_combo = NaN(length(subjs), nTest);
thresh_exp_combo = NaN(length(subjs), nTest);
thresh_ctrl_boot_combo = NaN(nBoot, nTest, length(subjs));
thresh_exp_boot_combo = NaN(nBoot, nTest, length(subjs));
thresh_ctrl_pred_combo = NaN(length(subjs), 181);
thresh_exp_pred_combo = NaN(length(subjs), 181);

for subj = 1:length(subjs)
    load(['../analyze/boot_psychometric_22.5_sub' num2str(subj) '.mat'], 'adaptor', 'test', 'thresh_ctrl', 'thresh_exp', 'thresh_ctrl_boot', 'thresh_exp_boot');
    load(['../model/fit_ctrl_4522.5_sub' num2str(subj) '.mat'], 'x', 'thresh_ctrl_pred_225');
    load(['../model/fit_2peak_4522.5_sub' num2str(subj) '.mat'], 'thresh_exp_pred_225');
    
    thresh_ctrl_combo(subj,:) = thresh_ctrl;
    thresh_exp_combo(subj,:) = thresh_exp;
    thresh_ctrl_boot_combo(:,:,subj) = thresh_ctrl_boot;
    thresh_exp_boot_combo(:,:,subj) = thresh_exp_boot;
    thresh_ctrl_pred_combo(subj,:) = thresh_ctrl_pred_225;
    thresh_exp_pred_combo(subj,:) = thresh_exp_pred_225;
end
thresh_ctrl_mean = mean(thresh_ctrl_combo,1);
thresh_exp_mean = mean(thresh_exp_combo,1);
thresh_ctrl_boot_mean = mean(thresh_ctrl_boot_combo,3);
thresh_exp_boot_mean = mean(thresh_exp_boot_combo,3);
thresh_ctrl_pred_mean = mean(thresh_ctrl_pred_combo,1);
thresh_exp_pred_mean = mean(thresh_exp_pred_combo,1);

thresh_ctrl_mean_confid = prctile(thresh_ctrl_boot_mean, [2.5, 50, 97.5], 1);
thresh_exp_mean_confid = prctile(thresh_exp_boot_mean, [2.5, 50, 97.5], 1);

subplot(1,2,2)
hold on
plot(x+adaptor(2), thresh_ctrl_pred_mean, 'LineWidth', 2, 'Color', ctrl_color);
plot(x+adaptor(2), thresh_exp_pred_mean, 'LineWidth', 2, 'Color', exp_color);
errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_ctrl_mean, thresh_ctrl_mean(1)], ...
    [thresh_ctrl_mean, thresh_ctrl_mean(1)]-[thresh_ctrl_mean_confid(1,:), thresh_ctrl_mean_confid(1,1)], ...
    [thresh_ctrl_mean_confid(3,:), thresh_ctrl_mean_confid(3,1)]-[thresh_ctrl_mean, thresh_ctrl_mean(1)], ...
    'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', ctrl_color, 'Color', ctrl_color)
errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_exp_mean, thresh_exp_mean(1)], ...
    [thresh_exp_mean, thresh_exp_mean(1)]-[thresh_exp_mean_confid(1,:), thresh_exp_mean_confid(1,1)], ...
    [thresh_exp_mean_confid(3,:), thresh_exp_mean_confid(3,1)]-[thresh_exp_mean, thresh_exp_mean(1)], ...
    'o', 'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', exp_color, 'Color', exp_color)

xlim([adaptor(2)-90-5 adaptor(2)+90+5])
ylim([0 8])
set(gca,'XTick',adaptor(2)-90:45:adaptor(2)+90, 'YTick',0:2:8)
xlabel('Test orientation (deg)')
ylabel('Discrimination threshold (deg)')
set(gca, 'FontSize', 18)



% figure(1)
% set(gcf,'Position',[00, 00, 1000, 400]);
% % saveas(1, ['thresh_data_fit_45_22.5_avesub.fig'])
% % saveas(1, ['thresh_data_fit_45_22.5_avesub.png'])
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, ['thresh_data_fit_45_22.5_avesub.pdf'], '-dpdf')
% % close(1)


