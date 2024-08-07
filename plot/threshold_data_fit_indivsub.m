% Fig.4
clear all;
ctrl_color = 0.6*[1, 1, 1];
subj_color = [0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];

subjs = 1:5;

figure(1)
set(gcf,'Position',[0, 0, 1000, 900]);

for subj = 1:length(subjs)
    %% 45 deg
    load(['../analyze/boot_psychometric_45_sub' num2str(subj) '.mat'], 'adaptor', 'test', 'thresh_ctrl', 'thresh_exp', 'thresh_ctrl_confid', 'thresh_exp_confid');
    load(['../model/fit_ctrl_4522.5_sub' num2str(subj) '.mat'], 'x', 'thresh_ctrl_pred_45');
    load(['../model/fit_2peak_4522.5_sub' num2str(subj) '.mat'], 'thresh_exp_pred_45');
    
    ymax = max([thresh_ctrl_confid(3,:), thresh_exp_confid(3,:), thresh_ctrl_pred_45, thresh_exp_pred_45]);

    subplot(4, 3, subj+floor((subj-1)/3)*3)
    hold on
    errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_ctrl, thresh_ctrl(1)], ...
        [thresh_ctrl, thresh_ctrl(1)]-[thresh_ctrl_confid(1,:), thresh_ctrl_confid(1,1)], ...
        [thresh_ctrl_confid(3,:), thresh_ctrl_confid(3,1)]-[thresh_ctrl, thresh_ctrl(1)], ...
        'o', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', ctrl_color, 'Color', ctrl_color)
    errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_exp, thresh_exp(1)], ...
        [thresh_exp, thresh_exp(1)]-[thresh_exp_confid(1,:), thresh_exp_confid(1,1)], ...
        [thresh_exp_confid(3,:), thresh_exp_confid(3,1)]-[thresh_exp, thresh_exp(1)], ...
        'o', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', subj_color(subj,:), 'Color', subj_color(subj,:))
    plot(x+adaptor(2), thresh_ctrl_pred_45, 'LineWidth', 2, 'Color', ctrl_color);
    plot(x+adaptor(2), thresh_exp_pred_45, 'LineWidth', 2, 'Color', subj_color(subj,:))

    xlim([adaptor(2)-90-5 adaptor(2)+90+5])
    set(gca,'XTick',adaptor(2)-90:45:adaptor(2)+90)
    xlabel('Test orientation (deg)')
    ylabel('Threshold (deg)')
    set(gca, 'FontSize', 9)

    %% 22.5 deg
    load(['../analyze/boot_psychometric_22.5_sub' num2str(subj) '.mat'], 'adaptor', 'test', 'thresh_ctrl', 'thresh_exp', 'thresh_ctrl_confid', 'thresh_exp_confid');
    load(['../model/fit_ctrl_4522.5_sub' num2str(subj) '.mat'], 'x', 'thresh_ctrl_pred_225');
    load(['../model/fit_2peak_4522.5_sub' num2str(subj) '.mat'], 'thresh_exp_pred_225');
    
    ymax = max([thresh_ctrl_confid(3,:), thresh_exp_confid(3,:), thresh_ctrl_pred_225, thresh_exp_pred_225, ymax]);

    subplot(4, 3, subj+floor((subj-1)/3)*3+3)
    hold on
    errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_ctrl, thresh_ctrl(1)], ...
        [thresh_ctrl, thresh_ctrl(1)]-[thresh_ctrl_confid(1,:), thresh_ctrl_confid(1,1)], ...
        [thresh_ctrl_confid(3,:), thresh_ctrl_confid(3,1)]-[thresh_ctrl, thresh_ctrl(1)], ...
        'o', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', ctrl_color, 'Color', ctrl_color)
    errorbar([-test+adaptor(2), -test(1)-180+adaptor(2)], [thresh_exp, thresh_exp(1)], ...
        [thresh_exp, thresh_exp(1)]-[thresh_exp_confid(1,:), thresh_exp_confid(1,1)], ...
        [thresh_exp_confid(3,:), thresh_exp_confid(3,1)]-[thresh_exp, thresh_exp(1)], ...
        'o', 'MarkerSize', 6, 'LineWidth', 1, 'MarkerEdgeColor', 'None', 'MarkerFaceColor', subj_color(subj,:), 'Color', subj_color(subj,:))
    plot(x+adaptor(2), thresh_ctrl_pred_225, 'LineWidth', 2, 'Color', ctrl_color);
    plot(x+adaptor(2), thresh_exp_pred_225, 'LineWidth', 2, 'Color', subj_color(subj,:))

    xlim([adaptor(2)-90-5 adaptor(2)+90+5])
    ylim([0 1.1*ymax])
    set(gca,'XTick',adaptor(2)-90:45:adaptor(2)+90)
    xlabel('Test orientation (deg)')
    ylabel('Threshold (deg)')
    set(gca, 'FontSize', 9)
    
    subplot(4, 3, subj+floor((subj-1)/3)*3)
    ylim([0 1.1*ymax])

end



% figure(1)
% set(gcf,'Position',[1500, 0, 1000, 900]);
% % saveas(1, ['thresh_data_fit_45_22.5_indivsub.fig'])
% % saveas(1, ['thresh_data_fit_45_22.5_indivsub.png'])
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, ['thresh_data_fit_45_22.5_indivsub.pdf'], '-dpdf')
% % close(1)



