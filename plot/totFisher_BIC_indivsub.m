% Fig.5b&c
clear all;

subjs = 1:5;
n_subj = length(subjs);

FIT_VERs = {'1peak', '2peak', '2peakF', '2peakK'};
n_ver = length(FIT_VERs);
vers_label = {'1-peak', ...
    '2-peak', ...
    '2-peak + Fisher', ...
    '2-peak + kernel'};

subj_color = [0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];
font = 16;

%% sqrt Fisher in ctrl & 2peakF
figure(1)
set(gcf,'Position',[0, 0, 960, 400]);
subplot(1,2,1)
hold on
ver_num = '06';
refline = plot([20 100], [20 100], 'k--', 'LineWidth', 1);
for j = 1:n_subj
    subj = subjs(j);
    load(['../model/fit_ctrl_4522.5_sub' num2str(subj) '.mat'], 'totFisher')
    load(['../model/fit_2peakF_4522.5_sub' num2str(subj) '.mat'], 'totFisher1')
    load(['../model/fit_boot_ctrl_sub' num2str(subj) '_boot200.mat'], 'Fisher0_confid')
    load(['../model/fit_boot_2peakF_sub' num2str(subj) '_boot200.mat'], 'Fisher1_confid')
    
    errorbar(totFisher, totFisher1, totFisher1-Fisher1_confid(1), Fisher1_confid(3)-totFisher1, totFisher-Fisher0_confid(1), Fisher0_confid(3)-totFisher, ...
        'o-', 'MarkerSize', 8, 'LineWidth', 2, 'Color', subj_color(j,:), 'MarkerFaceColor', subj_color(j,:))
end
xlim([20 100])
ylim([20 100])
set(gca,'XTick',20:20:100, 'YTick',20:20:100)
xlabel('Total Fisher - control')
ylabel('Total Fisher - oblique')
set(gca, 'FontSize', font)

%% BIC
subplot(1,2,2)
hold on
refline = plot([0.5 n_ver+0.5], [0 0], 'k--', 'LineWidth', 1);
for j = 1:n_subj
    subj = subjs(j);
    BIC_exp_combo = NaN(1, n_ver);
    for i = 1:n_ver
        load(['../model/fit_' FIT_VERs{i} '_4522.5_sub' num2str(subj) '.mat'], 'BIC_exp')
        BIC_exp_combo(i) = BIC_exp; 
    end
    relBIC_exp_combo = BIC_exp_combo - BIC_exp_combo(2);
    
    plot(1:n_ver, relBIC_exp_combo, 'o-', 'MarkerSize',8, 'LineWidth',2, 'Color',subj_color(j,:), 'MarkerFaceColor',subj_color(j,:), 'MarkerEdgeColor','none')
end
xlim([0.5 n_ver+0.5])
xticks(1:n_ver)
xticklabels(vers_label)
ylabel('Relative BIC')
set(gca, 'FontSize', font)


% figure(1)
% set(gcf,'Position',[0, 0, 960, 400]);
% % saveas(1, ['totFisher_ctrl_2peakF_BIC_indivsub.fig'])
% % saveas(1, ['totFisher_ctrl_2peakF_BIC_indivsub.png'])
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, ['totFisher_ctrl_2peakF_BIC_indivsub.pdf'], '-dpdf')

