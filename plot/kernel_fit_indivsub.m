% Fig.5a
clear all;
subj_color = [0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880;...
    0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840];

subjs = 1:5;

figure(1)
set(gcf,'Position',[00, 00, 960, 600]);

for subj = 1:length(subjs)
load(['../model/fit_2peak_4522.5_sub' num2str(subj) '.mat'], 'x', 'fit_kernel');
load(['../model/fit_2peakK_4522.5_sub' num2str(subj) '.mat'], 'fit_kernel_45', 'fit_kernel_225');

subplot(2,3,subj)
hold on
plot(x, fit_kernel, 'LineWidth', 2, 'Color', subj_color(subj,:));
plot(x, fit_kernel_45, '--', 'LineWidth', 2, 'Color', subj_color(subj,:));
plot(x, fit_kernel_225, ':', 'LineWidth', 2, 'Color', subj_color(subj,:));

xlim([-90 90])
ylim([0 0.02])
set(gca,'XTick',-90:45:90)
xlabel('Orientation relative to adaptor (deg)')
ylabel('normalized sqrt Fisher')
set(gca, 'FontSize', 12)

end



% figure(1)
% set(gcf,'Position',[00, 00, 960, 600]);
% % saveas(1, ['kernel_fit_2peak_2peakK_indivsub.fig'])
% % saveas(1, ['kernel_fit_2peak_2peakK_indivsub.png'])
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% print(fig, ['kernel_fit_2peak_2peakK_indivsub.pdf'], '-dpdf')
% % close(1)



