function fn_plotAttractorComparison(PR_norm, PR_norm_rr, win_labels, ...
    align_ER_fixed, align_ER_rr, chance_lvl, ...
    ev_recur_density, re_recur_density, cross_recur_density, recording_ID, figuresDir)
% FN_PLOTATTRACTORCOMPARISON  Three-panel comparison of fixed-window vs
% attractor(recurrence)-defined epochs: dimensionality (PR/N), Evoked-
% Recovery subspace alignment, and recurrence density.
%
%   FN_PLOTATTRACTORCOMPARISON(PR_norm, PR_norm_rr, win_labels, ...
%       align_ER_fixed, align_ER_rr, chance_lvl, ...
%       ev_recur_density, re_recur_density, cross_recur_density, recording_ID, figuresDir)

fig = figure('Name','Fixed vs RR-defined','Position',[100 100 1000 420]);

subplot(1,3,1);
data_pr = [PR_norm; PR_norm_rr];
bar(data_pr','grouped','EdgeColor','none');
set(gca,'XTick',1:3,'XTickLabel',win_labels,'FontSize',10);
ylabel('PR/N'); title('Dimensionality (PR/N)','FontWeight','normal');
legend({'Fixed','RR-defined'},'Location','best'); ylim([0 max(data_pr(:))*1.3+0.01]);

subplot(1,3,2);
bar([align_ER_fixed, align_ER_rr],'FaceColor','flat','EdgeColor','none', ...
    'CData',[0.30 0.50 0.80; 0.80 0.30 0.30]); hold on;
yline(chance_lvl,'r--','LineWidth',1.5,'Label',sprintf('Chance=%.3f',chance_lvl));
set(gca,'XTick',1:2,'XTickLabel',{'Fixed','RR-defined'});
ylabel('Alignment E->R'); title('Subspace alignment','FontWeight','normal'); ylim([0 1]);

subplot(1,3,3);
bar([ev_recur_density, re_recur_density, cross_recur_density],'FaceColor','flat','EdgeColor','none', ...
    'CData',[0.8 0.3 0.3; 0.3 0.5 0.8; 0.3 0.7 0.3]);
set(gca,'XTick',1:3,'XTickLabel',{'Evoked','Recovery',sprintf('Cross\n(Ev-in-Re)')}, 'FontSize',9);
xtickangle(0);
ylabel('Recurrence density'); title('Attractor comparison via recurrence density','FontWeight','normal'); ylim([0 1.05]);

sgtitle(sprintf('%s | Attractor epoch analysis', recording_ID),'FontWeight','normal');
exportgraphics(fig, fullfile(figuresDir, '05d_attractor_comparison.png'), 'Resolution', 550);

end
