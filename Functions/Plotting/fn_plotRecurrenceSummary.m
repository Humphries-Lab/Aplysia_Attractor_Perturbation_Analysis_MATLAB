function fn_plotRecurrenceSummary(traj_ev, t_ev_ax, traj_re, t_re_ax, R_ev, R_re, ...
    ev_recur_density, re_recur_density, cross_recur_density, recurs_ev_in_re, ...
    epsilon_rr, recording_ID, figuresDir)
% FN_PLOTRECURRENCESUMMARY  Four-panel recurrence-analysis figure: the
% Evoked and Recovery visual recurrence plots, the evoked-in-recovery
% cross-recurrence matrix, and a stem plot of which evoked moments recur
% in recovery.
%
%   FN_PLOTRECURRENCESUMMARY(traj_ev, t_ev_ax, traj_re, t_re_ax, R_ev, R_re, ...
%       ev_recur_density, re_recur_density, cross_recur_density, recurs_ev_in_re, ...
%       epsilon_rr, recording_ID, figuresDir)

fig = figure('Name','Recurrence analysis','Position',[50 50 1500 520]);

subplot(1,4,1);
imagesc(t_ev_ax, t_ev_ax, R_ev);
colormap(gca,[1 1 1; 0 0 0]); axis xy square;
xlabel('Time (s)'); ylabel('Time (s)');
title(sprintf('Evoked (visual RP)\nRecurrence density=%.3f', ev_recur_density),'FontWeight','normal');

subplot(1,4,2);
imagesc(t_re_ax, t_re_ax, R_re);
colormap(gca,[1 1 1; 0 0 0]); axis xy square;
xlabel('Time (s)'); ylabel('Time (s)');
title(sprintf('Recovery (visual RP)\nRecurrence density=%.3f', re_recur_density),'FontWeight','normal');

subplot(1,4,3);
n_ev = size(traj_ev,1); n_re = size(traj_re,1);
BLK = 200; R_cross = false(n_ev, n_re);
for i0 = 1:BLK:n_ev
    i1 = min(i0+BLK-1, n_ev); A = traj_ev(i0:i1,:);
    for j0 = 1:BLK:n_re
        j1 = min(j0+BLK-1, n_re); B = traj_re(j0:j1,:);
        AA = sum(A.^2,2); BB = sum(B.^2,2);
        Dij = sqrt(max(bsxfun(@plus,AA,BB') - 2*(A*B'), 0));
        R_cross(i0:i1, j0:j1) = Dij <= epsilon_rr;
    end
end
imagesc(t_ev_ax, t_re_ax, R_cross');
colormap(gca,[1 1 1; 0 0 0]); axis xy;
xlabel('Evoked time (s)'); ylabel('Recovery time (s)');
title(sprintf('Cross-recurrence (visual)\nEv-in-Re density=%.3f', cross_recur_density),'FontWeight','normal');

subplot(1,4,4);
stem(t_ev_ax, double(recurs_ev_in_re), 'k.', 'MarkerSize', 4);
xlabel('Evoked time (s)'); ylabel('Recurs in recovery? (0/1)');
title('Which evoked moments recur in recovery','FontWeight','normal');
ylim([-0.1 1.1]); box on;

sgtitle(sprintf('%s | Recurrence density analysis  eps=%.4f', recording_ID, epsilon_rr), 'FontWeight','normal');
exportgraphics(fig, fullfile(figuresDir, '05_recurrence_all.png'), 'Resolution', 500);

end
