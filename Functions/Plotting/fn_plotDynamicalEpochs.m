function fn_plotDynamicalEpochs(dyn, recording_ID, figuresDir)
% FN_PLOTDYNAMICALEPOCHS  Summary figure for fn_dynamicalEpochs.
%   (a) cycle x cycle subspace alignment   (b) cycle clusters on the time axis
%   (c) cycles per cluster per protocol epoch   (d) mean alignment between epochs
%
%   FN_PLOTDYNAMICALEPOCHS(dyn, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok
    fprintf('  Dynamical epochs: nothing to plot (%s)\n', dyn.message);
    return
end

nSeg = numel(dyn.segment_s); nG = dyn.nGroups; nE = numel(dyn.ep_labels);
cmap = lines(max(nG, 1));
bnd  = dyn.ep_bounds(2:end, 1);     % epoch boundaries

fig = figure('Name', 'Dynamical epochs', 'Position', [100 100 1400 800], 'Visible', 'off');

% (a) cycle-by-cycle alignment
subplot(2,2,1);
imagesc(dyn.similarity); clim([0 1]); colormap(gca, hot); colorbar; axis square; hold on;
for k = 1:numel(bnd)
    idx = find(dyn.cycle_mid >= bnd(k), 1, 'first');
    if ~isempty(idx) && idx > 1
        xline(idx - 0.5, 'w-', 'LineWidth', 1); yline(idx - 0.5, 'w-', 'LineWidth', 1);
    end
end
xlabel('Cycle #'); ylabel('Cycle #'); title('Cycle-by-cycle subspace alignment');

% (b) clusters on the time axis
subplot(2,2,2); hold on;
for c = 1:nSeg
    t1 = dyn.cycle_t(c,1); t2 = dyn.cycle_t(c,2);
    patch([t1 t2 t2 t1], [0 0 1 1], cmap(dyn.grps(c),:), 'EdgeColor', 'w');
end
xl = dyn.t_range; xlim(xl); ylim([0 1.15]);
for k = 1:numel(bnd)
    if bnd(k) > xl(1) && bnd(k) < xl(2), xline(bnd(k), 'k--', 'LineWidth', 1); end
end
for e = 1:nE
    a = max(dyn.ep_bounds(e,1), xl(1)); b = min(dyn.ep_bounds(e,2), xl(2));
    if b > a
        text((a+b)/2, 1.07, dyn.ep_labels{e}, 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
hL = gobjects(1, nG);
for g = 1:nG, hL(g) = patch(NaN, NaN, cmap(g,:)); end
legend(hL, arrayfun(@(g) sprintf('Cluster %d', g), 1:nG, 'UniformOutput', false), 'Location', 'eastoutside');
set(gca, 'YTick', []); xlabel('Time (s)'); title('Dynamical epochs (cycle clusters)'); box on;

% (c) cluster x protocol epoch counts
subplot(2,2,3);
imagesc(dyn.cluster_by_epoch); colormap(gca, flipud(gray)); axis tight;
for g = 1:nG
    for e = 1:nE
        text(e, g, sprintf('%d', dyn.cluster_by_epoch(g,e)), 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
    end
end
set(gca, 'XTick', 1:nE, 'XTickLabel', dyn.ep_labels, 'YTick', 1:nG, 'XTickLabelRotation', 30);
xlabel('Protocol epoch'); ylabel('Cluster'); title('Cycles per cluster in each protocol epoch');

% (d) alignment between protocol epochs
subplot(2,2,4);
imagesc(dyn.epoch_sim); clim([0 1]); colormap(gca, hot); colorbar; axis square;
for a = 1:nE
    for b = 1:nE
        if ~isnan(dyn.epoch_sim(a,b))
            text(b, a, sprintf('%.2f', dyn.epoch_sim(a,b)), 'HorizontalAlignment', 'center', ...
                'Color', 'w', 'FontWeight', 'bold');
        end
    end
end
set(gca, 'XTick', 1:nE, 'XTickLabel', dyn.ep_labels, 'YTick', 1:nE, 'YTickLabel', dyn.ep_labels, ...
    'XTickLabelRotation', 30);
title('Mean alignment between epochs (diagonal = within)');

sgtitle(sprintf('Dynamical epochs (jPCA cycles) | %s', recording_ID), 'Interpreter', 'none');
exportgraphics(fig, fullfile(figuresDir, '06_dynamical_epochs.png'), 'Resolution', 500);
close(fig);
end
