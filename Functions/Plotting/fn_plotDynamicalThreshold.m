function fn_plotDynamicalThreshold(dyn, recording_ID, figuresDir) %#ok<INUSD>
% FN_PLOTDYNAMICALTHRESHOLD  Change points and similarity threshold (publication layout).
%   A  Distribution of pairwise cycle alignment, same- vs different-cluster pairs,
%      with the shuffle-chance level (95th percentile) and the separating threshold.
%   B  Dynamical epochs (runs of one cluster) along time; black lines = change points;
%      grey strip = protocol epochs (used only to label, not to detect).
%   C  Alignment vs time between cycles, with a linear fit to same-epoch pairs.
%
%   FN_PLOTDYNAMICALTHRESHOLD(dyn, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok
    return
end
blue = [0.1 0.3 0.7]; gry = [0.5 0.5 0.5]; red = [0.8 0.1 0.1];
st = dyn.sim_thresh; pr = dyn.pairs;
nG = dyn.nGroups; cmap = lines(max(nG, 1));

fig = figure('Name', 'Dynamical epochs: threshold', 'Position', [100 100 1500 430], 'Visible', 'off');

% A - pairwise alignment
ax1 = subplot(1,3,1); hold on;
edges = 0:0.05:1;
w = pr.sim(pr.same_cluster); b = pr.sim(~pr.same_cluster);
hB = gobjects(0); hW = gobjects(0); top = 0;
if ~isempty(b)
    hB = histogram(b, edges, 'FaceColor', gry, 'EdgeColor', 'w', 'Normalization', 'probability');
    top = max(top, max(hB.Values));
end
if ~isempty(w)
    hW = histogram(w, edges, 'FaceColor', blue, 'EdgeColor', 'w', 'Normalization', 'probability', 'FaceAlpha', 0.6);
    top = max(top, max(hW.Values));
end
xlim([0 1]); ylim([0 max(top, eps) * 1.35]);
if ~isnan(st.chance95)
    xline(st.chance95, ':', 'Chance', 'Color', 'k', 'LineWidth', 1.5, 'FontSize', 9, ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
end
if ~isnan(st.separation)
    xline(st.separation, '--', 'Threshold', 'Color', red, 'LineWidth', 1.5, 'FontSize', 9, ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left');
end
hs = gobjects(0); ls = {};
if ~isempty(hW), hs(end+1) = hW; ls{end+1} = 'Same cluster'; end
if ~isempty(hB), hs(end+1) = hB; ls{end+1} = 'Different cluster'; end
if ~isempty(hs), legend(hs, ls, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9); end
if ~isnan(st.auc)
    text(0.97, 0.70, sprintf('AUC = %.2f', st.auc), 'Units', 'normalized', ...
        'HorizontalAlignment', 'right', 'FontSize', 9);
end
xlabel('Subspace alignment'); ylabel('Proportion of cycle pairs');

% B - dynamical epochs along time
ax2 = subplot(1,3,2); hold on;
for k = 1:numel(dyn.runs.cluster)
    t1 = dyn.runs.t(k,1); t2 = dyn.runs.t(k,2);
    patch([t1 t2 t2 t1], [0 0 1 1], cmap(dyn.runs.cluster(k),:), 'EdgeColor', 'none');
end
for k = 1:numel(dyn.change_t)
    xline(dyn.change_t(k), '-', 'Color', 'k', 'LineWidth', 0.75);
end
xl = dyn.t_range;
for e = 1:numel(dyn.ep_labels)                       % protocol-epoch strip (labels only)
    a = max(dyn.ep_bounds(e,1), xl(1)); bb = min(dyn.ep_bounds(e,2), xl(2));
    if bb > a
        shade = 0.93 - 0.08 * mod(e, 2);
        patch([a bb bb a], [1.06 1.06 1.22 1.22], shade * [1 1 1], 'EdgeColor', 'w');
        lab = dyn.ep_labels{e};
        if (bb - a) / (xl(2) - xl(1)) < 0.12, lab = [lab(1:min(4, numel(lab))) '.']; end
        text((a + bb) / 2, 1.14, lab, 'HorizontalAlignment', 'center', 'FontSize', 8);
    end
end
xlim(xl); ylim([0 1.25]);
hL = gobjects(1, nG);
for g = 1:nG, hL(g) = patch(NaN, NaN, cmap(g,:), 'EdgeColor', 'none'); end
lgd = legend(hL, arrayfun(@(g) sprintf('%d', g), 1:nG, 'UniformOutput', false), ...
    'Location', 'southoutside', 'Orientation', 'horizontal', 'Box', 'off', 'FontSize', 9);
title(lgd, 'Cluster');
set(gca, 'YTick', []); xlabel('Time (s)'); box on;

%C - alignment vs time between cycles
ax3 = subplot(1,3,3); hold on;
m0 = pr.type == 0; m1 = pr.type == 1; m2 = pr.type == 2;
hp = gobjects(0); lp = {};
if any(m0), hp(end+1) = plot(pr.lag_s(m0), pr.sim(m0), '.', 'Color', [0.8 0.8 0.8], 'MarkerSize', 10); lp{end+1} = 'Other pairs'; end
if any(m1), hp(end+1) = plot(pr.lag_s(m1), pr.sim(m1), '.', 'Color', blue, 'MarkerSize', 12); lp{end+1} = 'Same epoch'; end
if any(m2), hp(end+1) = plot(pr.lag_s(m2), pr.sim(m2), '.', 'Color', red, 'MarkerSize', 12); lp{end+1} = 'Evoked-Recovery'; end
hasFit = isfield(dyn.drift, 'ok') && dyn.drift.ok;
if hasFit
    xx = [0 max(pr.lag_s)];
    hp(end+1) = plot(xx, dyn.drift.intercept + dyn.drift.slope_per_min / 60 * xx, 'k--', 'LineWidth', 1.2);
    lp{end+1} = 'Linear fit (same epoch)';
end
if ~isnan(st.chance95)
    hp(end+1) = yline(st.chance95, ':', 'Color', 'k', 'LineWidth', 1.2); lp{end+1} = 'Chance';
end
if ~isempty(hp), legend(hp, lp, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9); end
xlabel('Time between cycles (s)'); ylabel('Subspace alignment'); ylim([0 1]);
if hasFit
    title(sprintf('Evoked-Recovery: %.2f observed, %.2f expected', dyn.drift.obs_ER, dyn.drift.pred_ER), ...
        'FontWeight', 'normal', 'FontSize', 10);
end

fn_formatAxes(fig, 'Arial', 10);
fn_letterAxes([ax1 ax2 ax3], 14, 'Arial');
exportgraphics(fig, fullfile(figuresDir, '06b_dynamical_threshold.png'), 'Resolution', 500);
close(fig);
end