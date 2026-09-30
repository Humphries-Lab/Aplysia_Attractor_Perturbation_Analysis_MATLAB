function fn_plotDynamicalCycleMetrics(cyc, dyn, recording_ID, figuresDir) %#ok<INUSL>
% FN_PLOTDYNAMICALCYCLEMETRICS  Per-cycle dimensionality, attractor movement, stimulus phase.
%
% Suggested caption:
%   "Cycles are recurrence-defined rotations of the jPCA plane; each point is one cycle,
%   plotted at its midpoint. Dashed lines mark the stimuli (P9, C2). (A) Participation ratio
%   (PR/N) of the population covariance within each cycle. (B) Subspace alignment index
%   (Elsayed et al. 2016; fraction of one cycle's variance captured by another cycle's
%   subspace, symmetrised) between each cycle and the previous cycle (blue) or the first
%   Evoked cycle (grey; triangle). (C) Fraction of the cycle elapsed when each stimulus was
%   applied (0 = cycle start, 1 = cycle end). Cycles are cut where the trajectory recurs, so
%   phase 0 is not a fixed point on the attractor."

if ~isfield(cyc, 'ok') || ~cyc.ok, return; end
colA = [0.8 0.1 0.1]; colB = [0.1 0.3 0.7]; colG = [0.5 0.5 0.5];
t = cyc.cycle_mid; nS = numel(cyc.stim);
xl = [cyc.cycle_t(1,1) cyc.cycle_t(end,2)];

fig = figure('Name', 'Cycle metrics', 'Position', [100 100 800 850], 'Visible', 'off', 'Color', 'w');
tl = tiledlayout(fig, 5, 1, 'TileSpacing', 'compact', 'Padding', 'loose');

% A - PR per cycle
axA = nexttile(tl, [2 1]); hold(axA, 'on');
plot(axA, t, cyc.PR_norm, '-o', 'Color', colA, 'MarkerFaceColor', colA, 'MarkerSize', 4, 'LineWidth', 1.25);
ylabel(axA, 'PR / N'); xlim(axA, xl); set(axA, 'XTickLabel', []);
title(axA, 'Dimensionality of each cycle', 'FontWeight', 'normal', 'FontSize', 10, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02 0]);
fn_stimLines(axA, cyc.stim);

% B - alignment per cycle
axB = nexttile(tl, [2 1]); hold(axB, 'on');
hR = plot(axB, t, cyc.align_ref,  '-s', 'Color', colG, 'MarkerFaceColor', colG, 'MarkerSize', 4, 'LineWidth', 1.25);
hP = plot(axB, t, cyc.align_prev, '-o', 'Color', colB, 'MarkerFaceColor', colB, 'MarkerSize', 4, 'LineWidth', 1.25);
hT = plot(axB, t(cyc.ref_cycle), 0.03, '^', 'MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colG, 'LineWidth', 1.25);
ylabel(axB, 'Subspace alignment'); ylim(axB, [0 1]); xlim(axB, xl);
xlabel(axB, 'Time (s)');
title(axB, 'Alignment of each cycle with other cycles', 'FontWeight', 'normal', 'FontSize', 10, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02 0]);
fn_stimLines(axB, cyc.stim);
lg = legend(axB, [hP hR hT], {'With previous cycle', 'With reference cycle', 'Reference cycle (first Evoked)'}, ...
    'Location', 'northeast', 'Box', 'off', 'Interpreter', 'none');
lg.FontSize = 9;

% C - phase of each stimulus in its cycle
axC = nexttile(tl); hold(axC, 'on');
for k = 1:nS
    s = cyc.stim(k);
    plot(axC, [0 1], [k k], '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 6);
    if ~isnan(s.phase_frac)
        plot(axC, s.phase_frac, k, 'v', 'MarkerSize', 9, 'MarkerFaceColor', colA, 'MarkerEdgeColor', 'k');
        text(axC, s.phase_frac, k - 0.38, sprintf('cycle %d', s.cycle), 'HorizontalAlignment', 'center', ...
            'FontSize', 9, 'Interpreter', 'none');
    end
end
set(axC, 'YTick', 1:nS, 'YTickLabel', {cyc.stim.name}, 'YDir', 'reverse', 'XTick', 0:0.25:1);
xlim(axC, [0 1]); ylim(axC, [0.2 nS+0.5]);
xlabel(axC, 'Phase within cycle (fraction elapsed)');
title(axC, 'Stimulus phase', 'FontWeight', 'normal', 'FontSize', 10, 'HorizontalAlignment', 'left', 'Units', 'normalized', 'Position', [0 1.02 0]);

fn_formatAxes(fig, 'Arial', 10);

% Panel letters placed in figure coordinates so they are never clipped or overlap labels
drawnow;
letters = {'A','B','C'}; axs = [axA axB axC];
for i = 1:3
    p = axs(i).OuterPosition;
    annotation(fig, 'textbox', [0.005, p(2)+p(4)-0.045, 0.05, 0.04], 'String', letters{i}, ...
        'EdgeColor', 'none', 'FontName', 'Arial', 'FontSize', 13, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top');
end

exportgraphics(fig, fullfile(figuresDir, '06d_dynamical_cycle_metrics.png'), 'Resolution', 500);
close(fig);
end

function fn_stimLines(ax, stim)
for k = 1:numel(stim)
    xline(ax, stim(k).t, 'k--', stim(k).name, 'LineWidth', 1, 'LabelOrientation', 'horizontal', ...
        'LabelVerticalAlignment', 'top', 'LabelHorizontalAlignment', 'left', 'FontSize', 9, 'Interpreter', 'none');
end
end