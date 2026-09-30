function fn_plotDynamicalControls(dyn, surr, recording_ID, figuresDir) %#ok<INUSD>
% FN_PLOTDYNAMICALCONTROLS  Controls for the dynamical-epoch analysis.
%   A  Segmentation contrast (within - between block alignment) for every possible split of the
%      cycle sequence; dashed lines = protocol landmarks; the maximum is the strongest unsupervised split.
%   B  Permutation null for the mean distance between protocol landmarks and the nearest change point.
%   C  Real data vs surrogates (neurons circularly shifted independently).
%
%   FN_PLOTDYNAMICALCONTROLS(dyn, surr, recording_ID, figuresDir)

if ~isfield(dyn, 'ok') || ~dyn.ok, return; end
blue = [0.1 0.3 0.7]; gry = [0.5 0.5 0.5]; red = [0.8 0.1 0.1];

fig = figure('Name', 'Dynamical epochs: controls', 'Position', [100 100 1500 430], 'Visible', 'off');

%% A - segmentation contrast
ax1 = subplot(1,3,1); hold on;
bd = dyn.boundary;
if bd.ok
    plot(bd.t, bd.contrast, '-', 'Color', blue, 'LineWidth', 1.5);
    plot(bd.best_t, bd.best_contrast, 'o', 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 7);
    yl = ylim; yl(2) = yl(2) + 0.12 * diff(yl); ylim(yl);
    lm = dyn.landmarks;
    for k = 1:numel(lm.t)
        if lm.t(k) > dyn.t_range(1) && lm.t(k) < dyn.t_range(2)
            xline(lm.t(k), '--', 'Color', red, 'LineWidth', 1.0);
            text(lm.t(k), yl(2) - mod(k, 2) * 0.07 * diff(yl), lm.short{k}, 'Color', red, ...
                'FontSize', 8, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'BackgroundColor', 'w', 'Margin', 0.5);
        end
    end
    xlim(dyn.t_range);
else
    text(0.5, 0.5, 'Too few cycles', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Split time (s)'); ylabel('Alignment contrast (within - between)');

%% B - landmark permutation
ax2 = subplot(1,3,2); hold on;
lt = dyn.landmark_test;
if lt.ok
    histogram(lt.null, 30, 'FaceColor', gry, 'EdgeColor', 'w', 'Normalization', 'probability');
    xline(lt.observed, '-', 'Color', red, 'LineWidth', 2);
    yl = ylim; ylim([0 yl(2) * 1.2]);
    text(0.97, 0.93, sprintf('p = %.3f', lt.p), 'Units', 'normalized', 'HorizontalAlignment', 'right', 'FontSize', 9);
else
    text(0.5, 0.5, 'No change points in window', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Mean distance, landmark to nearest change (s)'); ylabel('Proportion of permutations');

%% C - real vs surrogate
ax3 = subplot(1,3,3); hold on;
if isstruct(surr) && isfield(surr, 'real')
    names = {'Rotation (6 PCs)', 'Within - between'};
    sv = [surr.rot6(:), surr.contrast(:)];
    rv = [surr.real.rot6, surr.real.contrast];
    rng_state = rng; rng(1);
    for m = 1:2
        plot(m + 0.15 * (rand(size(sv,1),1) - 0.5), sv(:,m), 'o', 'MarkerFaceColor', gry, ...
            'MarkerEdgeColor', 'w', 'MarkerSize', 7);
    end
    rng(rng_state);
    hS = plot(NaN, NaN, 'o', 'MarkerFaceColor', gry, 'MarkerEdgeColor', 'w', 'MarkerSize', 7);
    hR = plot(1:2, rv, 'd', 'MarkerFaceColor', blue, 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'LineStyle', 'none');
    set(gca, 'XTick', 1:2, 'XTickLabel', names); xlim([0.5 2.5]);
    legend([hR hS], {'Data', sprintf('Surrogates (n=%d)', surr.n)}, 'Location', 'northeast', 'Box', 'off', 'FontSize', 9);
    ylabel('Value');
else
    text(0.5, 0.5, 'Surrogates not run', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off
end

fn_formatAxes(fig, 'Arial', 10);
fn_letterAxes([ax1 ax2 ax3], 14, 'Arial');
exportgraphics(fig, fullfile(figuresDir, '06c_dynamical_controls.png'), 'Resolution', 500);
close(fig);
end