%% VALIDATE_RECURRENCEDENSITY
%  Synthetic validation of fn_recurrenceDensity / fn_crossRecurrenceDensity
%  against two hand-built ground-truth scenarios:
%
%    Case 1 - four-phase simulation (baseline drift -> evoked limit cycle
%             -> transient -> recovery limit cycle), including a "1b"
%             variant where recovery settles onto a DIFFERENT attractor,
%             and a graded attractor-distance sweep.
%    Case 2 - a single stationary limit cycle throughout (no change; a
%             negative control), including an early-vs-late cross-check.
%
%  This is a unit test of the toolbox, not part of the analysis pipeline -
%  it uses no real recording data.

clear; close all; clc;

toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(toolboxRoot, 'Functions')));
run(fullfile(toolboxRoot, 'Config_UserSettings.m'));   % defines `cfg`
if ~exist(cfg.FIGURES_DIR,'dir'), mkdir(cfg.FIGURES_DIR); end

rng(cfg.validation_seed_recurrence);
fs_syn    = 50;        % Hz
f_osc_syn = 1/30;      % Hz (period = 30s, matches real data)
noise_amp = 0.05;

MAX_WIN_PTS  = cfg.MAX_WIN_PTS;
MAX_FULL_PTS = cfg.MAX_FULL_PTS;
MIN_LAG_S    = cfg.MIN_LAG_S;
ONSET_RATIO  = cfg.ONSET_RATIO;

fprintf('\n=== Synthetic Validation: Recurrence Density ===\n');

%% ---- Case 1: four-phase simulation with temporal changes ---------------
T_syn    = 600;
t_syn    = (0:T_syn*fs_syn-1)' / fs_syn;
t_P9_syn = 100; t_C2_syn = 400; t_ret_syn = 440;

Z_syn = zeros(T_syn*fs_syn, 2);
for k = 2:T_syn*fs_syn
    tk = t_syn(k);
    if tk < t_P9_syn
        Z_syn(k,:) = Z_syn(k-1,:) + 0.20*randn(1,2);
    elseif tk < t_C2_syn
        r = min(1, (tk-t_P9_syn)/20);
        Z_syn(k,:) = r*[cos(2*pi*f_osc_syn*tk), sin(2*pi*f_osc_syn*tk)] + noise_amp*randn(1,2);
    elseif tk < t_ret_syn
        Z_syn(k,:) = Z_syn(k-1,:) + 0.25*randn(1,2) + [3*exp(-(tk-t_C2_syn)/8), 0];
    else
        r = min(1, (tk-t_ret_syn)/20);
        Z_syn(k,:) = r*[cos(2*pi*f_osc_syn*tk), sin(2*pi*f_osc_syn*tk)] + noise_amp*randn(1,2);
    end
end

idx_ev_syn  = t_syn >= t_P9_syn & t_syn < t_C2_syn;
traj_ev_syn = Z_syn(idx_ev_syn,:); t_ev_syn = t_syn(idx_ev_syn);
ds_ev_syn   = max(1, floor(sum(idx_ev_syn)/MAX_WIN_PTS));
traj_ev_syn_ds = traj_ev_syn(1:ds_ev_syn:end,:); t_ev_syn_ds = t_ev_syn(1:ds_ev_syn:end);
eps_syn = fn_safeEpsilon(traj_ev_syn_ds);

idx_re_syn  = t_syn >= t_ret_syn;
traj_re_syn = Z_syn(idx_re_syn,:); t_re_syn = t_syn(idx_re_syn);
ds_re_syn   = max(1, floor(size(traj_re_syn,1)/MAX_WIN_PTS));
traj_re_syn_ds = traj_re_syn(1:ds_re_syn:end,:); t_re_syn_ds = t_re_syn(1:ds_re_syn:end);

% Full-recording sliding recurrence density (per-point, then binned)
ds_full_syn   = max(1, floor(size(Z_syn,1)/MAX_FULL_PTS));
idx_full_syn  = 1:ds_full_syn:size(Z_syn,1);
traj_full_syn = Z_syn(idx_full_syn,:); t_full_syn = t_syn(idx_full_syn);
buf_syn = 30; % one oscillation period as the "not enough future" buffer

[~, ~, recurs_full_syn, tested_full_syn] = fn_recurrenceDensity(traj_full_syn, t_full_syn, eps_syn, MIN_LAG_S, buf_syn);

win_f_syn = 30; step_f_syn = 5;
rr_syn_t = []; rr_syn_v = [];
for t0 = 0 : step_f_syn : T_syn-win_f_syn
    t1 = t0 + win_f_syn;
    inWin = t_full_syn >= t0 & t_full_syn < t1;
    tst   = inWin & tested_full_syn;
    rr_syn_t(end+1) = (t0+t1)/2; %#ok<*SAGROW>
    if any(tst), rr_syn_v(end+1) = mean(recurs_full_syn(tst)); else, rr_syn_v(end+1) = NaN; end
end

% Evoked / recovery self-recurrence density (own calibrated epsilon each),
% and evoked-in-recovery cross-recurrence density (shared, evoked eps).
eps_re_syn = fn_safeEpsilon(traj_re_syn_ds);
[ev_recur_density_syn, ~, ~, ~] = fn_recurrenceDensity(traj_ev_syn_ds, t_ev_syn_ds, eps_syn,   MIN_LAG_S, buf_syn);
[re_recur_density_syn, ~, ~, ~] = fn_recurrenceDensity(traj_re_syn_ds, t_re_syn_ds, eps_re_syn, MIN_LAG_S, buf_syn);
[cross_recur_density_syn, ~, ~] = fn_crossRecurrenceDensity(traj_ev_syn_ds, traj_re_syn_ds, eps_syn);
fprintf('  [Case 1: changing dynamics]  Evoked density (own eps)=%.3f  Recovery density (own eps)=%.3f  Cross density=%.3f\n', ...
    ev_recur_density_syn, re_recur_density_syn, cross_recur_density_syn);

[R_ev_syn, ~, ~] = fn_recurrencePlot(traj_ev_syn_ds, eps_syn, 'euc');
[R_re_syn, ~, ~] = fn_recurrencePlot(traj_re_syn_ds, eps_re_syn, 'euc');

n_ev_s = size(traj_ev_syn_ds,1); n_re_s = size(traj_re_syn_ds,1);
R_cross_syn = false(n_ev_s, n_re_s);
for ii = 1:n_ev_s
    d = sqrt(sum((traj_re_syn_ds - traj_ev_syn_ds(ii,:)).^2, 2));
    R_cross_syn(ii,:) = d <= eps_syn;
end

%% ---- Case 1b: recovery settles onto a DIFFERENT attractor ---------------
f_osc_diff = f_osc_syn * 2.3;
Z_diff = zeros(size(traj_re_syn));
for k = 1:size(Z_diff,1)
    tk = t_re_syn(k); r = min(1, (tk - t_re_syn(1))/20);
    Z_diff(k,:) = r*[3 + cos(2*pi*f_osc_diff*tk), 3 + sin(2*pi*f_osc_diff*tk)] + noise_amp*randn(1,2);
end
ds_diff = max(1, floor(size(Z_diff,1)/MAX_WIN_PTS));
traj_diff_ds = Z_diff(1:ds_diff:end,:); t_diff_ds = t_re_syn(1:ds_diff:end);

eps_diff_syn = fn_safeEpsilon(traj_diff_ds);
[re_diff_recur_density, ~, ~, ~] = fn_recurrenceDensity(traj_diff_ds, t_diff_ds, eps_diff_syn, MIN_LAG_S, buf_syn);
[cross_recur_density_diff, ~, ~] = fn_crossRecurrenceDensity(traj_ev_syn_ds, traj_diff_ds, eps_syn);
fprintf('  [Case 1b: DIFFERENT attractor in recovery]  Recovery self-density=%.3f (expect ~1) | Cross density=%.3f (expect ~0)\n', ...
    re_diff_recur_density, cross_recur_density_diff);

%% ---- Extension A: graded attractor-distance sweep (extends Case 1b) ----
fprintf('\nCase 1b extension: graded attractor-distance sweep...\n');
offsets_grid = [0 0.5 1 1.5 2 2.5 3 4];
cross_density_grid = zeros(size(offsets_grid));
for oi = 1:numel(offsets_grid)
    off = offsets_grid(oi);
    Z_off = zeros(size(traj_re_syn));
    for k = 1:size(Z_off,1)
        tk = t_re_syn(k); r = min(1, (tk - t_re_syn(1))/20);
        Z_off(k,:) = r*[off + cos(2*pi*f_osc_syn*tk), off + sin(2*pi*f_osc_syn*tk)] + noise_amp*randn(1,2);
    end
    ds_off = max(1, floor(size(Z_off,1)/MAX_WIN_PTS));
    traj_off_ds = Z_off(1:ds_off:end,:);
    [cross_density_grid(oi), ~, ~] = fn_crossRecurrenceDensity(traj_ev_syn_ds, traj_off_ds, eps_syn);
end
fprintf('  Cross-recurrence density vs. offset: %s\n', mat2str(round(cross_density_grid,2)));

%% ---- Case 2: no change - single stationary limit cycle throughout ------
T_syn2 = 600;
t_syn2 = (0:T_syn2*fs_syn-1)' / fs_syn;
Z_syn2 = [cos(2*pi*f_osc_syn*t_syn2), sin(2*pi*f_osc_syn*t_syn2)] + noise_amp*randn(numel(t_syn2),2);

ds_full_syn2  = max(1, floor(size(Z_syn2,1)/MAX_FULL_PTS));
idx_full_syn2 = 1:ds_full_syn2:size(Z_syn2,1);
traj_full_syn2 = Z_syn2(idx_full_syn2,:); t_full_syn2 = t_syn2(idx_full_syn2);
eps_syn2 = fn_safeEpsilon(traj_full_syn2);

[~, ~, recurs_full_syn2, tested_full_syn2] = fn_recurrenceDensity(traj_full_syn2, t_full_syn2, eps_syn2, MIN_LAG_S, buf_syn);

rr_syn2_t = []; rr_syn2_v = [];
for t0 = 0 : step_f_syn : T_syn2-win_f_syn
    t1 = t0 + win_f_syn;
    inWin = t_full_syn2 >= t0 & t_full_syn2 < t1;
    tst   = inWin & tested_full_syn2;
    rr_syn2_t(end+1) = (t0+t1)/2;
    if any(tst), rr_syn2_v(end+1) = mean(recurs_full_syn2(tst)); else, rr_syn2_v(end+1) = NaN; end
end
fprintf('  [Case 2: no change] mean recurrence density = %.3f (expect ~1, dropping only near the very end)\n', nanmean(rr_syn2_v));

% ---- Extension B: Case 2 cross-recurrence, early half vs. late half -----
fprintf('\nCase 2 extension: cross-recurrence (early vs. late half)...\n');
t_e1_start = 0;   t_e1_end = 250;
t_e2_start = 350; t_e2_end = T_syn2;

idx_e1_syn2 = t_syn2 >= t_e1_start & t_syn2 < t_e1_end;
idx_e2_syn2 = t_syn2 >= t_e2_start & t_syn2 <= t_e2_end;
traj_e1_full = Z_syn2(idx_e1_syn2,:); t_e1_full = t_syn2(idx_e1_syn2);
traj_e2_full = Z_syn2(idx_e2_syn2,:); t_e2_full = t_syn2(idx_e2_syn2);

ds_e1 = max(1, floor(size(traj_e1_full,1)/MAX_WIN_PTS));
traj_e1_ds = traj_e1_full(1:ds_e1:end,:); t_e1_ds = t_e1_full(1:ds_e1:end);
ds_e2 = max(1, floor(size(traj_e2_full,1)/MAX_WIN_PTS));
traj_e2_ds = traj_e2_full(1:ds_e2:end,:); t_e2_ds = t_e2_full(1:ds_e2:end);

eps_e1_syn2 = fn_safeEpsilon(traj_e1_ds);
eps_e2_syn2 = fn_safeEpsilon(traj_e2_ds);

[e1_recur_density_syn2, ~, ~, ~] = fn_recurrenceDensity(traj_e1_ds, t_e1_ds, eps_e1_syn2, MIN_LAG_S, buf_syn);
[e2_recur_density_syn2, ~, ~, ~] = fn_recurrenceDensity(traj_e2_ds, t_e2_ds, eps_e2_syn2, MIN_LAG_S, buf_syn);
[cross_density_syn2, ~, ~]       = fn_crossRecurrenceDensity(traj_e1_ds, traj_e2_ds, eps_e1_syn2);

[R_e1_syn2, ~, ~] = fn_recurrencePlot(traj_e1_ds, eps_e1_syn2, 'euc');
[R_e2_syn2, ~, ~] = fn_recurrencePlot(traj_e2_ds, eps_e2_syn2, 'euc');

n_e1_syn2 = size(traj_e1_ds,1); n_e2_syn2 = size(traj_e2_ds,1);
R_cross_syn2 = false(n_e1_syn2, n_e2_syn2);
for ii = 1:n_e1_syn2
    d = sqrt(sum((traj_e2_ds - traj_e1_ds(ii,:)).^2, 2));
    R_cross_syn2(ii,:) = d <= eps_e1_syn2;
end
fprintf('  [Case 2 cross-check] early self=%.3f  late self=%.3f  cross (early-in-late)=%.3f (all expected ~1)\n', ...
    e1_recur_density_syn2, e2_recur_density_syn2, cross_density_syn2);

%% =========================================================================
%  FIGURES
%  =========================================================================
col_evoked   = [0.00 0.45 0.70];
col_recovery = [0.85 0.35 0.05];
col_diff     = [0.00 0.60 0.50];
col_transit  = [0.65 0.65 0.65];
col_e1       = [0.25 0.45 0.65];
col_e2       = [0.75 0.55 0.20];

FONT = 'Arial'; FONT_SZ = 10; LETTER_SZ = 14; DPI = 500;

% --- FIGURE 1: Case 1 (changing dynamics) ---
FIG1_W = 1600; FIG1_H = 950;
MARGIN = 0.055; GAPX = 0.035; GAPY = 0.10;
rowH = [0.34 0.34 0.22];

fig1 = figure('Color','w','Units','pixels','Position',[50 50 FIG1_W FIG1_H]);

axA = axes('Parent',fig1,'Position', fn_axpos(1,1,1,2, rowH, MARGIN, GAPX, GAPY));
plot(axA, t_syn, Z_syn(:,1), 'Color', [0.15 0.15 0.15], 'LineWidth', 1); hold(axA,'on');
patch(axA, [t_P9_syn t_C2_syn t_C2_syn t_P9_syn], [-5 -5 6 6], col_evoked, 'EdgeColor','none','FaceAlpha',0.10);
patch(axA, [t_C2_syn t_ret_syn t_ret_syn t_C2_syn], [-5 -5 6 6], col_transit,'EdgeColor','none','FaceAlpha',0.10);
patch(axA, [t_ret_syn T_syn T_syn t_ret_syn], [-5 -5 6 6], col_recovery,'EdgeColor','none','FaceAlpha',0.10);
xlabel(axA,'Time (s)'); ylabel(axA,'Activity (a.u.)');
title(axA,'Simulated state trajectory','FontWeight','normal');
xlim(axA,[0 T_syn]); ylim(axA,[-5 6]); box(axA,'off');

axB = axes('Parent',fig1,'Position', fn_axpos(1,3,1,2, rowH, MARGIN, GAPX, GAPY));
plot(axB, rr_syn_t, rr_syn_v, 'Color',[0.1 0.1 0.1], 'LineWidth', 1.5); hold(axB,'on');
yline(axB, ONSET_RATIO, '--', sprintf('%.0f%% criterion', ONSET_RATIO*100), ...
      'Color', [0.5 0.5 0.5], 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
xline(axB, t_P9_syn, 'Color', col_evoked, 'LineWidth', 1.5);
xline(axB, t_C2_syn, 'Color', col_transit, 'LineWidth', 1.5);
xline(axB, t_ret_syn,'Color', col_recovery,'LineWidth', 1.5, 'LineStyle', ':');
xlabel(axB,'Time (s)'); ylabel(axB,'Recurrence density');
title(axB,'Sliding recurrence density','FontWeight','normal');
xlim(axB,[0 T_syn]); ylim(axB,[0 1.05]); box(axB,'off');

axC = axes('Parent',fig1,'Position', fn_axpos(2,1,1,1, rowH, MARGIN, GAPX, GAPY));
imagesc(axC, R_ev_syn); colormap(axC,[1 1 1; col_evoked]); axis(axC,'xy','square');
set(axC,'XTick',[],'YTick',[]); box(axC,'on');
title(axC, sprintf('Evoked (self)\ndensity = %.2f', ev_recur_density_syn), 'FontWeight','normal');

axD = axes('Parent',fig1,'Position', fn_axpos(2,2,1,1, rowH, MARGIN, GAPX, GAPY));
imagesc(axD, R_re_syn); colormap(axD,[1 1 1; col_recovery]); axis(axD,'xy','square');
set(axD,'XTick',[],'YTick',[]); box(axD,'on');
title(axD, sprintf('Recovery (self)\ndensity = %.2f', re_recur_density_syn), 'FontWeight','normal');

axE = axes('Parent',fig1,'Position', fn_axpos(2,3,1,1, rowH, MARGIN, GAPX, GAPY));
imagesc(axE, R_cross_syn'); colormap(axE,[1 1 1; 0.15 0.15 0.15]); axis(axE,'xy','square');
xlabel(axE,'Evoked'); ylabel(axE,'Recovery'); set(axE,'XTick',[],'YTick',[]); box(axE,'on');
title(axE, sprintf('Cross-recurrence\ndensity = %.2f', cross_recur_density_syn), 'FontWeight','normal');

axF = axes('Parent',fig1,'Position', fn_axpos(2,4,1,1, rowH, MARGIN, GAPX, GAPY));
bar_vals = [ev_recur_density_syn, re_recur_density_syn, re_diff_recur_density, cross_recur_density_syn, cross_recur_density_diff];
bar_colors = [col_evoked; col_recovery; col_diff; 0.35 0.35 0.35; 0.60 0.60 0.60];
bF = bar(axF, 1:5, bar_vals, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth', 0.8);
bF.CData = bar_colors;
set(axF,'XTick',1:5,'XTickLabel',{'Evoked','Recovery','Diff.','Cross(Ev-Rec)','Cross(Ev-Diff)'},'XTickLabelRotation',45);
ylabel(axF,'Density'); ylim(axF,[0 1.05]); box(axF,'off');
title(axF,'Summary','FontWeight','normal');

axG = axes('Parent',fig1,'Position', fn_axpos(3,2,1,2, rowH, MARGIN, GAPX, GAPY));
plot(axG, offsets_grid, cross_density_grid, '-o', 'Color',[0.25 0.30 0.35], ...
     'MarkerFaceColor',[0.25 0.30 0.35], 'LineWidth',1.5, 'MarkerSize',5);
xlabel(axG,'Spatial offset'); ylabel(axG,'Cross-recurrence density');
title(axG,'Attractor-distance sweep (sensitivity control)','FontWeight','normal');
ylim(axG,[0 1.05]); box(axG,'off');

fn_letterAxes([axA axB axC axD axE axF axG], LETTER_SZ, FONT);
fn_formatAxes(fig1, FONT, FONT_SZ);
drawnow;
exportgraphics(fig1, fullfile(cfg.FIGURES_DIR, 'validation_recurrence_case1_changing_dynamics.png'), 'Resolution', DPI);

% --- FIGURE 2: Case 2 (stationary control) ---
FIG2_W = 1600; FIG2_H = 720;
fig2 = figure('Color','w','Units','pixels','Position',[50 50 FIG2_W FIG2_H]);
tlo = tiledlayout(fig2, 2, 4, 'TileSpacing', 'normal', 'Padding', 'normal');

axH = nexttile(tlo, 1, [1 2]);
plot(axH, t_syn2, Z_syn2(:,1), 'Color', [0.15 0.15 0.15], 'LineWidth', 1); hold(axH,'on');
patch(axH, [t_e1_start t_e1_end t_e1_end t_e1_start], [-2 -2 2 2], col_e1, 'EdgeColor','none','FaceAlpha',0.10);
patch(axH, [t_e2_start t_e2_end t_e2_end t_e2_start], [-2 -2 2 2], col_e2, 'EdgeColor','none','FaceAlpha',0.10);
xlabel(axH,'Time (s)'); ylabel(axH,'Activity (a.u.)');
title(axH,'Stationary control trajectory','FontWeight','normal');
xlim(axH,[0 T_syn2]); ylim(axH,[-2 2]); box(axH,'off');

axI = nexttile(tlo, 3, [1 2]);
plot(axI, rr_syn2_t, rr_syn2_v, 'Color',[0.1 0.1 0.1], 'LineWidth', 1.5); hold(axI,'on');
yline(axI, ONSET_RATIO, '--', 'Color',[0.5 0.5 0.5], 'LineWidth',1.2);
xline(axI, t_e1_start, 'Color', col_e1, 'LineStyle',':'); xline(axI, t_e1_end, 'Color', col_e1, 'LineStyle',':');
xline(axI, t_e2_start, 'Color', col_e2, 'LineStyle',':'); xline(axI, t_e2_end,  'Color', col_e2, 'LineStyle',':');
xlabel(axI,'Time (s)'); ylabel(axI,'Recurrence density');
title(axI,'Sliding recurrence density','FontWeight','normal');
xlim(axI,[0 T_syn2]); ylim(axI,[0 1.05]); box(axI,'off');

axJ = nexttile(tlo, 5, [1 1]);
imagesc(axJ, R_e1_syn2); colormap(axJ,[1 1 1; col_e1]); axis(axJ,'xy','square');
set(axJ,'XTick',[],'YTick',[]); box(axJ,'on');
title(axJ, sprintf('Early (self)\ndensity = %.2f', e1_recur_density_syn2), 'FontWeight','normal');

axK = nexttile(tlo, 6, [1 1]);
imagesc(axK, R_e2_syn2); colormap(axK,[1 1 1; col_e2]); axis(axK,'xy','square');
set(axK,'XTick',[],'YTick',[]); box(axK,'on');
title(axK, sprintf('Late (self)\ndensity = %.2f', e2_recur_density_syn2), 'FontWeight','normal');

axL = nexttile(tlo, 7, [1 1]);
imagesc(axL, R_cross_syn2'); colormap(axL,[1 1 1; 0.15 0.15 0.15]); axis(axL,'xy','square');
xlabel(axL,'Early'); ylabel(axL,'Late'); set(axL,'XTick',[],'YTick',[]); box(axL,'on');
title(axL, sprintf('Cross-recurrence\ndensity = %.2f', cross_density_syn2), 'FontWeight','normal');

axM = nexttile(tlo, 8, [1 1]);
bar_vals2 = [e1_recur_density_syn2, e2_recur_density_syn2, cross_density_syn2];
bar_colors2 = [col_e1; col_e2; 0.35 0.35 0.35];
bM = bar(axM, 1:3, bar_vals2, 'FaceColor','flat', 'EdgeColor','k', 'LineWidth', 0.8, 'BarWidth', 0.6);
bM.CData = bar_colors2;
set(axM,'XTick',1:3,'XTickLabel',{'Early','Late','Cross'});
ylabel(axM,'Density'); ylim(axM,[0 1.05]); box(axM,'off');
title(axM,'Summary','FontWeight','normal');

drawnow;
fn_formatAxes(fig2, FONT, FONT_SZ);
legend(axH, 'off');
fn_letterAxes([axH axI axJ axK axL axM], LETTER_SZ, FONT);
drawnow;
exportgraphics(fig2, fullfile(cfg.FIGURES_DIR, 'validation_recurrence_case2_stationary_control.png'), 'Resolution', DPI);

fprintf('\n=== Synthetic Validation (Recurrence Density) complete. Figures saved:\n');
fprintf('    validation_recurrence_case1_changing_dynamics.png\n');
fprintf('    validation_recurrence_case2_stationary_control.png\n');
