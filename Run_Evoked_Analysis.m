%% RUN_EVOKED_ANALYSIS
%  Attractor analysis for a single P9
%  stimulus-evoked motor program with no perturbation (Baseline + Evoked
%  only, no Recovery epoch).
%
%  Companion to Run_Attractor_Analysis.m. Uses the same Functions/ folder
%  and Config_UserSettings.m, unmodified. The shared plotting functions
%  all expect a C2 landmark and/or a Recovery epoch, so each one is
%  wrapped in a local "_evoked" adapter (bottom of this file) that calls
%  the original function unchanged, edits the resulting figure (removes
%  the C2 marker, or the panels that have no Recovery data), and
%  re-exports over the same file. Functions/ itself is never edited.

clear; close all; clc;

thisFile = mfilename('fullpath');
if contains(thisFile, fullfile(tempdir))
    toolboxRoot = pwd;
else
    toolboxRoot = fileparts(thisFile);
end
addpath(genpath(fullfile(toolboxRoot, 'Functions')));

cfgFile = fullfile(toolboxRoot, 'Config_UserSettings.m');
if ~isfile(cfgFile)
    error('Run_Evoked_Analysis:ConfigNotFound', ...
        'Config_UserSettings.m not found at:\n  %s', cfgFile);
end
run(cfgFile);   % defines `cfg`

if ~exist(cfg.RESULTS_DIR,'dir'), mkdir(cfg.RESULTS_DIR); end
if ~exist(cfg.FIGURES_DIR,'dir'), mkdir(cfg.FIGURES_DIR); end
% Append the recording_ID to the base figures directory and create the subfolder
cfg.FIGURES_DIR = fullfile(cfg.FIGURES_DIR, cfg.recording_ID);
if ~exist(cfg.FIGURES_DIR,'dir'), mkdir(cfg.FIGURES_DIR); end

if ~exist(cfg.LOGS_DIR,'dir'), mkdir(cfg.LOGS_DIR); end
diary off;   % close any diary left open by a previous run that errored mid-way
logFile = fullfile(cfg.LOGS_DIR, sprintf('%s_log_%s.txt', cfg.recording_ID, datestr(now,'yyyymmdd_HHMMSS')));
diary(logFile);
diary on;
fprintf('Logging command window output to: %s\n', logFile);

%% Epoch timing
% t_P9/t_end are protocol landmarks and they are present in
% Config_UserSettings.m (cfg.t_P9_evoked / cfg.t_end_evoked) alongside
% every other per-recording setting, same as fn_getEpochTiming does for
% the 12min/20min protocols in Run_Attractor_Analysis.m.
t_P9  = cfg.t_P9_evoked;
t_end = cfg.t_end_evoked;

t_evoked_start = t_P9 + cfg.motor_buffer_s;
wins_s     = [0, t_P9; t_evoked_start, t_end];
win_labels = {'Baseline', 'Evoked'};
protocol_label = 'evoked';

fprintf('Epoch definitions:\n');
fprintf('  Baseline : 0 - %.0f s\n', t_P9);
fprintf('  Evoked   : %.0f - %.0f s  (P9 + %ds buffer, runs to end)\n', ...
    t_evoked_start, t_end, cfg.motor_buffer_s);

win_f  = round(cfg.slide_win_s  * cfg.fs);
step_f = round(cfg.slide_step_s * cfg.fs);

%% =========================================================================
%  SECTION 1 - LOAD + CLEAN + SMOOTH
%  =========================================================================
fprintf('\nSection 1: Loading and preprocessing...\n');

[peaks, nNeurons, nFrames] = fn_loadPeaksData(cfg.DATA_FILE, cfg.chunk_size);

% fn_loadPeaksData returns the file's full length, and nothing downstream
% crops to t_end [If that is not actually the end of the recording]

nFrames_full = nFrames;
frames_keep  = min(round(t_end * cfg.fs), nFrames_full);
if frames_keep < nFrames_full
    fprintf('  Cropping loaded file: %.1fs (%d frames) on disk -> using first %.1fs (%d frames) per t_end\n', ...
        nFrames_full/cfg.fs, nFrames_full, t_end, frames_keep);
    peaks   = peaks(:, 1:frames_keep);
    nFrames = frames_keep;
elseif frames_keep > nFrames_full
    t_end_old = t_end;
    t_end = nFrames_full / cfg.fs;
    warning('Run_Evoked_Analysis:TEndBeyondFile', ...
        't_end=%.1fs exceeds this file''s actual length (%.1fs) -- using the file''s full length instead.', ...
        t_end_old, t_end);
    wins_s(end,2) = t_end;
end

[peaks, good, rate_all, nN] = fn_qualityFilterNeurons(peaks, cfg.fs, cfg.min_rate, cfg.max_rate); %#ok<ASGLU>
fn_plotRateDistribution(rate_all, cfg.recording_ID, nNeurons, cfg.FIGURES_DIR);

fn_plotRaster_evoked(peaks, cfg.fs, nN, nFrames, t_P9, t_end, cfg.recording_ID, protocol_label, cfg.FIGURES_DIR);

sigma_s = fn_estimateKernelWidth(peaks, cfg.fs);
spike_conv = fn_gaussConvSpikes(peaks, sigma_s, cfg.fs);
clear peaks;
fprintf('  spike_conv: %.0f MB (single)\n', nN*nFrames*4/1e6);

fn_plotSmoothedActivity_evoked(spike_conv, cfg.fs, nN, nFrames, sigma_s, t_P9, t_end, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 2 - PCA + EPOCH CENTROIDS
%  =========================================================================
fprintf('\nSection 2: PCA...\n');

f_pca_start = round(t_evoked_start * cfg.fs);
[V, var_exp, mu_pca] = fn_computePCA(spike_conv, f_pca_start, nFrames, cfg.pca_block_size);
fprintf('  PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%\n', var_exp(1), var_exp(2), var_exp(3));
fn_plotScree(var_exp, nN, cfg.FIGURES_DIR);

[epoch_scores, centroids] = fn_epochScoresAndCentroids(spike_conv, wins_s, cfg.fs, V, mu_pca, win_labels); %#ok<ASGLU>

% fn_plotTrajectory3D splits one continuous trajectory at C2; it doesn't
% fit two separate Baseline/Evoked segments, so fn_plotTrajectory3D_evoked
% (below) builds the plot directly, reusing epoch_scores from above.
fn_plotTrajectory3D_evoked(epoch_scores, win_labels, var_exp, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 3 - WINDOWED PCA + PARTICIPATION RATIO
%  =========================================================================
fprintf('\nSection 3: Windowed PCA + participation ratio...\n');

win_frames = max(min(round(wins_s * cfg.fs), nFrames), 1);
WPCA = fn_windowedPCA(spike_conv, win_frames, []);
for w = 1:numel(WPCA)
    fprintf('  [%s]  top PC = %.1f%%  |  dims@95%% = %d\n', ...
        win_labels{w}, WPCA(w).var_explained(1), WPCA(w).n_dims_95pct);
end

PR_raw  = arrayfun(@(w) fn_participationRatio(WPCA(w).eigenvalues), 1:numel(WPCA));
PR_norm = PR_raw / nN;
fprintf('  PR/N - Baseline: %.3f  Evoked: %.3f\n', PR_norm);

[pr_t, pr_v] = fn_slidingParticipationRatio(spike_conv, win_f, step_f, cfg.fs, nN);

fn_plotPRTimeseries_evoked(pr_t, pr_v, t_P9, t_evoked_start, t_end, ...
    cfg.motor_buffer_s, cfg.slide_win_s, cfg.slide_step_s, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 4 - SUBSPACE ALIGNMENT (Baseline -> Evoked)
%  =========================================================================
fprintf('\nSection 4: Subspace alignment...\n');

[nDims, chance_lvl] = fn_chooseAlignmentK(cfg.K_ALIGN, WPCA(2), nN, cfg.align_var_thresh_pct);

align_mat = fn_epochAlignmentMatrix(WPCA, nDims);   % 2x2
align_BE  = align_mat(1,2);
align_BE_corrected = (align_BE - chance_lvl) / (1 - chance_lvl);
fprintf('  Baseline->Evoked = %.3f  (chance=%.3f, ratio=%.1fx)\n', align_BE, chance_lvl, align_BE/chance_lvl);
fprintf('  Baseline->Evoked (chance-corrected) = %.3f\n', align_BE_corrected);

fn_plotAlignmentMatrix(align_mat, win_labels, cfg.recording_ID, cfg.FIGURES_DIR);

K_range = max(2,floor(nDims/2)):min(floor(3*nDims/2), nN);
al_vs_k = fn_alignmentVsK(WPCA(1), WPCA(2), K_range);
fn_plotAlignmentVsK(K_range, al_vs_k, nN, nDims, cfg.recording_ID, cfg.FIGURES_DIR);

evoked_axes = WPCA(2).eigenvectors(:,1:nDims);
[al_t, al_v] = fn_slidingSubspaceAlignment(spike_conv, evoked_axes, win_f, step_f, cfg.fs, nDims);

fn_plotAlignmentTimeseries_evoked(al_t, al_v, chance_lvl, t_P9, t_evoked_start, t_end, ...
    cfg.motor_buffer_s, nDims, cfg.slide_win_s, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5a - RECURRENCE DENSITY ANALYSIS (Evoked only)
%  =========================================================================
fprintf('\nSection 5a: Recurrence density analysis...\n');

mu_global = mean(spike_conv, 2);

% MAX_FULL_PTS and EPS_PCTILE come from cfg, same fields Run_Attractor_Analysis.m
[traj_ep_global, ~] = fn_getEpochTrajectory(spike_conv, round(t_evoked_start*cfg.fs), round(t_end*cfg.fs), ...
    cfg.fs, V, mu_global, cfg.MAX_FULL_PTS);
epsilon_rr = prctile(pdist(traj_ep_global,'euclidean'), cfg.EPS_PCTILE);
fprintf('  eps = %.4f  (%gth pct of evoked pairwise distances)\n', epsilon_rr, cfg.EPS_PCTILE);

[traj_full, t_full] = fn_getEpochTrajectory(spike_conv, 1, nFrames, cfg.fs, V, mu_global, cfg.MAX_FULL_PTS);
[~, ~, recurs_full, tested_full] = fn_recurrenceDensity(traj_full, t_full, epsilon_rr, cfg.MIN_LAG_S, cfg.slide_win_s);
[rr_full_t, rr_full_v] = fn_slidingRecurrenceDensity(t_full, recurs_full, tested_full, 1, nFrames, cfg.fs, win_f, step_f); %#ok<ASGLU>

[traj_ev, t_ev_ax] = fn_getEpochTrajectory(spike_conv, round(t_evoked_start*cfg.fs), round(t_end*cfg.fs), ...
    cfg.fs, V, mu_global, cfg.MAX_WIN_PTS);
buf_ev = min(cfg.slide_win_s, 0.25*(t_ev_ax(end)-t_ev_ax(1)));

[ev_recur_density, ~, recurs_ev, ~] = fn_recurrenceDensity(traj_ev, t_ev_ax, epsilon_rr, cfg.MIN_LAG_S, buf_ev);
fprintf('  Evoked recurrence density (own eps): %.3f\n', ev_recur_density);

[R_ev, ~, ~] = fn_recurrencePlot(traj_ev, epsilon_rr, 'euc');
RQA_ev.RR = ev_recur_density; %#ok<STRNU>

% Baseline self-recurrence + Baseline-vs-Evoked cross-recurrence, using
% the evoked-calibrated epsilon so both sit on the same scale as
% ev_recur_density. Feeds panel 3 of the attractor-comparison figure.
[traj_base, t_base_ax] = fn_getEpochTrajectory(spike_conv, 1, round(t_P9*cfg.fs), cfg.fs, V, mu_global, cfg.MAX_WIN_PTS);
buf_base = min(cfg.slide_win_s, 0.25*(t_base_ax(end)-t_base_ax(1)));
[base_recur_density, ~, ~, ~] = fn_recurrenceDensity(traj_base, t_base_ax, epsilon_rr, cfg.MIN_LAG_S, buf_base);
fprintf('  Baseline recurrence density (evoked eps): %.3f\n', base_recur_density);

[cross_recur_density, recurs_base_in_ev, ~] = fn_crossRecurrenceDensity(traj_base, traj_ev, epsilon_rr); %#ok<ASGLU>
fprintf('  Cross-recurrence density (Baseline found in Evoked): %.3f\n', cross_recur_density);

% No Recovery epoch exists, so fn_plotRecurrenceSummary is called with
% traj_ev/t_ev_ax standing in for Recovery (well-formed dummy data);
% fn_plotRecurrenceSummary_evoked then deletes the resulting Recovery /
% cross-recurrence panels and replaces them with a "No Recovery epoch"
% label.
fn_plotRecurrenceSummary_evoked(traj_ev, t_ev_ax, R_ev, ev_recur_density, ...
    epsilon_rr, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5b - ATTRACTOR ONSET DETECTION (P9 -> end)
%  =========================================================================
fprintf('\nSection 5b: Attractor onset detection...\n');

onset = fn_detectAttractorEpoch(spike_conv, cfg.fs, V, mu_global, t_P9, t_end, ...
    epsilon_rr, cfg.MIN_LAG_S, cfg.slide_win_s, win_f, step_f, cfg.ONSET_RATIO, cfg.MAX_FULL_PTS, t_evoked_start);
t_attractor_onset = onset.t_lock;
onset_detected     = onset.detected;
if onset_detected
    fprintf('  Attractor ONSET: t=%.1f s (%.1f s after P9)\n', t_attractor_onset, t_attractor_onset - t_P9);
else
    fprintf('  Attractor ONSET: not detected -- fallback to %.0f s\n', t_attractor_onset);
end

fprintf('\n-- Epoch summary --\n');
idx_ev_lock = onset.win_t >= t_attractor_onset & onset.win_t < t_end;
if any(idx_ev_lock)
    fprintf('Evoked (%.1f-%.1f s): recurrence density = %.1f%% +/- %.1f%%\n', ...
        t_attractor_onset, t_end, nanmean(onset.win_v(idx_ev_lock))*100, nanstd(onset.win_v(idx_ev_lock))*100);
end

% fn_plotEpochDetection requires a return_ struct; return_placeholder has
% detected=false (skips Recovery shading) and empty win_t/win_v/
% win_v_fixed so its unconditional plot() calls are no-ops.
return_placeholder.detected    = false;
return_placeholder.t_lock      = NaN;
return_placeholder.win_t       = [];
return_placeholder.win_v       = [];
return_placeholder.win_v_fixed = [];
fn_plotEpochDetection_evoked(onset, return_placeholder, t_P9, t_end, cfg.ONSET_RATIO, cfg.recording_ID, cfg.FIGURES_DIR);

onset_win_t = onset.win_t; onset_win_v = onset.win_v; onset_win_v_fixed = onset.win_v_fixed; %#ok<NASGU>
eps_onset = onset.epsilon_within; %#ok<NASGU>

%% =========================================================================
%  SECTION 5c - PR & ALIGNMENT ON THE RR-DEFINED EVOKED WINDOW
%  =========================================================================
fprintf('\nSection 5c: PR and alignment on the RR-defined evoked window...\n');

wins_s_rr     = [0, t_P9; t_attractor_onset, t_end];
win_labels_rr = {'Baseline','Attractor-evoked'};

win_frames_rr = max(min(round(wins_s_rr * cfg.fs), nFrames), 1);
WPCA_rr       = fn_windowedPCA(spike_conv, win_frames_rr, []);

PR_raw_rr  = arrayfun(@(w) fn_participationRatio(WPCA_rr(w).eigenvalues), 1:numel(WPCA_rr));
PR_norm_rr = PR_raw_rr / nN;
fprintf('  PR/N [RR-defined] - Baseline: %.3f  Evoked: %.3f\n', PR_norm_rr);

align_mat_rr = fn_epochAlignmentMatrix(WPCA_rr, nDims);
fprintf('  Alignment B->E [RR-defined]: %.3f  |  [fixed]: %.3f  |  Chance: %.3f\n', ...
    align_mat_rr(1,2), align_mat(1,2), chance_lvl);

% fn_plotAttractorComparison hard-codes 3 epochs in panel 1's XTick, so
% fn_plotAttractorComparison_evoked pads Evoked as a placeholder 3rd bar
% there and removes it afterward. Panel 3 needs no padding -- it already
% has 3 real values (Evoked, Baseline, cross-recurrence) -- only its tick
% labels are corrected.
fn_plotAttractorComparison_evoked(PR_norm, PR_norm_rr, win_labels, ...
    align_mat(1,2), align_mat_rr(1,2), chance_lvl, ...
    ev_recur_density, base_recur_density, cross_recur_density, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SAVE
%  =========================================================================
save_path   = fullfile(cfg.RESULTS_DIR, sprintf('%s_results.mat', cfg.recording_ID));
nN_saved    = nN;
nDims_align = nDims;
recording_ID = cfg.recording_ID; fs = cfg.fs; %#ok<NASGU>
ONSET_RATIO = cfg.ONSET_RATIO; %#ok<NASGU>

save(save_path, ...
    'recording_ID','fs', ...
    't_P9','t_end', ...
    'wins_s','wins_s_rr','win_labels','win_labels_rr','nN','nN_saved', ...
    't_attractor_onset','onset_detected', ...
    'ONSET_RATIO', ...
    'WPCA','PR_raw','PR_norm','pr_t','pr_v', ...
    'WPCA_rr','PR_raw_rr','PR_norm_rr', ...
    'align_mat','align_mat_rr','al_t','al_v','nDims','nDims_align','chance_lvl', ...
    'align_BE','align_BE_corrected', ...
    'RQA_ev', ...
    'epsilon_rr', ...
    'rr_full_t','rr_full_v', ...
    'onset_win_t','onset_win_v','onset_win_v_fixed', ...
    '-v7.3');

fprintf('\nSaved: %s\nDone.\n', save_path);
diary off;

%% =========================================================================
%  LOCAL FUNCTIONS
%  Each wrapper calls the original Functions/Plotting/ function
%  unmodified, then edits the resulting figure before re-exporting over
%  the same file path. Functions/ is never touched.
% =========================================================================

function fn_plotTrajectory3D_evoked(epoch_scores, win_labels, var_exp, recording_ID, figuresDir)
% Baseline/Evoked trajectories with a P9 marker at their boundary.
% epoch_scores{1}/{2} are already PCA-projected by
% fn_epochScoresAndCentroids; the P9 point is Baseline's last sample
% (Baseline window = [0, t_P9]).
baseline = epoch_scores{1};
evoked   = epoch_scores{2};
p9_point = baseline(end, :);

fig = figure('Name', 'Population trajectory', 'Position', [200 100 900 800]);
plot3(baseline(:,1), baseline(:,2), baseline(:,3), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.0); hold on;
plot3(evoked(:,1), evoked(:,2), evoked(:,3), 'r-', 'LineWidth', 1.2);
scatter3(p9_point(1), p9_point(2), p9_point(3), 100, 'gs', 'filled');

xlabel(sprintf('PC1 (%.1f%%)', var_exp(1)));
ylabel(sprintf('PC2 (%.1f%%)', var_exp(2)));
zlabel(sprintf('PC3 (%.1f%%)', var_exp(3)));
legend({win_labels{1}, win_labels{2}, 'P9 onset'}, 'Location', 'best');
grid on; view([-35 25]);
title(sprintf('Population trajectory | %s', recording_ID));
exportgraphics(fig, fullfile(figuresDir, '02_trajectory_3D.png'), 'Resolution', 500);
end

function fn_deleteXlinesAtValue(fig, targetValue)
% Deletes ConstantLine (xline/yline) objects matching targetValue. Uses
% findall, not findobj -- fn_plotEpochDetection draws its xlines with
% HandleVisibility='off', which findobj silently skips.
allLines = findall(fig, 'Type', 'ConstantLine');
for k = 1:numel(allLines)
    if isequal(allLines(k).Value, targetValue)
        delete(allLines(k));
    end
end
end

function fn_plotRaster_evoked(peaks, fs, nN, nFrames, t_P9, t_end, recording_ID, protocol_label, figuresDir)
fn_plotRaster(peaks, fs, nN, nFrames, t_P9, t_end, recording_ID, protocol_label, figuresDir);
fig = gcf;
fn_deleteXlinesAtValue(fig, t_end);
exportgraphics(fig, fullfile(figuresDir, '01b_raster.png'), 'Resolution', 500);
end

function fn_plotSmoothedActivity_evoked(spike_conv, fs, nN, nFrames, sigma_s, t_P9, t_end, recording_ID, figuresDir)
fn_plotSmoothedActivity(spike_conv, fs, nN, nFrames, sigma_s, t_P9, t_end, recording_ID, figuresDir);
fig = gcf;
fn_deleteXlinesAtValue(fig, t_end);
exportgraphics(fig, fullfile(figuresDir, '01c_smoothed.png'), 'Resolution', 500);
end

function fn_plotPRTimeseries_evoked(pr_t, pr_v, t_P9, t_evoked_start, t_end, motor_buffer_s, slide_win_s, slide_step_s, recording_ID, figuresDir)
% t_end stands in for both t_C2 and t_recovery_start (both xlines land at
% the same value and are deleted together). No explicit xlim here, so
% re-autoscale after deleting the marker.
fn_plotPRTimeseries(pr_t, pr_v, t_P9, t_evoked_start, t_end, t_end, motor_buffer_s, 0, slide_win_s, slide_step_s, recording_ID, figuresDir);
fig = gcf;
fn_deleteXlinesAtValue(fig, t_end);
xlim(findobj(fig, 'Type', 'axes'), 'auto');
exportgraphics(fig, fullfile(figuresDir, '03_PR_timeseries.png'), 'Resolution', 500);
end

function fn_plotAlignmentTimeseries_evoked(al_t, al_v, chance_lvl, t_P9, t_evoked_start, t_end, motor_buffer_s, nDims, slide_win_s, recording_ID, figuresDir)
% Same pattern as fn_plotPRTimeseries_evoked.
fn_plotAlignmentTimeseries(al_t, al_v, chance_lvl, t_P9, t_evoked_start, t_end, t_end, motor_buffer_s, 0, nDims, slide_win_s, recording_ID, figuresDir);
fig = gcf;
fn_deleteXlinesAtValue(fig, t_end);
xlim(findobj(fig, 'Type', 'axes'), 'auto');
exportgraphics(fig, fullfile(figuresDir, '04c_alignment_timeseries.png'), 'Resolution', 500);
end

function fn_plotEpochDetection_evoked(onset, return_placeholder, t_P9, t_end, onsetRatio, recording_ID, figuresDir)
% t_end stands in for t_C2 (removed below) and is also the Evoked
% shading's real right edge (xlim is explicit, so the deletion doesn't
% affect the axis). return_placeholder carries no data, so its "C2->end"
% line is empty; DisplayNames are corrected since the search actually ran
% P9->end, not P9->C2.
fn_plotEpochDetection(onset, return_placeholder, t_P9, t_end, t_end, onsetRatio, recording_ID, figuresDir);
fig = gcf;
fn_deleteXlinesAtValue(fig, t_end);

ax = findobj(fig, 'Type', 'axes');
allLines = findobj(ax, 'Type', 'Line');
for k = 1:numel(allLines)
    c  = allLines(k).Color;
    ls = allLines(k).LineStyle;
    if isequal(c, [0 0 1])
        delete(allLines(k));
    elseif isequal(c, [1 0 0]) && strcmp(ls, '-')
        allLines(k).DisplayName = 'P9->end (within-epoch eps)';
    elseif isequal(c, [1 0 0]) && strcmp(ls, ':')
        allLines(k).DisplayName = 'P9->end (fixed evoked eps, comparison)';
    end
end
legend(ax, 'Location', 'northeastoutside');

exportgraphics(fig, fullfile(figuresDir, '05b_epoch_detection.png'), 'Resolution', 550);
end

function fn_plotRecurrenceSummary_evoked(traj_ev, t_ev_ax, R_ev, ev_recur_density, epsilon_rr, recording_ID, figuresDir)
% Calls the original with traj_ev/t_ev_ax/R_ev reused as the stand-in
% Recovery arguments, purely so its internal sizing/loop logic runs on
% well-formed inputs. Panels 2-4 (Recovery / cross-recurrence /
% evoked-in-recovery) are then replaced with a "not applicable" label.
% Panel 1 (Evoked) is real throughout.
fn_plotRecurrenceSummary(traj_ev, t_ev_ax, traj_ev, t_ev_ax, R_ev, R_ev, ...
    ev_recur_density, NaN, NaN, false(size(t_ev_ax)), epsilon_rr, recording_ID, figuresDir);

fig = gcf;
ax_all = findobj(fig, 'Type', 'axes');
pos = cell2mat(get(ax_all, 'Position'));
[~, order] = sort(pos(:,1));
ax_all = ax_all(order);

for k = 2:4
    cla(ax_all(k));
    axis(ax_all(k), 'off');
    text(ax_all(k), 0.5, 0.5, 'No Recovery epoch', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center', ...
        'FontSize', 11, 'Color', [0.4 0.4 0.4]);
end

exportgraphics(fig, fullfile(figuresDir, '05a_recurrence_all.png'), 'Resolution', 500);
end

function fn_plotAttractorComparison_evoked(PR_norm, PR_norm_rr, win_labels, ...
    align_ER_fixed, align_ER_rr, chance_lvl, ev_recur_density, base_recur_density, cross_recur_density, recording_ID, figuresDir)
% Panel 1's XTick=1:3 layout is hard-coded, so a placeholder 3rd "epoch"
% (Evoked's PR duplicated, label '(n/a)') is padded in and removed after.
% Panel 3 needs no padding -- its 3 bars are all real -- only its tick
% labels are corrected. Panel 2 is untouched.
PR_norm_padded    = [PR_norm, PR_norm(end)];
PR_norm_rr_padded = [PR_norm_rr, PR_norm_rr(end)];
win_labels_padded = [win_labels, {'(n/a)'}];

fn_plotAttractorComparison(PR_norm_padded, PR_norm_rr_padded, win_labels_padded, ...
    align_ER_fixed, align_ER_rr, chance_lvl, ...
    ev_recur_density, base_recur_density, cross_recur_density, recording_ID, figuresDir);

fig = gcf;
ax_all = findobj(fig, 'Type', 'axes');
pos = cell2mat(get(ax_all, 'Position'));
[~, order] = sort(pos(:,1));
ax_all = ax_all(order);

bars1 = findobj(ax_all(1), 'Type', 'Bar');
for k = 1:numel(bars1)
    yd = bars1(k).YData; yd(3) = NaN; bars1(k).YData = yd;
end
set(ax_all(1), 'XTick', 1:2, 'XTickLabel', win_labels, 'XLim', [0.4 2.6]);

set(ax_all(3), 'XTickLabel', {'Evoked', 'Baseline', sprintf('Cross\n(Base->Ev)')});
title(ax_all(3), 'Baseline vs Evoked recurrence density', 'FontWeight', 'normal');

exportgraphics(fig, fullfile(figuresDir, '05c_attractor_comparison.png'), 'Resolution', 550);
end
