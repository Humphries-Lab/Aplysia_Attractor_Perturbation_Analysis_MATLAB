%% RUN_ATTRACTOR_ANALYSIS
clear; close all; clc;

thisFile = mfilename('fullpath');
if contains(thisFile, fullfile(tempdir))
    toolboxRoot = pwd;   % assumes MATLAB's current folder is set to the toolbox root
else
    toolboxRoot = fileparts(thisFile);
end

addpath(genpath(fullfile(toolboxRoot, 'Functions')));

cfgFile = fullfile(toolboxRoot, 'Config_UserSettings.m');
if ~isfile(cfgFile)
    error('Run_Attractor_Analysis:ConfigNotFound', ...
        ['Config_UserSettings.m not found at:\n  %s\n' ...
         'toolboxRoot resolved to:\n  %s\n' ...
         'Check MATLAB''s Current Folder is set to the toolbox root, ' ...
         'or run the full script (not "Run Section").'], cfgFile, toolboxRoot);
end
run(cfgFile);   % defines `cfg`

if ~exist(cfg.RESULTS_DIR,'dir'), mkdir(cfg.RESULTS_DIR); end
if ~exist(cfg.FIGURES_DIR,'dir'), mkdir(cfg.FIGURES_DIR); end

%% =========================================================================
%  Epoch timing 
%  =========================================================================
timing = fn_getEpochTiming(cfg.protocol, cfg.fs, cfg.motor_buffer_s, cfg.recovery_delay_s);
t_P9 = timing.t_P9; t_C2 = timing.t_C2; t_end = timing.t_end;
t_evoked_start   = timing.t_evoked_start;
t_recovery_start = timing.t_recovery_start;
wins_s     = timing.wins_s;
win_labels = timing.win_labels;

win_f  = round(cfg.slide_win_s  * cfg.fs);
step_f = round(cfg.slide_step_s * cfg.fs);

%% =========================================================================
%  SECTION 1 - LOAD + CLEAN + SMOOTH
%  =========================================================================
fprintf('\nSection 1: Loading and preprocessing...\n');

[peaks, nNeurons, nFrames] = fn_loadPeaksData(cfg.DATA_FILE, cfg.chunk_size);

[peaks, good, rate_all, nN] = fn_qualityFilterNeurons(peaks, cfg.fs, cfg.min_rate, cfg.max_rate); %#ok<ASGLU>
fn_plotRateDistribution(rate_all, cfg.recording_ID, nNeurons, cfg.FIGURES_DIR);
fn_plotRaster(peaks, cfg.fs, nN, nFrames, t_P9, t_C2, cfg.recording_ID, cfg.protocol, cfg.FIGURES_DIR);

sigma_s = fn_estimateKernelWidth(peaks, cfg.fs);
spike_conv = fn_gaussConvSpikes(peaks, sigma_s, cfg.fs);
clear peaks;
fprintf('  spike_conv: %.0f MB (single)\n', nN*nFrames*4/1e6);

fn_plotSmoothedActivity(spike_conv, cfg.fs, nN, nFrames, sigma_s, t_P9, t_C2, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 2 - PCA + 3-D POPULATION TRAJECTORY
%  =========================================================================
fprintf('\nSection 2: PCA...\n');

f_pca_start = round(t_evoked_start * cfg.fs);
[V, var_exp, mu_pca] = fn_computePCA(spike_conv, f_pca_start, nFrames, cfg.pca_block_size);
fprintf('  PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%\n', var_exp(1), var_exp(2), var_exp(3));

scores  = fn_projectOntoPCs(spike_conv, f_pca_start, nFrames, V, mu_pca, 3, cfg.pca_block_size);
t_pca   = (f_pca_start:nFrames-1) / cfg.fs;
fn_plotScree(var_exp, nN, cfg.FIGURES_DIR);
fn_plotTrajectory3D(scores, t_pca, t_C2, var_exp, cfg.recording_ID, cfg.FIGURES_DIR);

[epoch_scores, centroids] = fn_epochScoresAndCentroids(spike_conv, wins_s, cfg.fs, V, mu_pca, win_labels); %#ok<ASGLU>

%% =========================================================================
%  SECTION 2b - WINDOWED PCA
%  =========================================================================
fprintf('\nSection 3: Windowed PCA + participation ratio...\n');

win_frames = max(min(round(wins_s * cfg.fs), nFrames), 1);
WPCA = fn_windowedPCA(spike_conv, win_frames, []);
for w = 1:numel(WPCA)
    fprintf('  [%s]  top PC = %.1f%%  |  dims@95%% = %d\n', ...
        win_labels{w}, WPCA(w).var_explained(1), WPCA(w).n_dims_95pct);
end

%% =========================================================================
%  SECTION 3 - DIMENSIONALITY (PARTICIPATION RATIO)
%  =========================================================================

PR_raw  = arrayfun(@(w) fn_participationRatio(WPCA(w).eigenvalues), 1:numel(WPCA));
PR_norm = PR_raw / nN;
fprintf('  PR/N - Baseline: %.3f  Evoked: %.3f  Recovery: %.3f\n', PR_norm);

[pr_t, pr_v] = fn_slidingParticipationRatio(spike_conv, win_f, step_f, cfg.fs, nN);
fn_plotPRTimeseries(pr_t, pr_v, t_P9, t_evoked_start, t_C2, t_recovery_start, ...
    cfg.motor_buffer_s, cfg.recovery_delay_s, cfg.slide_win_s, cfg.slide_step_s, ...
    cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 4 - SUBSPACE ALIGNMENT
%
%  A(V1,V2) = (1/K) trace(Q2'Q1 Q1'Q2) = mean squared cosine of principal
%  angles between two K-dim subspaces. Chance level for random K-frames in
%  N dims is K/N. See Cunningham & Yu (2014) Nat Neurosci; Semedo et al.
%  (2019) Neuron. K is chosen data-driven, as the number of Evoked-epoch
%  dimensions explaining cfg.align_var_thresh_pct percent variance.
%  =========================================================================
fprintf('\nSection 4: Subspace alignment...\n');

[nDims, chance_lvl] = fn_chooseAlignmentK(cfg.K_ALIGN, WPCA(2), nN, cfg.align_var_thresh_pct);

align_mat = fn_epochAlignmentMatrix(WPCA, nDims);
align_ER  = align_mat(2,3);
align_ER_corrected = (align_ER - chance_lvl) / (1 - chance_lvl);
fprintf('  Evoked->Recovery = %.3f  (chance=%.3f, ratio=%.1fx)\n', align_ER, chance_lvl, align_ER/chance_lvl);
fprintf('  Evoked->Recovery (chance-corrected) = %.3f\n', align_ER_corrected);

fn_plotAlignmentMatrix(align_mat, win_labels, cfg.recording_ID, cfg.FIGURES_DIR);

K_range = max(2,floor(nDims/2)):min(floor(3*nDims/2), nN);
al_vs_k = fn_alignmentVsK(WPCA(2), WPCA(3), K_range);
fn_plotAlignmentVsK(K_range, al_vs_k, nN, nDims, cfg.recording_ID, cfg.FIGURES_DIR);

evoked_axes = WPCA(2).eigenvectors(:,1:nDims);
[al_t, al_v] = fn_slidingSubspaceAlignment(spike_conv, evoked_axes, win_f, step_f, cfg.fs, nDims);
fn_plotAlignmentTimeseries(al_t, al_v, chance_lvl, t_P9, t_evoked_start, t_C2, t_recovery_start, ...
    cfg.motor_buffer_s, cfg.recovery_delay_s, nDims, cfg.slide_win_s, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5 - RECURRENCE ANALYSIS (fixed-window epochs)
%  =========================================================================
fprintf('\nSection 5: Recurrence density analysis...\n');

mu_global = mean(spike_conv, 2);

% Calibrate epsilon from the Evoked epoch (density-matched percentile)
[traj_ep_global, ~] = fn_getEpochTrajectory(spike_conv, round(t_evoked_start*cfg.fs), round(t_C2*cfg.fs), ...
    cfg.fs, V, mu_global, 1500);
epsilon_rr = prctile(pdist(traj_ep_global,'euclidean'), 10);
fprintf('  eps = %.4f  (10th pct of evoked pairwise distances)\n', epsilon_rr);

% Full-recording per-point recurrence, binned into the sliding-window grid
[traj_full, t_full] = fn_getEpochTrajectory(spike_conv, 1, nFrames, cfg.fs, V, mu_global, cfg.MAX_FULL_PTS);
fprintf('  Computing per-point recurrence density over full recording (%d pts)...\n', numel(t_full));
[~, ~, recurs_full, tested_full] = fn_recurrenceDensity(traj_full, t_full, epsilon_rr, cfg.MIN_LAG_S, cfg.slide_win_s);
[rr_full_t, rr_full_v] = fn_slidingRecurrenceDensity(t_full, recurs_full, tested_full, 1, nFrames, cfg.fs, win_f, step_f);
fprintf('  Computed %d recurrence-density windows.\n', numel(rr_full_t));

% Evoked / recovery trajectories for self- and cross-recurrence
[traj_ev, t_ev_ax] = fn_getEpochTrajectory(spike_conv, round(t_evoked_start*cfg.fs), round(t_C2*cfg.fs), ...
    cfg.fs, V, mu_global, cfg.MAX_WIN_PTS);
[traj_re, t_re_ax] = fn_getEpochTrajectory(spike_conv, round(t_recovery_start*cfg.fs), nFrames, ...
    cfg.fs, V, mu_global, cfg.MAX_WIN_PTS);

buf_ev = min(cfg.slide_win_s, 0.25*(t_ev_ax(end)-t_ev_ax(1)));
buf_re = min(cfg.slide_win_s, 0.25*(t_re_ax(end)-t_re_ax(1)));

% Self recurrence density: each epoch uses ITS OWN calibrated epsilon
% ("is this epoch internally on an attractor" is a per-epoch question).
epsilon_re = fn_safeEpsilon(traj_re);
fprintf('  eps: evoked(own)=%.4f  recovery(own)=%.4f\n', epsilon_rr, epsilon_re);

[ev_recur_density, ~, recurs_ev, ~] = fn_recurrenceDensity(traj_ev, t_ev_ax, epsilon_rr, cfg.MIN_LAG_S, buf_ev); %#ok<ASGLU>
[re_recur_density, ~, recurs_re, ~] = fn_recurrenceDensity(traj_re, t_re_ax, epsilon_re, cfg.MIN_LAG_S, buf_re); %#ok<ASGLU>
fprintf('  Evoked   recurrence density (own eps): %.3f\n', ev_recur_density);
fprintf('  Recovery recurrence density (own eps): %.3f\n', re_recur_density);

[R_ev, ~, ~] = fn_recurrencePlot(traj_ev, epsilon_rr, 'euc');
[R_re, ~, ~] = fn_recurrencePlot(traj_re, epsilon_re, 'euc');

% Cross-recurrence: does each evoked point recur ANYWHERE in recovery?
% (Shared epsilon, since this asks a "same reference scale" question.)
fprintf('  Computing cross-recurrence density (evoked recurring in recovery)...\n');
[cross_recur_density, recurs_ev_in_re, ~] = fn_crossRecurrenceDensity(traj_ev, traj_re, epsilon_rr);
fprintf('  Cross-recurrence density (evoked found in recovery): %.3f\n', cross_recur_density);

RQA_ev.RR = ev_recur_density; %#ok<STRNU>
RQA_re.RR = re_recur_density; %#ok<STRNU>
cross_rr_mean = cross_recur_density; %#ok<NASGU>

fn_plotRecurrenceSummary(traj_ev, t_ev_ax, traj_re, t_re_ax, R_ev, R_re, ...
    ev_recur_density, re_recur_density, cross_recur_density, recurs_ev_in_re, ...
    epsilon_rr, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5c - EPOCH DETECTION VIA RECURRENCE DENSITY
%  Detect the first time the population "locks onto" an attractor,
%  searching P9->C2 for onset and C2->end for return, each using a
%  within-search-window calibrated epsilon (fixed evoked eps kept only as
%  a plotted comparison).
%  =========================================================================
fprintf('\nSection 5c: Attractor epoch detection based on recurrence density...\n');

onset = fn_detectAttractorEpoch(spike_conv, cfg.fs, V, mu_global, t_P9, t_C2, ...
    epsilon_rr, cfg.MIN_LAG_S, cfg.slide_win_s, win_f, step_f, cfg.ONSET_RATIO, cfg.MAX_FULL_PTS, t_evoked_start);
t_attractor_onset = onset.t_lock;
onset_detected     = onset.detected;
if onset_detected
    fprintf('  Attractor ONSET: t=%.1f s (%.1f s after P9)\n', t_attractor_onset, t_attractor_onset - t_P9);
else
    fprintf('  Attractor ONSET: not detected -- fallback to %.0f s\n', t_attractor_onset);
end

return_ = fn_detectAttractorEpoch(spike_conv, cfg.fs, V, mu_global, t_C2, t_end, ...
    epsilon_rr, cfg.MIN_LAG_S, cfg.slide_win_s, win_f, step_f, cfg.ONSET_RATIO, cfg.MAX_FULL_PTS, ...
    t_C2 + (cfg.slide_win_s/2));
t_attractor_return = return_.t_lock;
return_detected     = return_.detected;
if return_detected
    fprintf('  Attractor RETURN: t=%.1f s (%.1f s after C2)\n', t_attractor_return, t_attractor_return - t_C2);
else
    fprintf('  Attractor RETURN: not detected -- fallback to %.0f s\n', t_attractor_return);
end

fprintf('\n-- Epoch summary --\n');
idx_ev_lock  = onset.win_t  >= t_attractor_onset  & onset.win_t  < t_C2;
idx_rec_lock = return_.win_t >= t_attractor_return;
if any(idx_ev_lock)
    fprintf('Evoked  (%.1f-%.1f s): recurrence density = %.1f%% +/- %.1f%%\n', ...
        t_attractor_onset, t_C2, nanmean(onset.win_v(idx_ev_lock))*100, nanstd(onset.win_v(idx_ev_lock))*100);
end
if any(idx_rec_lock)
    fprintf('Recovery (%.1f-end):   recurrence density = %.1f%% +/- %.1f%%\n', ...
        t_attractor_return, nanmean(return_.win_v(idx_rec_lock))*100, nanstd(return_.win_v(idx_rec_lock))*100);
end
fprintf('Cross-recurrence density: %.3f\n', cross_recur_density);

fn_plotEpochDetection(onset, return_, t_P9, t_C2, t_end, cfg.ONSET_RATIO, cfg.recording_ID, cfg.FIGURES_DIR);

% Preserve original variable names for the saved results file
onset_win_t = onset.win_t; onset_win_v = onset.win_v; onset_win_v_fixed = onset.win_v_fixed; %#ok<NASGU>
return_win_t = return_.win_t; return_win_v = return_.win_v; return_win_v_fixed = return_.win_v_fixed; %#ok<NASGU>
eps_onset = onset.epsilon_within; eps_return = return_.epsilon_within; %#ok<NASGU>

%% =========================================================================
%  SECTION 5d - PR & SUBSPACE ALIGNMENT ON ATTRACTOR-DEFINED EPOCHS
%  =========================================================================
fprintf('\nSection 5d: PR and alignment on attractor-defined epochs...\n');

wins_s_rr     = [0, t_P9; t_attractor_onset, t_C2; t_attractor_return, t_end];
win_labels_rr = {'Baseline','Attractor-evoked','Attractor-recovery'};

win_frames_rr = max(min(round(wins_s_rr * cfg.fs), nFrames), 1);
WPCA_rr       = fn_windowedPCA(spike_conv, win_frames_rr, []);

PR_raw_rr  = arrayfun(@(w) fn_participationRatio(WPCA_rr(w).eigenvalues), 1:numel(WPCA_rr));
PR_norm_rr = PR_raw_rr / nN;
fprintf('  PR/N [RR-defined] - Baseline: %.3f  Evoked: %.3f  Recovery: %.3f\n', PR_norm_rr);

align_mat_rr = fn_epochAlignmentMatrix(WPCA_rr, nDims);
fprintf('  Alignment E->R [RR-defined]: %.3f  |  [fixed]: %.3f  |  Chance: %.3f\n', ...
    align_mat_rr(2,3), align_mat(2,3), chance_lvl);

fn_plotAttractorComparison(PR_norm, PR_norm_rr, win_labels, ...
    align_mat(2,3), align_mat_rr(2,3), chance_lvl, ...
    ev_recur_density, re_recur_density, cross_recur_density, cfg.recording_ID, cfg.FIGURES_DIR);
    

%% =========================================================================
%  SAVE
%  =========================================================================
save_path   = fullfile(cfg.RESULTS_DIR, sprintf('%s_results.mat', cfg.recording_ID));
nN_saved    = nN;
nDims_align = nDims;
recording_ID = cfg.recording_ID; protocol = cfg.protocol; fs = cfg.fs; %#ok<NASGU>
N_SIGMA = cfg.N_SIGMA; N_CONSEC = cfg.N_CONSEC; ONSET_RATIO = cfg.ONSET_RATIO; %#ok<NASGU>
t_C2_saved = t_C2; %#ok<NASGU>

save(save_path, ...
    'recording_ID','protocol','fs', ...
    't_P9','t_C2','t_C2_saved','t_end','cfg', ...
    'wins_s','wins_s_rr','win_labels','win_labels_rr','nN','nN_saved', ...
    't_attractor_onset','t_attractor_return','onset_detected','return_detected', ...
    'N_SIGMA','N_CONSEC','ONSET_RATIO', ...
    'WPCA','PR_raw','PR_norm','pr_t','pr_v', ...
    'WPCA_rr','PR_raw_rr','PR_norm_rr', ...
    'align_mat','align_mat_rr','al_t','al_v','nDims','nDims_align','chance_lvl', ...
    'align_ER','align_ER_corrected', ...
    'RQA_ev','RQA_re','cross_rr_mean','cross_recur_density','recurs_ev_in_re', ...
    'epsilon_rr','epsilon_re','eps_onset','eps_return', ...
    'rr_full_t','rr_full_v','dom_period', ...
    'onset_win_t','onset_win_v','onset_win_v_fixed', ...
    'return_win_t','return_win_v','return_win_v_fixed', ...
    '-v7.3');

fprintf('\nSaved: %s\nDone.\n', save_path);
