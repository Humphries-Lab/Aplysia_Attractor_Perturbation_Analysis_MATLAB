%% RUN_SENSITIZATION_ANALYSIS
%  Attractor / population-dynamics analysis for a SENSITIZATION protocol:
%  one 25-minute recording built by concatenating five 5-minute files,
%  each containing one P9 stimulus at t=90s (file-relative). This gives
%  1 Baseline and 5 Evoked epochs in a single continuous recording.
%
%  Companion to Run_Evoked_Analysis.m (single-stimulus case) and
%  Run_Attractor_Analysis.m (Baseline+Evoked+Recovery case). Uses the
%  same Functions/ folder and Config_UserSettings.m, unmodified.

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
    error('Run_Sensitization_Analysis:ConfigNotFound', ...
        'Config_UserSettings.m not found at:\n  %s', cfgFile);
end
run(cfgFile);   % defines `cfg`

if ~isfile(cfg.DATA_FILE)
    [~, dataName, dataExt] = fileparts(cfg.DATA_FILE);
    found = fn_findShallow(toolboxRoot, [dataName dataExt]);
    if isempty(found)
        nearby = fn_listMatFilesShallow(toolboxRoot);
        if isempty(nearby)
            listMsg = '(no .mat files found directly under the project folder or its subfolders)';
        else
            listMsg = strjoin(nearby, newline);
        end
        error('Run_Sensitization_Analysis:DataFileNotFound', ...
            ['"%s" not found under %s or its subfolders.\n' ...
             '.mat files that DO exist there:\n%s\n' ...
             'Compare the name above to cfg.DATA_FILE in Config_UserSettings.m -- ' ...
             'this is almost always a spelling/typo mismatch.'], ...
            [dataName dataExt], toolboxRoot, listMsg);
    elseif numel(found) > 1
        warning('Run_Sensitization_Analysis:MultipleMatches', ...
            'Multiple files named "%s" found -- using the first match:\n  %s', ...
            [dataName dataExt], found{1});
    end
    cfg.DATA_FILE = found{1};
end
fprintf('Using DATA_FILE: %s\n', cfg.DATA_FILE);

if ~exist(cfg.RESULTS_DIR,'dir'), mkdir(cfg.RESULTS_DIR); end
cfg.FIGURES_DIR = fullfile(cfg.FIGURES_DIR, cfg.recording_ID);
if ~exist(cfg.FIGURES_DIR,'dir'), mkdir(cfg.FIGURES_DIR); end

% Clear this recording's old PNGs so every figure in the folder comes from THIS run.
% A warning here means the file is locked (open in an image viewer / syncing) and could not be replaced.
oldPngs = dir(fullfile(cfg.FIGURES_DIR, '*.png'));
for k = 1:numel(oldPngs)
    f = fullfile(oldPngs(k).folder, oldPngs(k).name);
    delete(f);
    if isfile(f)
        warning('Run_Sensitization_Analysis:LockedFile', 'Could not delete %s -- close it in any viewer/sync app and re-run.', f);
    end
end

% Command-window logging (same as Run_Evoked_Analysis.m)
if ~exist(cfg.LOGS_DIR,'dir'), mkdir(cfg.LOGS_DIR); end
diary off;   % close any diary left open by a previous run that errored mid-way
logFile = fullfile(cfg.LOGS_DIR, sprintf('%s_log_sensitization.txt', cfg.recording_ID));
diary(logFile);
diary on;
fprintf('Logging command window output to: %s\n', logFile);

%% Epoch timing -- sensitization protocol (edit per recording)
nStim         = 5;      % number of P9 stimuli / sub-files
FILELEN_S     = 300;    % length of each concatenated sub-file (s)
T_P9_INFILE_S = 90;     % P9 onset within each sub-file, file-relative (s)

t_end_nominal = FILELEN_S * nStim;   % expected total length if all nStim stimuli are complete

win_f  = round(cfg.slide_win_s  * cfg.fs);
step_f = round(cfg.slide_step_s * cfg.fs);

%% =========================================================================
%  SECTION 1 - LOAD + CLEAN + SMOOTH
%  ========================================================================
fprintf('\nSection 1: Loading and preprocessing...\n');

[peaks, nNeurons, nFrames] = fn_loadPeaksData(cfg.DATA_FILE, cfg.chunk_size);

nFrames_full = nFrames;
frames_keep  = min(round(t_end_nominal * cfg.fs), nFrames_full);
if frames_keep < nFrames_full
    fprintf('  Cropping loaded file: %.1fs (%d frames) on disk -> using first %.1fs (%d frames) for %d stimuli\n', ...
        nFrames_full/cfg.fs, nFrames_full, t_end_nominal, frames_keep, nStim);
    peaks   = peaks(:, 1:frames_keep);
    nFrames = frames_keep;
elseif frames_keep > nFrames_full
    fprintf('  File on disk (%.1fs) is shorter than the expected %d-stimulus length (%.1fs).\n', ...
        nFrames_full/cfg.fs, nStim, t_end_nominal);
end
t_end = nFrames / cfg.fs;

fileBounds_nominal = 0:FILELEN_S:FILELEN_S*nStim;
t_P9_all_nominal   = fileBounds_nominal(1:end-1) + T_P9_INFILE_S;

nStim_expected = nStim;
nStim = sum(t_P9_all_nominal + cfg.motor_buffer_s < t_end);
if nStim < nStim_expected
    warning('Run_Sensitization_Analysis:TruncatedRecording', ...
        ['Recording (%.1fs) is too short for all %d expected stimuli -- only %d complete ' ...
         'stimulus/stimuli fit. Dropping the incomplete trailing one(s) rather than analysing ' ...
         'an invalid partial window.'], t_end, nStim_expected, nStim);
end
if nStim < 1
    error('Run_Sensitization_Analysis:NoCompleteStimuli', ...
        'Recording (%.1fs) is too short to contain even one complete stimulus (needs >= %.1fs).', ...
        t_end, T_P9_INFILE_S + cfg.motor_buffer_s);
end

fileBounds      = fileBounds_nominal(1:nStim+1);
fileBounds(end) = min(fileBounds(end), t_end);
t_P9_all           = t_P9_all_nominal(1:nStim);
t_evoked_start_all = t_P9_all + cfg.motor_buffer_s;

% ONE genuine baseline, followed by nStim evoked periods
wins_s     = zeros(1+nStim, 2);
win_labels = cell(1, 1+nStim);

wins_s(1, :)   = [fileBounds(1), t_P9_all(1)];
win_labels{1}  = 'Baseline';
for i = 1:nStim
    wins_s(i+1, :)   = [t_evoked_start_all(i), fileBounds(i+1)];
    win_labels{i+1}  = sprintf('Evoked_%d', i);
end
baseline_idx = 1;
evoked_idx   = 2 : (1+nStim);

protocol_label = sprintf('sensitization (1x Baseline, %dx P9)', nStim);

fprintf('Epoch definitions (%d stimuli):\n', nStim);
fprintf('  Baseline: %.0f-%.0f s\n', fileBounds(1), t_P9_all(1));
for i = 1:nStim
    fprintf('  File %d [%.0f-%.0f s] | Evoked_%d: %.0f-%.0f s  (P9 @ %.0fs + %ds buffer)\n', ...
        i, fileBounds(i), fileBounds(i+1), i, ...
        t_evoked_start_all(i), fileBounds(i+1), t_P9_all(i), cfg.motor_buffer_s);
end

[peaks, good, rate_all, nN] = fn_qualityFilterNeurons(peaks, cfg.fs, cfg.min_rate, cfg.max_rate); %#ok<ASGLU>
fn_plotRateDistribution(rate_all, cfg.recording_ID, nNeurons, cfg.FIGURES_DIR);

fn_plotRaster_sensitization(peaks, cfg.fs, nN, nFrames, t_P9_all, fileBounds, cfg.recording_ID, protocol_label, cfg.FIGURES_DIR);

sigma_s = fn_estimateKernelWidth(peaks, cfg.fs);
spike_conv = fn_gaussConvSpikes(peaks, sigma_s, cfg.fs);
clear peaks;
fprintf('  spike_conv: %.0f MB (single)\n', nN*nFrames*4/1e6);

fn_plotSmoothedActivity_sensitization(spike_conv, cfg.fs, nN, nFrames, sigma_s, t_P9_all, fileBounds, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 2 - PCA + EPOCH CENTROIDS
%  ========================================================================
fprintf('\nSection 2: PCA...\n');

% Same convention as Run_Evoked_Analysis.m: PCA basis is fit from the
% start of the (first) evoked period to the end, not on the baseline.
f_pca_start = max(round(t_evoked_start_all(1) * cfg.fs), 1);
[V, var_exp, mu_pca] = fn_computePCA(spike_conv, f_pca_start, nFrames, cfg.pca_block_size);
fprintf('  PC1=%.1f%%  PC2=%.1f%%  PC3=%.1f%%\n', var_exp(1), var_exp(2), var_exp(3));
fn_plotScree(var_exp, nN, cfg.FIGURES_DIR);

[epoch_scores, centroids] = fn_epochScoresAndCentroids(spike_conv, wins_s, cfg.fs, V, mu_pca, win_labels); %#ok<ASGLU>

fn_plotTrajectory3D_sensitization(epoch_scores, baseline_idx, evoked_idx, var_exp, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 3 - WINDOWED PCA + PARTICIPATION RATIO (per epoch)
%  ========================================================================
fprintf('\nSection 3: Windowed PCA + participation ratio...\n');

win_frames = max(min(round(wins_s * cfg.fs), nFrames), 1);
WPCA = fn_windowedPCA(spike_conv, win_frames, []);
for w = 1:numel(WPCA)
    fprintf('  [%s]  top PC = %.1f%%  |  dims@95%% = %d\n', ...
        win_labels{w}, WPCA(w).var_explained(1), WPCA(w).n_dims_95pct);
end

PR_raw  = arrayfun(@(w) fn_participationRatio(WPCA(w).eigenvalues), 1:numel(WPCA));
PR_norm = PR_raw / nN;
fprintf('  PR/N per epoch: '); fprintf('%.3f  ', PR_norm); fprintf('\n');

fn_plotPRvsStimulus(PR_norm, baseline_idx, evoked_idx, cfg.recording_ID, cfg.FIGURES_DIR);

[pr_t, pr_v] = fn_slidingParticipationRatio(spike_conv, win_f, step_f, cfg.fs, nN);
fn_plotPRTimeseries_sensitization(pr_t, pr_v, t_P9_all, fileBounds, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 4 - SUBSPACE ALIGNMENT ACROSS ALL EPOCHS
%  ========================================================================
fprintf('\nSection 4: Subspace alignment...\n');

% Dimensionality/chance level chosen from the FIRST evoked epoch.
[nDims, chance_lvl] = fn_chooseAlignmentK(cfg.K_ALIGN, WPCA(evoked_idx(1)), nN, cfg.align_var_thresh_pct);

align_mat = fn_epochAlignmentMatrix(WPCA, nDims);
fn_plotAlignmentMatrix(align_mat, win_labels, cfg.recording_ID, cfg.FIGURES_DIR);

align_to_first        = nan(1, nStim);
align_consecutive     = nan(1, max(nStim-1, 0));
align_evoked_to_base1 = nan(1, nStim);

for i = 1:nStim
    align_to_first(i)        = align_mat(evoked_idx(1), evoked_idx(i));
    align_evoked_to_base1(i) = align_mat(baseline_idx,  evoked_idx(i));
end
for i = 1:nStim-1
    align_consecutive(i)     = align_mat(evoked_idx(i), evoked_idx(i+1));
end

fprintf('  Chance level = %.3f\n', chance_lvl);
fprintf('  Alignment Evoked_1 -> Evoked_i: ');     fprintf('%.3f  ', align_to_first); fprintf('\n');
fprintf('  Alignment Evoked_i -> Evoked_i+1: ');   fprintf('%.3f  ', align_consecutive); fprintf('\n');
fprintf('  Alignment Baseline -> Evoked_i: ');   fprintf('%.3f  ', align_evoked_to_base1); fprintf('\n');

fn_plotAlignmentToFirst(align_to_first, chance_lvl, cfg.recording_ID, cfg.FIGURES_DIR);
fn_plotAlignmentConsecutive(align_consecutive, chance_lvl, cfg.recording_ID, cfg.FIGURES_DIR);
fn_plotAlignmentToBaseline1(align_evoked_to_base1, chance_lvl, cfg.recording_ID, cfg.FIGURES_DIR);

% Sliding alignment to the FIRST evoked subspace
evoked_axes_1 = WPCA(evoked_idx(1)).eigenvectors(:,1:nDims);
[al_t, al_v] = fn_slidingSubspaceAlignment(spike_conv, evoked_axes_1, win_f, step_f, cfg.fs, nDims);
fn_plotAlignmentTimeseries_sensitization(al_t, al_v, chance_lvl, t_P9_all, fileBounds, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5a - RECURRENCE ANALYSIS PER STIMULUS
%  ========================================================================
fprintf('\nSection 5a: Recurrence density analysis...\n');

mu_global = mean(spike_conv, 2);

% Same epsilon convention as Run_Evoked_Analysis.m (cfg.MAX_FULL_PTS / cfg.EPS_PCTILE),
% calibrated on Evoked_1.
[traj_ev1_global, ~] = fn_getEpochTrajectory(spike_conv, round(t_evoked_start_all(1)*cfg.fs), round(fileBounds(2)*cfg.fs), ...
    cfg.fs, V, mu_global, cfg.MAX_FULL_PTS);
epsilon_rr = prctile(pdist(traj_ev1_global,'euclidean'), cfg.EPS_PCTILE);
fprintf('  eps = %.4f  (%gth pct of Evoked_1 pairwise distances)\n', epsilon_rr, cfg.EPS_PCTILE);

traj_evoked = cell(1, nStim);
t_evoked_ax = cell(1, nStim);
R_evoked    = cell(1, nStim);
ev_recur_density = nan(1, nStim);

% Baseline recurrence (computed once)
[traj_base1, t_base1_ax] = fn_getEpochTrajectory(spike_conv, max(round(fileBounds(1)*cfg.fs),1), round(t_P9_all(1)*cfg.fs), ...
    cfg.fs, V, mu_global, cfg.MAX_WIN_PTS);
buf_base1 = min(cfg.slide_win_s, 0.25*(t_base1_ax(end)-t_base1_ax(1)));
[base_recur_density_1, ~, ~, ~] = fn_recurrenceDensity(traj_base1, t_base1_ax, epsilon_rr, cfg.MIN_LAG_S, buf_base1);
fprintf('  Baseline RR=%.3f\n', base_recur_density_1);

for i = 1:nStim
    [traj_ev, t_ev_ax] = fn_getEpochTrajectory(spike_conv, round(t_evoked_start_all(i)*cfg.fs), round(fileBounds(i+1)*cfg.fs), ...
        cfg.fs, V, mu_global, cfg.MAX_WIN_PTS);
    buf_ev = min(cfg.slide_win_s, 0.25*(t_ev_ax(end)-t_ev_ax(1)));
    [ev_recur_density(i), ~, ~, ~] = fn_recurrenceDensity(traj_ev, t_ev_ax, epsilon_rr, cfg.MIN_LAG_S, buf_ev);
    [R_ev, ~, ~] = fn_recurrencePlot(traj_ev, epsilon_rr, 'euc');

    traj_evoked{i} = traj_ev; t_evoked_ax{i} = t_ev_ax; R_evoked{i} = R_ev;
    fprintf('  Stim %d: Evoked RR=%.3f\n', i, ev_recur_density(i));
end

cross_recur_to_first = nan(1, nStim);
for i = 1:nStim
    if i == 1
        cross_recur_to_first(i) = ev_recur_density(1);
    else
        [cross_recur_to_first(i), ~, ~] = fn_crossRecurrenceDensity(traj_evoked{1}, traj_evoked{i}, epsilon_rr);
    end
    fprintf('  Stim %d: cross-recurrence to Evoked_1 = %.3f\n', i, cross_recur_to_first(i));
end

fn_plotRecurrenceVsStimulus(ev_recur_density, base_recur_density_1, cross_recur_to_first, cfg.recording_ID, cfg.FIGURES_DIR);

fn_plotRecurrenceSummary_sensitization(traj_evoked, t_evoked_ax, R_evoked, ev_recur_density, epsilon_rr, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5b - ATTRACTOR ONSET DETECTION
%  ========================================================================
fprintf('\nSection 5b: Attractor onset detection...\n');

t_attractor_onset_all = nan(1, nStim);
onset_detected_all    = false(1, nStim);
onset_latency_all     = nan(1, nStim);
onset_all             = cell(1, nStim);

for i = 1:nStim
    onset_i = fn_detectAttractorEpoch(spike_conv, cfg.fs, V, mu_global, t_P9_all(i), fileBounds(i+1), ...
        epsilon_rr, cfg.MIN_LAG_S, cfg.slide_win_s, win_f, step_f, cfg.ONSET_RATIO, cfg.MAX_FULL_PTS, t_evoked_start_all(i));
    onset_all{i} = onset_i;
    t_attractor_onset_all(i) = onset_i.t_lock;
    onset_detected_all(i)    = onset_i.detected;
    onset_latency_all(i)     = onset_i.t_lock - t_P9_all(i);
    if onset_i.detected
        fprintf('  Stim %d: ONSET t=%.1fs (%.1fs after P9)\n', i, onset_i.t_lock, onset_latency_all(i));
    else
        fprintf('  Stim %d: ONSET not detected -- fallback %.1fs\n', i, onset_i.t_lock);
    end
end

fn_plotOnsetLatencyVsStimulus(onset_latency_all, onset_detected_all, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SECTION 5c - PR & ALIGNMENT ON THE ATTRACTOR-DEFINED EVOKED WINDOWS
%  ========================================================================
fprintf('\nSection 5c: PR and alignment on attractor-defined evoked windows...\n');

wins_s_rr     = zeros(1+nStim, 2);
win_labels_rr = cell(1, 1+nStim);
wins_s_rr(1,:)   = [fileBounds(1), t_P9_all(1)];
win_labels_rr{1} = 'Baseline';

for i = 1:nStim
    wins_s_rr(i+1,  :)   = [t_attractor_onset_all(i), fileBounds(i+1)];
    win_labels_rr{i+1}   = sprintf('Attractor-evoked_%d', i);
end
win_frames_rr = max(min(round(wins_s_rr * cfg.fs), nFrames), 1);
WPCA_rr       = fn_windowedPCA(spike_conv, win_frames_rr, []);

PR_raw_rr  = arrayfun(@(w) fn_participationRatio(WPCA_rr(w).eigenvalues), 1:numel(WPCA_rr));
PR_norm_rr = PR_raw_rr / nN;
fprintf('  PR/N [RR-defined] per epoch: '); fprintf('%.3f  ', PR_norm_rr); fprintf('\n');

align_mat_rr = fn_epochAlignmentMatrix(WPCA_rr, nDims);
align_B1_E_rr = nan(1, nStim);
for i = 1:nStim
    align_B1_E_rr(i) = align_mat_rr(baseline_idx, evoked_idx(i));
end
fprintf('  Alignment Baseline -> Attractor-evoked_i: '); fprintf('%.3f  ', align_B1_E_rr); fprintf('\n');
fn_plotAttractorMetrics(PR_norm_rr, align_B1_E_rr, chance_lvl, baseline_idx, evoked_idx, cfg.recording_ID, cfg.FIGURES_DIR);

%% =========================================================================
%  SAVE
%  ========================================================================
save_path   = fullfile(cfg.RESULTS_DIR, sprintf('%s_results_sensitization.mat', cfg.recording_ID));
nN_saved    = nN;
nDims_align = nDims;
recording_ID = cfg.recording_ID; fs = cfg.fs; %#ok<NASGU>
ONSET_RATIO = cfg.ONSET_RATIO; %#ok<NASGU>

save(save_path, ...
    'recording_ID','fs','nStim','fileBounds','t_P9_all','t_evoked_start_all','t_end', ...
    'wins_s','win_labels','baseline_idx','evoked_idx','nN','nN_saved', ...
    'WPCA','PR_raw','PR_norm','pr_t','pr_v', ...
    'align_mat','nDims','nDims_align','chance_lvl', ...
    'align_to_first','align_evoked_to_base1','align_consecutive','al_t','al_v', ...
    'epsilon_rr','ev_recur_density','base_recur_density_1','cross_recur_to_first', ...
    't_attractor_onset_all','onset_detected_all','onset_latency_all','onset_all', ...
    'wins_s_rr','win_labels_rr','WPCA_rr','PR_raw_rr','PR_norm_rr','align_mat_rr','align_B1_E_rr', ...
    'ONSET_RATIO', ...
    '-v7.3');

fprintf('\nSaved: %s\nDone.\n', save_path);
diary off;


%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function c = fn_colEvoked(), c = [0.8 0.1 0.1]; end   % evoked data (same red as Run_Evoked_Analysis)
function c = fn_colBlue(),   c = [0.1 0.3 0.7]; end   % alignment / recurrence / onset plots
function c = fn_colBase(),   c = [0.5 0.5 0.5]; end   % baseline / chance reference (gray)

function fn_labelEvoked(x, y, c)
    if nargin < 3, c = fn_colEvoked(); end
    % In-figure label at the right end of the evoked line (mirrors the Baseline yline label)
    text(x + 0.06, y, 'Evoked', 'Color', c, 'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'middle');
end

function matches = fn_findShallow(rootDir, filename)
    matches = {};
    if isempty(rootDir) || ~isfolder(rootDir), return; end
    candidates = {rootDir};
    d = dir(rootDir); d = d([d.isdir]); d = d(~ismember({d.name}, {'.','..'}));
    for k = 1:numel(d), candidates{end+1} = fullfile(rootDir, d(k).name); end %#ok<AGROW>
    for c = 1:numel(candidates)
        listing = dir(candidates{c}); listing = listing(~[listing.isdir]);
        hit = strcmpi({listing.name}, filename);
        if any(hit)
            names = {listing(hit).name};
            for j = 1:numel(names)
                matches{end+1} = fullfile(candidates{c}, names{j}); %#ok<AGROW>
            end
        end
    end
end

function names = fn_listMatFilesShallow(rootDir)
    names = {};
    if isempty(rootDir) || ~isfolder(rootDir), return; end
    candidates = {rootDir};
    d = dir(rootDir); d = d([d.isdir]); d = d(~ismember({d.name}, {'.','..'}));
    for k = 1:numel(d), candidates{end+1} = fullfile(rootDir, d(k).name); end %#ok<AGROW>
    for c = 1:numel(candidates)
        listing = dir(fullfile(candidates{c}, '*.mat'));
        for j = 1:numel(listing)
            names{end+1} = fullfile(listing(j).folder, listing(j).name); %#ok<AGROW>
        end
    end
end

function fn_plotRaster_sensitization(peaks, fs, nN, nFrames, t_P9_all, fileBounds, recording_ID, protocol_label, figuresDir)
    t_axis = (0:nFrames-1) / fs;
    fig = figure('Name', 'Raster (sensitization)', 'Position', [100 100 1400 700], 'Visible', 'off');
    hold on;
    [rr, cc] = find(peaks);
    plot(t_axis(cc), rr, 'k.', 'MarkerSize', 3);
    hP9 = xline(t_P9_all(1), 'r-', 'LineWidth', 1.2);
    for k = 2:numel(t_P9_all)
        xline(t_P9_all(k), 'r-', 'LineWidth', 1.2);
    end
    hBound = xline(fileBounds(2), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    for k = 3:numel(fileBounds)-1
        xline(fileBounds(k), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
    xlim([0, fileBounds(end)]);
    xticks(0:250:fileBounds(end));
    ylim([0.5, nN+0.5]);
    xlabel('Time (s)'); ylabel('Neuron #');
    title(sprintf('Raster | %s | %s', recording_ID, protocol_label), 'Interpreter', 'none');
    legend([hP9, hBound], {'P9', 'File boundary'}, 'Location', 'northeastoutside');
    box on;
    exportgraphics(fig, fullfile(figuresDir, '01b_raster_sensitization.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotSmoothedActivity_sensitization(spike_conv, fs, nN, nFrames, sigma_s, t_P9_all, fileBounds, recording_ID, figuresDir)
    t_axis = (0:nFrames-1) / fs;
    fig = figure('Name', 'Smoothed population activity (sensitization)', 'Position', [100 100 1400 500], 'Visible', 'off');
    imagesc(t_axis([1 end]), [1 nN], spike_conv);
    colormap(hot); colorbar; hold on;
    for k = 1:numel(t_P9_all)
        xline(t_P9_all(k), 'w-', 'P9', 'LineWidth', 1.2, 'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'top');
    end
    for k = 2:numel(fileBounds)-1
        xline(fileBounds(k), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
    xlim([0, fileBounds(end)]); xticks(0:250:fileBounds(end)); ylim([0.5, nN+0.5]);
    xlabel('Time (s)'); ylabel('Neuron #');
    title(sprintf('Gaussian-smoothed activity (sigma=%.1f ms) | %s', sigma_s*1000, recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '01c_smoothed_sensitization.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotTrajectory3D_sensitization(epoch_scores, baseline_idx, evoked_idx, var_exp, recording_ID, figuresDir)
    nStim = numel(evoked_idx);

    cmap = [0.85 0.10 0.10;   % Evoked_1  red
            0.10 0.35 0.90;   % Evoked_2  blue
            0.10 0.65 0.20;   % Evoked_3  green
            0.95 0.55 0.00;   % Evoked_4  orange
            0.55 0.15 0.70];  % Evoked_5  purple
    if nStim > size(cmap, 1), cmap = lines(nStim); end
    markers = {'o', 's', 'd', '^', 'v'};

    fig = figure('Name', 'Population trajectory', 'Position', [200 100 950 850], 'Visible', 'off');
    hold on;

    % Baseline path (gray)
    b = epoch_scores{baseline_idx};
    hB = plot3(b(:,1), b(:,2), b(:,3), '-', 'Color', fn_colBase(), 'LineWidth', 0.75);

    % Evoked paths, one color per stimulus, marker at the start of each
    hEv = gobjects(1, nStim);
    for i = 1:nStim
        e = epoch_scores{evoked_idx(i)};
        hEv(i) = plot3(e(:,1), e(:,2), e(:,3), '-', 'Color', cmap(i,:), 'LineWidth', 1.5);
        m_shape = markers{mod(i-1, numel(markers))+1};
        plot3(e(1,1), e(1,2), e(1,3), m_shape, 'MarkerSize', 9, ...
              'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
    end

    xlabel(sprintf('PC1 (%.1f%%)', var_exp(1)));
    ylabel(sprintf('PC2 (%.1f%%)', var_exp(2)));
    zlabel(sprintf('PC3 (%.1f%%)', var_exp(3)));
    grid on; view([-35 25]);
    title(sprintf('Population trajectory across sensitization | %s', recording_ID), 'Interpreter', 'none');
    legend([hB, hEv], [{'Baseline'}, arrayfun(@(i) sprintf('Evoked_%d', i), 1:nStim, 'UniformOutput', false)], ...
        'Location', 'eastoutside');

    exportgraphics(fig, fullfile(figuresDir, '02_trajectory_3D_sensitization.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotPRvsStimulus(PR_norm, baseline_idx, evoked_idx, recording_ID, figuresDir)
    nStim = numel(evoked_idx);
    fig = figure('Name', 'PR vs stimulus', 'Position', [200 100 700 450], 'Visible', 'off');
    plot(1:nStim, PR_norm(evoked_idx), 's-', 'Color', fn_colEvoked(), 'LineWidth', 1.5); hold on;
    yline(PR_norm(baseline_idx), '--', 'Baseline', 'Color', fn_colBase(), 'LineWidth', 1.3, 'LabelHorizontalAlignment', 'left');
    fn_labelEvoked(nStim, PR_norm(evoked_idx(end)));
    xlim([0.8, nStim + 0.9]); xticks(1:nStim); grid on;
    xlabel('Stimulus #'); ylabel('PR / N');
    title(sprintf('Participation ratio across repeated P9 stimuli | %s', recording_ID), 'Interpreter', 'none');
    legend off;
    exportgraphics(fig, fullfile(figuresDir, '03b_PR_vs_stimulus.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotPRTimeseries_sensitization(pr_t, pr_v, t_P9_all, fileBounds, recording_ID, figuresDir)
    fig = figure('Name', 'PR timeseries (sensitization)', 'Position', [100 100 1400 450], 'Visible', 'off');
    plot(pr_t, pr_v, 'k-', 'LineWidth', 1.0); hold on;
    for k = 1:numel(t_P9_all)
        xline(t_P9_all(k), 'r-', 'LineWidth', 1.2);
    end
    for k = 2:numel(fileBounds)-1
        xline(fileBounds(k), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
    xlim([0, fileBounds(end)]); xticks(0:250:fileBounds(end));
    xlabel('Time (s)'); ylabel('Participation ratio (PR/N)');
    title(sprintf('Sliding participation ratio | %s', recording_ID), 'Interpreter', 'none');
    box on;
    exportgraphics(fig, fullfile(figuresDir, '03_PR_timeseries_sensitization.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotAlignmentToFirst(align_to_first, chance_lvl, recording_ID, figuresDir)
    nStim = numel(align_to_first);
    fig = figure('Name', 'Alignment to Evoked_1', 'Position', [200 100 700 450], 'Visible', 'off');
    plot(1:nStim, align_to_first, 'd-', 'Color', fn_colBlue(), 'LineWidth', 1.5); hold on;
    yline(chance_lvl, ':', 'Chance', 'Color', fn_colBase(), 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
    xlabel('Stimulus #'); ylabel('Subspace alignment to Evoked_1');
    xticks(1:nStim); xlim([0.8, nStim+0.2]); grid on;
    title(sprintf('Evoked-subspace drift from first stimulus | %s', recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '04d_alignment_to_first.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotAlignmentToBaseline1(align_evoked_to_base1, chance_lvl, recording_ID, figuresDir)
    nStim = numel(align_evoked_to_base1);
    fig = figure('Name', 'Alignment to Baseline', 'Position', [200 100 700 450], 'Visible', 'off');
    plot(1:nStim, align_evoked_to_base1, 'd-', 'Color', fn_colBlue(), 'LineWidth', 1.5); hold on;
    yline(chance_lvl, ':', 'Chance', 'Color', fn_colBase(), 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
    xlabel('Stimulus #'); ylabel('Subspace alignment (Baseline \rightarrow Evoked_i)');
    xticks(1:nStim); xlim([0.8, nStim+0.2]); grid on;
    title(sprintf('Evoked alignment relative to initial Baseline | %s', recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '04e_alignment_to_baseline1.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotAlignmentConsecutive(align_consecutive, chance_lvl, recording_ID, figuresDir)
    % One point per stimulus PAIR, placed exactly on an integer tick and
    % labelled "1->2", "2->3", ... so points no longer float between ticks.
    nPairs = numel(align_consecutive);
    fig = figure('Name', 'Consecutive Evoked alignment', 'Position', [200 100 700 450], 'Visible', 'off');
    plot(1:nPairs, align_consecutive, 'd-', 'Color', fn_colBlue(), 'LineWidth', 1.5); hold on;
    yline(chance_lvl, ':', 'Chance', 'Color', fn_colBase(), 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
    pairLabels = arrayfun(@(k) sprintf('%d->%d', k, k+1), 1:nPairs, 'UniformOutput', false);
    xticks(1:nPairs); xticklabels(pairLabels);
    xlim([0.5, nPairs+0.5]); grid on;
    xlabel('Stimulus pair (i -> i+1)'); ylabel('Subspace alignment (Evoked_i \rightarrow Evoked_{i+1})');
    title(sprintf('Step-to-step Evoked alignment | %s', recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '04f_alignment_consecutive.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotAlignmentTimeseries_sensitization(al_t, al_v, chance_lvl, t_P9_all, fileBounds, recording_ID, figuresDir)
    fig = figure('Name', 'Alignment timeseries', 'Position', [100 100 1400 450], 'Visible', 'off');
    plot(al_t, al_v, 'k-', 'LineWidth', 1.0); hold on;
    yline(chance_lvl, ':', 'Color', fn_colBase(), 'LineWidth', 1.0);
    for k = 1:numel(t_P9_all)
        xline(t_P9_all(k), 'r-', 'LineWidth', 1.2);
    end
    for k = 2:numel(fileBounds)-1
        xline(fileBounds(k), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.0);
    end
    xlim([0, fileBounds(end)]); xticks(0:250:fileBounds(end));
    xlabel('Time (s)'); ylabel('Alignment to Evoked_1 subspace');
    title(sprintf('Sliding subspace alignment | %s', recording_ID), 'Interpreter', 'none');
    box on;
    exportgraphics(fig, fullfile(figuresDir, '04c_alignment_timeseries_sensitization.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotRecurrenceVsStimulus(ev_rr, base_rr_1, cross_to_first, recording_ID, figuresDir)
    nStim = numel(ev_rr);
    fig = figure('Name', 'Recurrence vs stimulus', 'Position', [200 100 1100 450], 'Visible', 'off');

    subplot(1,2,1);
    plot(1:nStim, ev_rr, 's-', 'Color', fn_colBlue(), 'LineWidth', 1.4); hold on;
    yline(base_rr_1, '--', 'Baseline', 'Color', fn_colBase(), 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
    fn_labelEvoked(nStim, ev_rr(end), fn_colBlue());
    xlim([0.8, nStim + 0.9]);
    xlabel('Stimulus #'); ylabel('Self-recurrence density');
    xticks(1:nStim); grid on;
    legend off;   % legend removed: baseline is already labelled on the plot
    title('Within-epoch stereotypy');

    subplot(1,2,2);
    plot(1:nStim, cross_to_first, 'd-', 'Color', fn_colBlue(), 'LineWidth', 1.6);
    xlabel('Stimulus #'); ylabel('Cross-recurrence to Evoked_1');
    xticks(1:nStim); grid on;
    title('Drift from the first evoked response');

    sgtitle(sprintf('Recurrence density across repeated P9 stimuli | %s', recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '05b_recurrence_vs_stimulus.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotOnsetLatencyVsStimulus(onset_latency_all, onset_detected_all, recording_ID, figuresDir)
    nStim = numel(onset_latency_all);
    fig = figure('Name', 'Onset latency vs stimulus', 'Position', [200 100 700 450], 'Visible', 'off');
    plot(1:nStim, onset_latency_all, 'o-', 'Color', fn_colBlue(), 'LineWidth', 1.4); hold on;
    notDet = find(~onset_detected_all);
    if ~isempty(notDet)
        plot(notDet, onset_latency_all(notDet), 'kx', 'MarkerSize', 10, 'LineWidth', 1.5);
    end
    xlabel('Stimulus #'); ylabel('Attractor onset latency (s, from P9)');
    xticks(1:nStim); grid on;
    title(sprintf('Onset latency across repeated P9 stimuli | %s', recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '05e_onset_latency_vs_stimulus.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotRecurrenceSummary_sensitization(traj_evoked, t_evoked_ax, R_evoked, ev_recur_density, epsilon_rr, recording_ID, figuresDir)
    nStim = numel(traj_evoked);
    fig = figure('Name', 'Recurrence plots per stimulus', 'Position', [100 100 260*nStim 320], 'Visible', 'off');
    for i = 1:nStim
        subplot(1, nStim, i);
        imagesc(t_evoked_ax{i}, t_evoked_ax{i}, R_evoked{i});
        axis square; set(gca, 'YDir', 'normal');
        colormap(gca, flipud(gray));
        title(sprintf('Evoked_%d\nRR=%.2f', i, ev_recur_density(i)));
        if i == 1, xlabel('Time (s)'); ylabel('Time (s)'); end
    end
    sgtitle(sprintf('Recurrence plots (eps=%.4f) | %s', epsilon_rr, recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '05_recurrence_all_sensitization.png'), 'Resolution', 500);
    close(fig);
end

function fn_plotAttractorMetrics(PR_norm_rr, align_B1_E_rr, chance_lvl, baseline_idx, evoked_idx, recording_ID, figuresDir)
    nStim = numel(align_B1_E_rr);
    fig = figure('Name', 'PR and Alignment of RR-defined Epochs', 'Position', [200 100 1100 450], 'Visible', 'off');

    subplot(1,2,1);
    plot(1:nStim, PR_norm_rr(evoked_idx), 's-', 'Color', fn_colEvoked(), 'LineWidth', 1.5); hold on;
    yline(PR_norm_rr(baseline_idx), '--', 'Baseline', 'Color', fn_colBase(), 'LineWidth', 1.3, 'LabelHorizontalAlignment', 'left');
    fn_labelEvoked(nStim, PR_norm_rr(evoked_idx(end)));
    xlim([0.8, nStim + 0.9]);
    xlabel('Stimulus #'); ylabel('PR / N');
    xticks(1:nStim); grid on;
    title('Participation ratio');

    subplot(1,2,2);
    plot(1:nStim, align_B1_E_rr, 'd-', 'Color', fn_colBlue(), 'LineWidth', 1.5); hold on;
    yline(chance_lvl, ':', 'Chance', 'Color', fn_colBase(), 'LineWidth', 1.0, 'LabelHorizontalAlignment', 'left');
    xlabel('Stimulus #'); ylabel('Alignment (Baseline \rightarrow Evoked_i)');
    xticks(1:nStim); grid on;
    title('Subspace alignment');

    sgtitle(sprintf('PR and Subspace alignment on RR-defined epochs | %s', recording_ID), 'Interpreter', 'none');
    exportgraphics(fig, fullfile(figuresDir, '05d_attractor_metrics.png'), 'Resolution', 500);
    close(fig);
end
