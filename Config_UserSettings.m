%% CONFIG_USERSETTINGS
%  Single, user-edited source of all parameters for the Aplysia attractor
%  perturbation analysis pipeline. No parameter should be hard-coded
%  anywhere else in the toolbox - if you need to change something, it
%  should be changeable from here.
%
%  Usage: edit the values below, save, then run Run_Attractor_Analysis.m
%  (which calls this script automatically).

cfg = struct();

% ---- Paths --------------------------------------------------------------
cfg.DATA_FILE   = fullfile(toolboxRoot, 'YourDataFile.mat');
cfg.RESULTS_DIR = fullfile(toolboxRoot, 'Results');
cfg.FIGURES_DIR = fullfile(toolboxRoot, 'Figures');

% ---- Recording metadata --------------------------------------------------
cfg.recording_ID = 'RecID';   % like Jan01/Feb25 etc. 
cfg.protocol     = '12min';   % '12min' or '20min'
cfg.fs           = 1629;      % sampling rate (fps)

% ---- Data loading ---------------------------------------------------------
cfg.chunk_size = 50;          % neurons loaded per chunk (memory control)

% ---- Neuron quality filter ------------------------------------------------
cfg.min_rate = 0.01;          % Hz, minimum mean event rate to keep a neuron
cfg.max_rate = 50;            % Hz, maximum mean event rate to keep a neuron

% ---- Epoch timing --------------------------------------------------------
cfg.motor_buffer_s   = 30;    % s, buffer after P9 before calling it "Evoked"
cfg.recovery_delay_s = 60;    % s, delay after C2 before calling it "Recovery"

% ---- Sliding window (PR, alignment, recurrence time series) --------------
cfg.slide_win_s  = 60;        % s
cfg.slide_step_s =  5;        % s

% ---- PCA / subspace alignment --------------------------------------------
cfg.K_ALIGN        = 'auto';  % 'auto' = dims at cfg.align_var_thresh_pct evoked variance; or integer
cfg.align_var_thresh_pct = 80;% percent evoked variance for 'auto' K choice
cfg.pca_block_size = 5000;    % columns per block for covariance accumulation

% ---- Recurrence-density analysis -----------------------------------------
cfg.MAX_WIN_PTS  = 500;       % downsample target, per-epoch trajectories (plots)
cfg.MAX_FULL_PTS = 4000;      % downsample target, full-recording pass
cfg.MIN_LAG_S    = 2.0;       % Theiler-window-style minimum lag
cfg.ONSET_RATIO  = 0.90;      % recurrence-density criterion for "locked on"

% ---- Principled RR threshold (documentation only; not currently used by
%      the epoch-detection logic, which is percentile/ratio based - kept
%      here so N_SIGMA/N_CONSEC stay user-editable if that logic is
%      reinstated) --------------------------------------------------------
%   threshold = baseline_RR_mean + N_SIGMA x baseline_RR_std
%   N_SIGMA = 2   -> ~97.7% confidence above the baseline null (one-sided)
%   N_SIGMA = 1.5 -> more sensitive (pick up earlier attractor lock-in)
%   N_SIGMA = 3   -> conservative (only count very strong attractors)
cfg.N_SIGMA  = 2;
cfg.N_CONSEC = 1;

% ---- FFT / dominant oscillation period ------------------------------------
cfg.fft_cutoff_period_s = 250;  % s, high-pass cutoff (suppress drift slower than this)

% ---- Synthetic validation (only used by the Validation/ scripts) --------
cfg.validation_seed_alignment = 1;
cfg.validation_seed_recurrence = 42;
cfg.validation_nReps = 50;
cfg.validation_tol   = 0.05;
