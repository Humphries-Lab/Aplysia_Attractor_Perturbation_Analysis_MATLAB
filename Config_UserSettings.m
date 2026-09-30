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
% cfg.DATA_FILE_INPUT accepts either:
%   - a bare filename, e.g. 'Sep1225.mat'  -> resolved relative to toolboxRoot
%   - a full path anywhere on disk, e.g. 'C:\Data\Sep1225.mat', or
%     '/home/user/data/Sep1225.mat' -> used as-is
cfg.DATA_FILE_INPUT = 'data.mat';   % <-- EDIT THIS
[pathPart, ~, ~] = fileparts(cfg.DATA_FILE_INPUT);
if isempty(pathPart)
    cfg.DATA_FILE = fullfile(toolboxRoot, cfg.DATA_FILE_INPUT);
else
    cfg.DATA_FILE = cfg.DATA_FILE_INPUT;
end

cfg.RESULTS_DIR = fullfile(toolboxRoot, 'Results');
cfg.FIGURES_DIR = fullfile(toolboxRoot, 'Figures');
cfg.LOGS_DIR = fullfile(toolboxRoot, 'Logs');

% ---- Recording metadata --------------------------------------------------
cfg.recording_ID = 'Sep12';   % EDIT AS REQUIRED 
cfg.protocol     = '12min';   % EDIT AS REQUIRED --> '12min' or '20min' etc. [or '25min' (for concatenated files])
cfg.fs           = 1629;      % sampling rate (fps)

% ---- Data loading ---------------------------------------------------------
cfg.chunk_size = 50;          % neurons loaded per chunk to save memory

% ---- Neuron quality filter ------------------------------------------------
cfg.min_rate = 0.01;          % Hz, minimum mean event rate to keep a neuron
cfg.max_rate = 50;            % Hz, maximum mean event rate to keep a neuron

% ---- Epoch timing --------------------------------------------------------
cfg.t_P9_evoked  = 120;   % s, stimulus onset --> EDIT AS REQUIRED
cfg.t_end_evoked = 420;   % s, end of recording/analysis window --> EDIT AS REQUIRED
cfg.motor_buffer_s   = 30;    % s, buffer after P9 before calling it "Evoked"
cfg.recovery_delay_s = 60;    % s, delay after C2 before calling it "Recovery"

% ---- Sliding window (PR, alignment, recurrence time series) --------------
cfg.slide_win_s  = 60;        % s
cfg.slide_step_s =  5;        % s

% ---- PCA / subspace alignment --------------------------------------------
cfg.K_ALIGN        = 'auto';  % 'auto' = dims at cfg.align_var_thresh_pct evoked variance; or integer
cfg.align_var_thresh_pct = 80; % percent evoked variance for 'auto' K choice
cfg.pca_block_size = 5000;    % columns per block for covariance accumulation

% ---- Recurrence-density analysis -----------------------------------------
cfg.MAX_WIN_PTS  = 500;       % downsample target, per-epoch trajectories (plots)
cfg.MAX_FULL_PTS = 4000;      % downsample target, full-recording pass
cfg.MIN_LAG_S    = 2.0;       % Theiler-window-style minimum lag
cfg.ONSET_RATIO  = 0.90;      % recurrence-density criterion for "locked on"
cfg.EPS_PCTILE   = 10;        % percentile of pairwise distances used to calibrate epsilon_rr

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

% ---- Synthetic validation (only used by the Validation/ scripts) --------
cfg.validation_seed_alignment = 1;
cfg.validation_seed_recurrence = 42;
cfg.validation_nReps = 50;
cfg.validation_tol   = 0.05;

% ---- Dynamical epochs (cycle-by-cycle jPCA alignment; Run_Attractor_Analysis Section 5d) ----
%      Uses Andrea Colins Rodriguez's Dynamical_epochs + jPCA_Aplysia code
%      (added as a git submodule in Dynamical_epochs/):
%      git clone --recurse-submodules <this repo URL>
cfg.DYN_ENABLE       = true;   % false = skip Section 5d
cfg.DYN_CODE_DIR     = fullfile(toolboxRoot, 'Dynamical_epochs');
cfg.DYN_WINDOW       = 'full'; % 'full' = whole recording: epochs are found WITHOUT protocol times (they only label the cycles afterwards)
                               % 'post_evoked' = from t_evoked_start to the end (skips the P9 transient; baseline not analysed)
cfg.DYN_REC_SWEEP    = [5 10 15]; % recurrence thresholds (%) re-run to test robustness; [] = skip the sweep
cfg.DYN_BIN_MS       = 50;     % ms, bin size for jPCA / cycle detection
cfg.DYN_MIN_SEG_S    = 10;     % s, minimum cycle duration
cfg.DYN_REC_THRESH   = 10;     % %, recurrence threshold (~5 young, ~10 older animals)
cfg.DYN_SEED         = 1;      % rng seed (clustering is stochastic)
cfg.DYN_NNULL        = 2000;   % shuffle samples for the chance level of cycle-to-cycle alignment
cfg.DYN_NSURR        = 5;      % surrogate runs (each neuron circularly shifted independently); 0 = skip. Each run repeats the jPCA cycle search.
cfg.DYN_NPERM        = 10000;  % permutations for the landmark-vs-change-point test
