function cohort = fn_cohortAnalysis(RESULTS_DIR, recording_IDs, varargin)
%FN_COHORTANALYSIS  Wrapper for cohort-level (across-recording/animal)
%  analysis of outputs produced by Run_Attractor_Analysis.
%
%   cohort = FN_COHORTANALYSIS(RESULTS_DIR, recording_IDs)
%   cohort = FN_COHORTANALYSIS(..., 'Name', Value, ...)
%
%  PURPOSE
%  -------
%  Run_Attractor_Analysis.m is a SINGLE-RECORDING pipeline: it runs on
%  one animal/session at a time and saves one '<recording_ID>_results.mat'
%  per run. This function is the entry point for the step that comes
%  AFTER that -- pooling those per-recording .mat files across animals
%  so the results are ready for cohort-level statistical testing.
%
%  This function is a WRAPPER / SCAFFOLD ONLY. It does not decide or run
%  any statistical test. It only:
%    1. Loads every requested '<recording_ID>_results.mat' from
%       RESULTS_DIR
%    2. Extracts a fixed set of scalar summary metrics per recording
%       into one tidy row-per-animal table
%    3. Returns that table (plus the untouched raw per-recording
%       structs) so the experimentalist can apply whichever
%       statistical test/model fits their design, sample size, and
%       hypotheses -- see "STATISTICAL ANALYSIS" below.
%
%  Which test to run is deliberately NOT hard-coded here: it depends on
%  experimental-design choices (paired vs unpaired structure, number of
%  animals/recordings per animal, whether distributional assumptions
%  hold, planned vs exploratory comparisons, etc.) that only the
%  experimentalist can make.
%
%  INPUTS
%  ------
%  RESULTS_DIR    - path to the folder containing '<ID>_results.mat'
%                    files (same as cfg.RESULTS_DIR in
%                    Config_UserSettings.m)
%  recording_IDs  - cell array of recording_ID strings to pool, e.g.
%                    {'Animal1_Trial1','Animal2_Trial1','Animal3_Trial1'}
%                    Pass {} or omit to auto-discover every
%                    '*_results.mat' file present in RESULTS_DIR.
%
%  NAME-VALUE PARAMETERS
%  ----------------------
%  'Metrics'    - cellstr, which summary fields to pull out of each
%                 result file (default: the list returned by the local
%                 function fn_defaultCohortMetrics(), below). Restrict
%                 this to only the metrics relevant to your question.
%  'SaveTable'  - true/false, also write the pooled table to
%                 <RESULTS_DIR>/CohortAnalysis_Summary.mat (default: true)
%
%  OUTPUT
%  ------
%  cohort.summaryTable - MATLAB table, one row per recording. Columns =
%                        recording_ID, protocol, + requested Metrics.
%  cohort.raw          - 1xN struct array, the untouched contents of
%                        every loaded '*_results.mat' file, indexed in
%                        the same order as summaryTable's rows (for
%                        anything not captured as a scalar summary,
%                        e.g. full pr_t/pr_v or al_t/al_v timeseries).
%  cohort.missing      - cellstr of recording_IDs that were requested
%                        but had no matching '*_results.mat' file, or
%                        whose file failed to load.
%
%  =====================================================================
%  STATISTICAL ANALYSIS -- LEFT TO THE EXPERIMENTALIST
%  =====================================================================
%  cohort.summaryTable is deliberately just tidy data. Below is a menu of
%  what it is set up for, and the parameters that would need to be
%  decided before running any test. No test is chosen or run here.
%
%  Run_Attractor_Analysis computes every metric under TWO epoch
%  definitions, and the default table pulls both:
%    - FIXED epochs      : Baseline / Evoked / Recovery, from protocol
%                           timing alone (t_P9, t_C2, motor/recovery
%                           delays) -- suffix-free column names below.
%    - ATTRACTOR-DEFINED  : Baseline / Attractor-evoked / Attractor-
%      (RR-based) epochs    recovery, from when recurrence density
%                           actually crosses cfg.ONSET_RATIO -- '_rr'
%                           suffix on column names below.
%  Keeping both lets you check whether a result depends on how the
%  epoch boundary was drawn (see item 6).
%
%  1) SUBSPACE ALIGNMENT vs CHANCE
%     (columns: align_ER, align_ER_corrected, chance_lvl  |
%      RR-defined: align_ER_rr)
%     - Question: is evoked->recovery alignment above the chance level
%       expected for random K-dim subspaces, and is that consistent
%       across animals?
%     - Candidate approaches: one-sample t-test or Wilcoxon signed-rank
%       of align_ER_corrected against 0; a permutation/shuffle test
%       using each animal's own chance_lvl as the null; bootstrap CI on
%       the cohort mean of align_ER_corrected; the same tests repeated on
%       align_ER_rr to see if attractor-defined epochs change the
%       conclusion.
%     - Parameters to decide: parametric vs non-parametric (depends on
%       normality and N); one- vs two-tailed; align_ER_corrected already
%       removes each animal's own chance level before pooling -- decide
%       whether that per-animal correction is the right normalization
%       for your question, or whether align_ER and chance_lvl should be
%       modeled jointly instead (e.g. as a ratio, or with chance_lvl as
%       a covariate). chance_lvl is shared by the fixed and RR-defined
%       alignment values (K and N don't change, only the window does).
%
%  2) DIMENSIONALITY (PARTICIPATION RATIO) ACROSS EPOCHS
%     (columns: PR_norm_Baseline, PR_norm_Evoked, PR_norm_Recovery |
%      RR-defined: PR_norm_rr_Baseline, PR_norm_rr_AttractorEvoked,
%      PR_norm_rr_AttractorRecovery)
%     - Question: does normalized participation ratio change across
%       epochs, and does that pattern hold across animals -- under
%       either epoch definition?
%     - Candidate approaches: repeated-measures ANOVA or a linear
%       mixed-effects model (epoch as fixed effect, animal/recording as
%       random intercept) for the 3-epoch comparison; a paired t-test
%       or Wilcoxon signed-rank for a single pairwise contrast (e.g.
%       Evoked vs Recovery only); repeat on the PR_norm_rr_* columns.
%     - Parameters to decide: fixed-effects vs mixed-effects (mixed
%       effects needed if there are multiple recordings per animal);
%       sphericity assumption for RM-ANOVA; multiple-comparison
%       correction across epoch pairs (Tukey, Bonferroni, or FDR).
%
%  3) RECURRENCE DENSITY / RATE
%     (columns: RQA_ev_RR, RQA_re_RR -- fixed-epoch, own-epsilon self-
%      recurrence for Evoked/Recovery; cross_recur_density -- shared-
%      epsilon cross-recurrence, evoked found in recovery)
%     - Question: do evoked / recovery / cross recurrence densities
%       differ from each other, or from a null?
%     - Candidate approaches: paired test (evoked vs recovery, within
%       animal); comparison against a surrogate/shuffled null; effect
%       size (Cohen's d or rank-biserial correlation) alongside any
%       p-value.
%     - Parameters to decide: whether recurrence densities (bounded on
%       [0,1], often skewed) need a variance-stabilizing transform
%       before a parametric test; whether the null should be analytic
%       or resampled per animal. Note there is no separate "RR-defined"
%       recurrence density column by default -- RQA_ev_RR/RQA_re_RR
%       already come from the attractor-locked portion of each epoch
%       when onset_detected/return_detected is true (see
%       Run_Attractor_Analysis.m Section 5). A recurrence density
%       computed on the *fixed*-window-only portion (ignoring the
%       detected lock time) is not saved as a scalar; if you want that
%       comparison, derive it from cohort.raw(i).onset_win_v /
%       .return_win_v (full per-window recurrence-density timeseries,
%       still available per recording) restricted to the fixed-window
%       indices yourself.
%
%  4) ATTRACTOR ONSET / RETURN TIMING
%     (columns: t_attractor_onset, t_attractor_return, onset_detected,
%      return_detected)
%     - Question: is onset/return latency consistent across animals,
%       and how reliably does detection succeed?
%     - Candidate approaches: descriptive statistics (median/IQR) of
%       latency given the typically small N in this kind of recording;
%       report detection rate (mean of onset_detected/return_detected)
%       as a proportion with a binomial CI; only move to inferential
%       statistics on latency itself if N supports it.
%     - Parameters to decide: how to handle non-detections (drop,
%       impute, or treat as right-censored -- censored/survival methods
%       may be more appropriate than dropping non-detections outright).
%
%  5) CROSS-CUTTING CONSIDERATIONS (all of the above)
%     - N animals/recordings available, and whether that N supports the
%       assumptions of a parametric test at all.
%     - Family-wise correction if several of the metrics/epochs above
%       are tested together (Bonferroni, Holm, or Benjamini-Hochberg
%       FDR).
%     - Whether comparisons are planned/confirmatory or exploratory --
%       worth deciding, and ideally recording, before looking at the
%       pooled table.
%
%  6) FIXED vs ATTRACTOR-DEFINED EPOCHS -- AGREEMENT / SENSITIVITY
%     (compare any unsuffixed column to its '_rr' counterpart, e.g.
%      PR_norm_Evoked vs PR_norm_rr_AttractorEvoked, or align_ER vs
%      align_ER_rr)
%     - Question: does the result depend on whether epoch boundaries
%       come from fixed protocol timing or from detected attractor
%       lock-on? If both give the same answer, that is itself evidence
%       the result isn't an artifact of the epoch-timing choice.
%     - Candidate approaches: paired test between the fixed and '_rr'
%       version of the same metric (within animal); a Bland-Altman-
%       style agreement plot (difference vs mean of the two
%       definitions) rather than a significance test, since the goal is
%       agreement, not detecting a difference; correlation between the
%       two across animals as a weaker consistency check.
%     - Parameters to decide: this is usually a robustness/sensitivity
%       analysis rather than a hypothesis test in its own right -- worth
%       deciding up front whether a formal test is even the right tool,
%       versus reporting both epoch definitions alongside each other.
%
%  This function stops at producing cohort.summaryTable. Choosing (and
%  running) the right test from the menu above -- or a different one
%  entirely -- is a scientific judgment call for the experimentalist,
%  not something this pipeline should decide on their behalf.
%
%  EXAMPLE
%  -------
%    ids = {'Animal1_Trial1','Animal2_Trial1','Animal3_Trial1'};
%    cohort = fn_cohortAnalysis(cfg.RESULTS_DIR, ids);
%
%    % cohort.summaryTable is now ready for, e.g.:
%    %   [~,p] = ttest(cohort.summaryTable.align_ER_corrected);
%    % or a linear mixed-effects model / permutation test / etc. --
%    % see the menu above for what to consider before picking one.

% ---- parse optional Name-Value parameters ---------------------------------
p = inputParser;
addParameter(p, 'Metrics', fn_defaultCohortMetrics(), @iscellstr);
addParameter(p, 'SaveTable', true, @islogical);
parse(p, varargin{:});
metrics    = p.Results.Metrics;
saveTable  = p.Results.SaveTable;

% ---- resolve which files to load ------------------------------------------
if nargin < 2 || isempty(recording_IDs)
    matFiles = dir(fullfile(RESULTS_DIR, '*_results.mat'));
    recording_IDs = erase({matFiles.name}, '_results.mat');
end
recording_IDs = recording_IDs(:)';

% ---- load each recording's results and pull out summary metrics -----------
rows    = {};
rawList = {};
missing = {};

for i = 1:numel(recording_IDs)
    thisID = recording_IDs{i};
    fpath  = fullfile(RESULTS_DIR, sprintf('%s_results.mat', thisID));

    if ~isfile(fpath)
        missing{end+1} = thisID; %#ok<AGROW>
        continue;
    end

    try
        S = load(fpath);
    catch
        missing{end+1} = thisID; %#ok<AGROW>
        continue;
    end

    row = table();
    row.recording_ID = string(thisID);
    row.protocol      = string(fn_getFieldOrNaN(S, 'protocol'));

    for m = 1:numel(metrics)
        row.(metrics{m}) = fn_getFieldOrNaN(S, metrics{m});
    end

    rows{end+1}    = row; %#ok<AGROW>
    rawList{end+1} = S;   %#ok<AGROW>
end

if isempty(rows)
    cohort.summaryTable = table();
else
    cohort.summaryTable = vertcat(rows{:});
end
cohort.raw     = [rawList{:}];
cohort.missing = missing;

if ~isempty(missing)
    fprintf('fn_cohortAnalysis: %d recording(s) skipped (no/unreadable results file):\n', numel(missing));
    fprintf('  %s\n', strjoin(missing, ', '));
end

% ---- optionally persist the pooled table -----------------------------------
if saveTable && ~isempty(cohort.summaryTable)
    outPath = fullfile(RESULTS_DIR, 'CohortAnalysis_Summary.mat');
    summaryTable = cohort.summaryTable; %#ok<NASGU>
    save(outPath, 'summaryTable');
    fprintf('fn_cohortAnalysis: pooled %d recording(s) -> %s\n', ...
        height(cohort.summaryTable), outPath);
end

end

% =============================================================================
function metrics = fn_defaultCohortMetrics()
%FN_DEFAULTCOHORTMETRICS  Default scalar summary fields pulled from each
%  '<recording_ID>_results.mat' into the cohort-level summary table.
%  Edit this list (or pass 'Metrics' to fn_cohortAnalysis) to add or
%  drop fields -- anything saved by Run_Attractor_Analysis's SAVE block
%  can be requested here, as long as it is scalar (per-window timeseries
%  like pr_t/pr_v stay in cohort.raw instead).

metrics = { ...
    ... % -- subspace alignment: fixed epochs, then attractor(RR)-defined --
    'align_ER', 'align_ER_corrected', 'chance_lvl', 'nDims_align', ...
    'align_ER_rr', ...
    ... % -- participation ratio: fixed epochs, then attractor(RR)-defined --
    'PR_norm_Baseline', 'PR_norm_Evoked', 'PR_norm_Recovery', ...
    'PR_norm_rr_Baseline', 'PR_norm_rr_AttractorEvoked', 'PR_norm_rr_AttractorRecovery', ...
    ... % -- recurrence rate: self (own-epsilon) per fixed epoch, and cross --
    'RQA_ev_RR', 'RQA_re_RR', 'cross_recur_density', ...
    'epsilon_rr', 'epsilon_re', ...
    ... % -- attractor onset/return timing & detection --
    't_attractor_onset', 't_attractor_return', ...
    'onset_detected', 'return_detected' ...
};
end

% =============================================================================
function val = fn_getFieldOrNaN(S, fieldName)
%FN_GETFIELDORNAN  Small helper: look up fieldName in struct S, handling
%  the non-scalar cases used in the saved results (RQA_ev.RR / RQA_re.RR,
%  PR_norm(1:3), PR_norm_rr(1:3), align_mat_rr(2,3)) as flattened named
%  fields, and returning NaN if the field is absent so a missing metric
%  never breaks the pool.

switch fieldName
    case 'RQA_ev_RR'
        val = fn_safeGet(S, {'RQA_ev','RR'});
    case 'RQA_re_RR'
        val = fn_safeGet(S, {'RQA_re','RR'});
    case 'PR_norm_Baseline'
        val = fn_safeIndex(S, 'PR_norm', 1);
    case 'PR_norm_Evoked'
        val = fn_safeIndex(S, 'PR_norm', 2);
    case 'PR_norm_Recovery'
        val = fn_safeIndex(S, 'PR_norm', 3);
    case 'PR_norm_rr_Baseline'
        val = fn_safeIndex(S, 'PR_norm_rr', 1);
    case 'PR_norm_rr_AttractorEvoked'
        val = fn_safeIndex(S, 'PR_norm_rr', 2);
    case 'PR_norm_rr_AttractorRecovery'
        val = fn_safeIndex(S, 'PR_norm_rr', 3);
    case 'align_ER_rr'
        % align_mat_rr is nW x nW; (2,3) = Attractor-evoked -> Attractor-
        % recovery, matching how align_ER reads align_mat(2,3).
        val = fn_safeMatrixIndex(S, 'align_mat_rr', 2, 3);
    otherwise
        if isfield(S, fieldName)
            val = S.(fieldName);
        else
            val = NaN;
        end
end
end

function val = fn_safeGet(S, fieldChain)
val = NaN;
if isfield(S, fieldChain{1}) && isfield(S.(fieldChain{1}), fieldChain{2})
    val = S.(fieldChain{1}).(fieldChain{2});
end
end

function val = fn_safeIndex(S, fieldName, idx)
val = NaN;
if isfield(S, fieldName) && numel(S.(fieldName)) >= idx
    val = S.(fieldName)(idx);
end
end

function val = fn_safeMatrixIndex(S, fieldName, row, col)
val = NaN;
if isfield(S, fieldName)
    M = S.(fieldName);
    if size(M,1) >= row && size(M,2) >= col
        val = M(row, col);
    end
end
end