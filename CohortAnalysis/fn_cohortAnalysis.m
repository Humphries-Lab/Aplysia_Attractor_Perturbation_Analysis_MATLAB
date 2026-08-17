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
%  per run. This function pools those per-recording .mat files across animals
%  so the results are ready for cohort-level statistical testing.
%
%  It does not decide or run any statistical test. It only:
%    1. Loads every requested '<recording_ID>_results.mat' from
%       RESULTS_DIR
%    2. Extracts a fixed set of scalar summary metrics per recording
%       into one tidy row-per-animal table
%    3. Returns that table (plus the untouched raw per-recording
%       structs) 
%
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
%  cohort.raw          - 1xN cell array, the untouched contents of
%                        every loaded '*_results.mat' file, indexed in
%                        the same order as summaryTable's rows (for
%                        anything not captured as a scalar summary,
%                        e.g. full pr_t/pr_v or al_t/al_v timeseries).
%  cohort.missing      - cellstr of recording_IDs that were requested
%                        but had no matching '*_results.mat' file, or
%                        whose file failed to load.


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
% Stored as a cell array, not a struct array: different '<ID>_results.mat'
% files can carry slightly different sets of saved fields (e.g. if they
% were produced by different versions of Run_Attractor_Analysis.m), and
% struct-array concatenation requires an identical field set across all
% elements. A cell array has no such requirement. Index as cohort.raw{i}.
cohort.raw     = rawList;
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
