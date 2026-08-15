function group = fn_groupLevelAnalysis(RESULTS_DIR, recording_IDs, varargin)
%FN_GROUPLEVELANALYSIS  Wrapper for group-level (across-recording/animal)
%  analysis of outputs produced by Run_Attractor_Analysis.
%
%   group = FN_GROUPLEVELANALYSIS(RESULTS_DIR, recording_IDs)
%   group = FN_GROUPLEVELANALYSIS(..., 'Name', Value, ...)
%
%  PURPOSE
%  -------
%  Run_Attractor_Analysis.m is a SINGLE-RECORDING pipeline: it runs on
%  one animal/session at a time and saves one '<recording_ID>_results.mat'
%  per run. This function is the entry point for the step that comes
%  AFTER that -- pooling those per-recording .mat files across animals
%  so the results are ready for group-level statistical testing.
%
%  This function is a WRAPPER ONLY. It does not decide or run
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
%  
