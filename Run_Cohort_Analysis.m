%% RUN_COHORT_ANALYSIS
%  Run this AFTER Run_Attractor_Analysis.m has been run for every animal.
%  Pools the listed '<recording_ID>_results.mat' files into one table,
%  ready for you to run your own statistics on.

addpath('CohortAnalysis');

RESULTS_DIR = 'Results';   % <-- EDIT if your results are saved somewhere else

% List exactly which recordings belong in this cohort -
recording_IDs = {'Animal1_Trial1', 'Animal2_Trial1', 'Animal3_Trial1'};   % <-- EDIT (don't include _results.mat)

% To pool every '*_results.mat' file in RESULTS_DIR instead, comment the
% line above out and uncomment this one:
% recording_IDs = {};

cohort = fn_cohortAnalysis(RESULTS_DIR, recording_IDs);
cohort.summaryTable
