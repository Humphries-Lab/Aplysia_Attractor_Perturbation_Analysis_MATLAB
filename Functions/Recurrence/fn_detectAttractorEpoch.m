function epoch = fn_detectAttractorEpoch(spike_conv, fs, V, mu, searchStart_s, searchEnd_s, ...
    epsilon_fixed, minLagS, maxWindowS, win_f, step_f, onsetRatio, maxFullPts, fallback_s)
% FN_DETECTATTRACTOREPOCH  Locate the first time within a search window at
% which the population trajectory "locks onto" an attractor, defined as
% the first sliding window whose recurrence density reaches onsetRatio.
%
% Two epsilon thresholds are evaluated for comparison: a threshold
% calibrated WITHIN the search window itself (primary, used for
% detection) and a fixed threshold carried over from elsewhere (e.g. the
% Evoked epoch), kept only as a plotted comparison.
%
%   epoch = FN_DETECTATTRACTOREPOCH(spike_conv, fs, V, mu, searchStart_s, searchEnd_s, ...
%               epsilon_fixed, minLagS, maxWindowS, win_f, step_f, onsetRatio, maxFullPts, fallback_s)
%
% INPUTS
%   spike_conv     - nN x nFrames smoothed activity
%   fs             - sampling rate (fps)
%   V, mu          - shared PCA basis and mean (top-3 PCs used for trajectory)
%   searchStart_s, searchEnd_s - time window (s) to search within
%   epsilon_fixed  - fixed epsilon (e.g. calibrated from the Evoked epoch),
%                    kept only as a comparison trace
%   minLagS        - minimum lag (s) for a genuine recurrence (see fn_recurrenceDensity)
%   maxWindowS     - don't test points with less than this much future left
%   win_f, step_f  - sliding window length / step, in frames
%   onsetRatio     - recurrence-density threshold (e.g. 0.90) defining "locked on"
%   maxFullPts     - downsample target for the per-point recurrence test
%   fallback_s     - time (s) to report if no window reaches onsetRatio
%
% OUTPUT
%   epoch - struct with fields:
%     .t_lock          - detected lock-on time (s), or fallback_s if not detected
%     .detected        - logical, whether a threshold-crossing window was found
%     .epsilon_within  - epsilon calibrated within the search window
%     .win_t           - window-centre times (s)
%     .win_v           - within-epoch-epsilon recurrence density per window
%     .win_v_fixed     - fixed-epsilon recurrence density per window (comparison)

f0 = max(round(searchStart_s * fs), 1);
f1 = min(round(searchEnd_s   * fs), size(spike_conv,2));

[traj, t_axis] = fn_getEpochTrajectory(spike_conv, f0, f1, fs, V, mu, maxFullPts);

epsilon_within = prctile(pdist(traj,'euclidean'), 10);
fprintf('  Search eps: within-epoch=%.4f | fixed=%.4f\n', epsilon_within, epsilon_fixed);

[~, ~, recurs_local, tested]  = fn_recurrenceDensity(traj, t_axis, epsilon_within, minLagS, maxWindowS);
[~, ~, recurs_fixed, ~]       = fn_recurrenceDensity(traj, t_axis, epsilon_fixed,  minLagS, maxWindowS);

[win_t, win_v]       = fn_slidingRecurrenceDensity(t_axis, recurs_local, tested, f0, f1, fs, win_f, step_f);
[~,     win_v_fixed] = fn_slidingRecurrenceDensity(t_axis, recurs_fixed, tested, f0, f1, fs, win_f, step_f);

idx = find(win_v >= onsetRatio, 1, 'first');
detected = ~isempty(idx);
if detected
    t_lock = win_t(idx);
else
    t_lock = fallback_s;
end

epoch.t_lock         = t_lock;
epoch.detected       = detected;
epoch.epsilon_within = epsilon_within;
epoch.win_t          = win_t;
epoch.win_v          = win_v;
epoch.win_v_fixed    = win_v_fixed;

end
