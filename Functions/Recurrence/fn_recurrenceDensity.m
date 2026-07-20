function [pRecur, delay_s, recurs, tested] = fn_recurrenceDensity(traj, t_s, epsilon, minLagS, maxWindowS)
% FN_RECURRENCEDENSITY  Per-point forward-looking recurrence density.
%   (Refactored from Bruno et al 2017 eLife, recurrence density,
%   Fig 3G / Methods.) Asks, per point: "does this point ever come back
%   close to where it is now, at some LATER time?" For a bounded,
%   stationary/periodic trajectory the answer is yes for almost every
%   point (pRecur -> 1); it drops below 1 for transient/drifting points
%   that never return.
%
% INPUTS
%   traj       - N x D matrix, trajectory points in time order
%   t_s        - N x 1 vector, times (s) matching traj (must be increasing)
%   epsilon    - recurrence threshold (Euclidean distance, same units as traj)
%   minLagS    - (optional) minimum delay (s) required to count as a genuine
%                return, on top of automatic removal of the leading
%                contiguous "still nearby" run (default: 0)
%   maxWindowS - (optional) don't test points that have less than this many
%                seconds of trajectory remaining after them (default: 0)
%
% OUTPUTS
%   pRecur   - scalar in [0 1]: fraction of TESTED points that recur at least once
%   delay_s  - N x 1, shortest recurrence delay (s) per point (NaN if none / untested)
%   recurs   - N x 1 logical, true for tested points that do recur
%   tested   - N x 1 logical, true for points with enough future to be tested
%
% NOTE: cost is O(N^2). Downsample traj/t_s before calling this on long
% recordings.
 
if nargin < 4 || isempty(minLagS),    minLagS = 0; end
if nargin < 5 || isempty(maxWindowS), maxWindowS = 0; end
 
N       = size(traj, 1);
delay_s = nan(N, 1);
recurs  = false(N, 1);
tested  = false(N, 1);
 
for i = 1:N-1
    if (t_s(end) - t_s(i)) < maxWindowS
        continue;                                  % not enough future left to fairly test
    end
    tested(i) = true;
 
    futureIdx = (i+1):N;                            % all later time-points
    if isempty(futureIdx), continue; end
 
    d      = sqrt(sum((traj(futureIdx,:) - traj(i,:)).^2, 2));
    within = d <= epsilon;
 
    if ~any(within)
        continue;                                   % never comes back close
    end
 
    % Remove the leading contiguous run of "within epsilon" points: this is
    % just the local tangent of the trajectory immediately after i (it
    % hasn't left the epsilon-ball yet), not a genuine return.
    firstGap = find(~within, 1, 'first');
    if isempty(firstGap)
        continue;                                   % stays within epsilon for the rest of the recording
    end
 
    candidates = within;
    candidates(1:firstGap-1) = false;                % drop the leading run
 
    % also enforce the explicit minimum lag, if any survives it
    if minLagS > 0
        tooSoon = (t_s(futureIdx) - t_s(i)) < minLagS;
        candidates(tooSoon) = false;
    end
 
    relIdx = find(candidates, 1, 'first');
    if isempty(relIdx)
        continue;                                   % no genuine return survives
    end
 
    hitIdx     = futureIdx(relIdx);
    recurs(i)  = true;
    delay_s(i) = t_s(hitIdx) - t_s(i);
end
 
if any(tested)
    pRecur = mean(recurs(tested));
else
    pRecur = NaN;
end
 
end
