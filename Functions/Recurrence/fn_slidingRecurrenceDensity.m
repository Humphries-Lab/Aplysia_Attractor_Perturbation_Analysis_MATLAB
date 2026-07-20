function [win_t, win_v] = fn_slidingRecurrenceDensity(t_axis, recurs, tested, f0, f1, fs, win_f, step_f)
% FN_SLIDINGRECURRENCEDENSITY  Bin a per-point recurrence test (from
% fn_recurrenceDensity) into a sliding-window time series, i.e. the
% fraction of TESTED points recurring within each window.
%
%   [win_t, win_v] = FN_SLIDINGRECURRENCEDENSITY(t_axis, recurs, tested, f0, f1, fs, win_f, step_f)
%
% INPUTS
%   t_axis  - N x 1 times (s) of the (downsampled) trajectory points
%   recurs  - N x 1 logical, per-point recurrence flag (from fn_recurrenceDensity)
%   tested  - N x 1 logical, per-point tested flag (from fn_recurrenceDensity)
%   f0, f1  - frame range (in the ORIGINAL, full-fs frame grid) over which
%             to slide the window
%   fs      - sampling rate (fps)
%   win_f   - window length, in frames
%   step_f  - step size, in frames
%
% OUTPUTS
%   win_t - 1 x nWin window-centre times (s)
%   win_v - 1 x nWin fraction of tested points recurring in each window
%           (NaN where no point in that window was tested)

win_t = [];
win_v = [];
t_axis = t_axis(:)';
recurs = recurs(:)';
tested = tested(:)';

for s = f0:step_f:(f1-win_f)
    t0 = s / fs;
    t1 = (s + win_f) / fs;
    t_centre = (t0 + t1) / 2;
    inWin = t_axis >= t0 & t_axis < t1;
    testedInWin = inWin & tested;
    win_t(end+1) = t_centre; %#ok<AGROW>
    if any(testedInWin)
        win_v(end+1) = mean(recurs(testedInWin)); %#ok<AGROW>
    else
        win_v(end+1) = NaN; %#ok<AGROW>
    end
end

end
