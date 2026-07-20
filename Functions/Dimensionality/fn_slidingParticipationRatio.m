function [pr_t, pr_v] = fn_slidingParticipationRatio(spike_conv, win_f, step_f, fs, nN)
% FN_SLIDINGPARTICIPATIONRATIO  Normalised participation ratio (PR/N) in a
% sliding window across the whole recording.
%
%   [pr_t, pr_v] = FN_SLIDINGPARTICIPATIONRATIO(spike_conv, win_f, step_f, fs, nN)
%
% INPUTS
%   spike_conv - nN x nFrames smoothed activity
%   win_f      - window length, in frames
%   step_f     - step size, in frames
%   fs         - sampling rate (fps)
%   nN         - number of neurons (for normalisation of PR)
%
% OUTPUTS
%   pr_t - 1 x nWin window-centre times (s)
%   pr_v - 1 x nWin PR/N values

nFrames = size(spike_conv, 2);
pr_t = [];
pr_v = [];
for s = 1:step_f:(nFrames-win_f+1)
    seg   = double(spike_conv(:, s:s+win_f-1));
    seg   = seg - mean(seg, 2);
    ev_s  = sort(eig((seg*seg')/(win_f-1)), 'descend');
    pr_t(end+1) = (s + win_f/2) / fs; %#ok<AGROW>
    pr_v(end+1) = fn_participationRatio(ev_s) / nN; %#ok<AGROW>
end

end
