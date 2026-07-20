function [traj, t_axis] = fn_getEpochTrajectory(spike_conv, f0, f1, fs, V, mu, maxPts)
% FN_GETEPOCHTRAJECTORY  Extract a frame window, downsample it to at most
% maxPts points, and project onto the top-3 PCs of a shared basis. Used
% throughout the recurrence-density analyses, which all operate on a
% low-dimensional (3D) trajectory for tractable pairwise-distance costs.
%
%   [traj, t_axis] = FN_GETEPOCHTRAJECTORY(spike_conv, f0, f1, fs, V, mu, maxPts)
%
% INPUTS
%   spike_conv - nN x nFrames smoothed activity
%   f0, f1     - frame window (inclusive)
%   fs         - sampling rate (fps)
%   V          - nN x nN PCA eigenvector basis
%   mu         - nN x 1 mean used to centre the data (must match V)
%   maxPts     - target number of points after downsampling
%
% OUTPUTS
%   traj   - nPts x 3 trajectory (top-3 PC scores)
%   t_axis - nPts x 1 times (s), matching traj

f0 = max(f0, 1);
f1 = min(f1, size(spike_conv,2));
ds = max(1, floor((f1-f0+1) / maxPts));
idx = f0:ds:f1;

seg  = double(spike_conv(:, idx)) - double(mu);
traj = (V(:,1:3)' * seg)';
t_axis = idx(:) / fs;

end
