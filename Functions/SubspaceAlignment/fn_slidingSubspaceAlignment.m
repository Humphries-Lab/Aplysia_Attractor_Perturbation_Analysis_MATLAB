function [al_t, al_v] = fn_slidingSubspaceAlignment(spike_conv, reference_axes, win_f, step_f, fs, nDims)
% FN_SLIDINGSUBSPACEALIGNMENT  Alignment of the instantaneous (sliding
% window) activity subspace with a fixed reference subspace, over time.
%
%   [al_t, al_v] = FN_SLIDINGSUBSPACEALIGNMENT(spike_conv, reference_axes, win_f, step_f, fs, nDims)
%
% INPUTS
%   spike_conv     - nN x nFrames smoothed activity
%   reference_axes - nN x nDims reference subspace (e.g. Evoked-epoch axes)
%   win_f          - window length, in frames
%   step_f         - step size, in frames
%   fs             - sampling rate (fps)
%   nDims          - subspace dimension K (must match size(reference_axes,2))
%
% OUTPUTS
%   al_t - 1 x nWin window-centre times (s)
%   al_v - 1 x nWin alignment of each window's top-K subspace with reference_axes

nFrames = size(spike_conv, 2);
al_t = [];
al_v = [];
for s = 1:step_f:(nFrames-win_f+1)
    seg = double(spike_conv(:, s:s+win_f-1));
    seg = seg - mean(seg, 2);
    [Vs, ~] = eig((seg*seg')/(win_f-1));
    Vs = Vs(:, end:-1:1);   % ascending -> descending eigenvalue order
    al_t(end+1) = (s + win_f/2) / fs; %#ok<AGROW>
    al_v(end+1) = fn_subspaceAlignment(reference_axes, Vs(:,1:nDims)); %#ok<AGROW>
end

end
