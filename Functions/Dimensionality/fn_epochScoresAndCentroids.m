function [epoch_scores, centroids] = fn_epochScoresAndCentroids(spike_conv, wins_s, fs, V, mu, win_labels)
% FN_EPOCHSCORESANDCENTROIDS  Project each epoch onto the top-3 PCs of a
% shared basis, and report each epoch's centroid in that space.
%
%   [epoch_scores, centroids] = FN_EPOCHSCORESANDCENTROIDS(spike_conv, wins_s, fs, V, mu, win_labels)
%
% INPUTS
%   spike_conv - nN x nFrames smoothed activity
%   wins_s     - nW x 2 [start end] epoch windows, in seconds
%   fs         - sampling rate (fps)
%   V          - nN x nN PCA eigenvector basis (from fn_computePCA)
%   mu         - nN x 1 mean used to centre the data (must match V)
%   win_labels - 1 x nW cell array of epoch names, for logging
%
% OUTPUTS
%   epoch_scores - 1 x nW cell array, each element T_w x 3 (PC1-3 scores)
%   centroids    - 1 x nW cell array, each element 1 x 3 (mean PC1-3)

nFrames = size(spike_conv, 2);
nW = size(wins_s, 1);
epoch_scores = cell(1, nW);
centroids    = cell(1, nW);

for w = 1:nW
    f0 = max(round(wins_s(w,1)*fs), 1);
    f1 = min(round(wins_s(w,2)*fs), nFrames);
    seg_raw = double(spike_conv(:, f0:f1)) - double(mu);
    epoch_scores{w} = (V(:,1:3)' * seg_raw)';
    centroids{w}    = mean(epoch_scores{w}, 1);
    if nargin >= 6 && ~isempty(win_labels)
        fprintf('  [%s] %d frames, centroid: [%.4f %.4f %.4f]\n', ...
            win_labels{w}, size(epoch_scores{w},1), centroids{w});
    end
end

end
