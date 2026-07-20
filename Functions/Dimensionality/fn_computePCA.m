function [V, var_exp, mu, ev] = fn_computePCA(spike_conv, f_start, f_end, blockSize)
% FN_COMPUTEPCA  PCA on a specified frame window of the smoothed
% population activity, using block-accumulated covariance so the full
% (nN x nN) covariance never requires materialising the whole segment
% at once.
%
%   [V, var_exp, mu, ev] = FN_COMPUTEPCA(spike_conv, f_start, f_end, blockSize)
%
% INPUTS
%   spike_conv - nN x nFrames smoothed activity (single/double)
%   f_start    - first frame (inclusive) of the window used to build the
%                covariance (e.g. onset of the Evoked epoch, to exclude
%                Baseline from the PCA basis)
%   f_end      - last frame (inclusive) of the window (e.g. nFrames)
%   blockSize  - (optional) columns processed per block (default 5000)
%
% OUTPUTS
%   V       - nN x nN eigenvectors, columns sorted by descending eigenvalue
%   var_exp - nN x 1, percent variance explained per component
%   mu      - nN x 1, mean activity over the window (needed to project
%             any other segment onto the same basis)
%   ev      - nN x 1, eigenvalues, descending

if nargin < 4 || isempty(blockSize)
    blockSize = 5000;
end

nN = size(spike_conv, 1);
nT = f_end - f_start + 1;

mu = mean(spike_conv(:, f_start:f_end), 2);
C  = zeros(nN, nN);
for b0 = f_start:blockSize:f_end
    b1  = min(b0 + blockSize - 1, f_end);
    seg = double(spike_conv(:, b0:b1)) - double(mu);
    C   = C + seg*seg';
end
C = C / (nT - 1);

[V, D]     = eig(C);
[ev, idx]  = sort(diag(D), 'descend');
V          = V(:, idx);
var_exp    = 100 * ev / sum(ev);

end
