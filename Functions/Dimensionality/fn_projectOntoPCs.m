function scores = fn_projectOntoPCs(spike_conv, f_start, f_end, V, mu, nDims, blockSize)
% FN_PROJECTONTOPCS  Project a frame window of activity onto the top
% nDims columns of a (previously computed) PCA basis.
%
%   scores = FN_PROJECTONTOPCS(spike_conv, f_start, f_end, V, mu, nDims, blockSize)
%
% INPUTS
%   spike_conv - nN x nFrames smoothed activity
%   f_start,f_end - frame window (inclusive) to project
%   V          - nN x nN (or nN x >=nDims) eigenvector basis, from fn_computePCA
%   mu         - nN x 1 mean used to centre the data (must match the basis)
%   nDims      - number of leading components to project onto
%   blockSize  - (optional) columns processed per block (default 5000)
%
% OUTPUT
%   scores - nDims x (f_end-f_start+1) single-precision projected scores

if nargin < 7 || isempty(blockSize)
    blockSize = 5000;
end

Vk    = V(:, 1:nDims);
nT    = f_end - f_start + 1;
scores = zeros(nDims, nT, 'single');
col   = 1;
for b0 = f_start:blockSize:f_end
    b1  = min(b0 + blockSize - 1, f_end);
    blk = b1 - b0 + 1;
    seg = double(spike_conv(:, b0:b1)) - double(mu);
    scores(:, col:col+blk-1) = single(Vk' * seg);
    col = col + blk;
end

end
