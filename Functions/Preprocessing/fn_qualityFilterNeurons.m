function [peaks_out, good, rate_all, nN] = fn_qualityFilterNeurons(peaks, fs, min_rate, max_rate)
% FN_QUALITYFILTERNEURONS  Drop silent/pathological neurons by mean event rate.
%
%   [peaks_out, good, rate_all, nN] = FN_QUALITYFILTERNEURONS(peaks, fs, min_rate, max_rate)
%
% INPUTS
%   peaks    - nNeurons x nFrames logical event matrix
%   fs       - sampling rate (fps)
%   min_rate - minimum mean event rate (Hz) to keep a neuron
%   max_rate - maximum mean event rate (Hz) to keep a neuron
%
% OUTPUTS
%   peaks_out - filtered event matrix (nN x nFrames)
%   good      - nNeurons x 1 logical mask of kept neurons
%   rate_all  - nNeurons x 1 mean event rate (Hz), computed on the FULL
%               (pre-filter) population
%   nN        - number of neurons kept

nNeurons = size(peaks, 1);
nFrames  = size(peaks, 2);

rate_all = sum(peaks, 2) / (nFrames/fs);
fprintf('\n  Rate distribution (%d neurons):\n', nNeurons);
fprintf('  Silent=0: %d  |  <0.01Hz: %d  |  0.01-0.1Hz: %d  |  >0.1Hz: %d\n', ...
    sum(rate_all==0), sum(rate_all>0 & rate_all<0.01), ...
    sum(rate_all>=0.01 & rate_all<0.1), sum(rate_all>=0.1));

good      = rate_all >= min_rate & rate_all <= max_rate;
peaks_out = peaks(good, :);
nN        = sum(good);

fprintf('  Kept %d / %d neurons (%.3f-%.0f Hz)\n', nN, nNeurons, min_rate, max_rate);

end
