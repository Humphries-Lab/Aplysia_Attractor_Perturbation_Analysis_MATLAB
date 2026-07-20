function sigma_s = fn_estimateKernelWidth(peaks, fs)
% FN_ESTIMATEKERNELWIDTH  Data-driven Gaussian smoothing width from the
% median inter-spike interval (ISI) pooled across all neurons.
%
%   sigma_s = FN_ESTIMATEKERNELWIDTH(peaks, fs)
%
% INPUTS
%   peaks - nN x nFrames logical event matrix
%   fs    - sampling rate (fps)
%
% OUTPUT
%   sigma_s - Gaussian kernel standard deviation, in seconds (median ISI)

nN  = size(peaks, 1);
ISI = [];
for k = 1:nN
    sp = find(peaks(k,:));
    if numel(sp) > 1
        ISI = [ISI, diff(sp)/fs]; %#ok<AGROW>
    end
end
sigma_s = median(ISI);
fprintf('  Gaussian sigma = %.4f s  (%.1f frames)  [median ISI]\n', sigma_s, sigma_s*fs);

end
