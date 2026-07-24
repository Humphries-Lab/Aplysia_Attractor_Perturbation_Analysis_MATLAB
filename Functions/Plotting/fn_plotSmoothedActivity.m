function fn_plotSmoothedActivity(spike_conv, fs, nN, nFrames, sigma_s, t_P9, t_C2, recording_ID, figuresDir)
% FN_PLOTSMOOTHEDACTIVITY  Heatmap of the Gaussian-smoothed population
% activity, with P9/C2 landmark lines.
%
%   FN_PLOTSMOOTHEDACTIVITY(spike_conv, fs, nN, nFrames, sigma_s, t_P9, t_C2, recording_ID, figuresDir)

fig = figure('Name','Smoothed Activity','Position',[50 50 1400 400]);
imagesc((0:nFrames-1)/fs, 1:nN, spike_conv); colormap(hot); axis tight; colorbar;
xlabel('Time (s)'); ylabel('Neuron #');
title(sprintf('Gaussian-smoothed activity (sigma=%.1f ms) | %s', sigma_s*1000, recording_ID));
hold on;
xline(t_P9,'w-','LineWidth',2,'Label','P9');
xline(t_C2,'c-','LineWidth',2,'Label','C2');
exportgraphics(fig, fullfile(figuresDir,'01c_smoothed.png'), 'Resolution', 500);

end
