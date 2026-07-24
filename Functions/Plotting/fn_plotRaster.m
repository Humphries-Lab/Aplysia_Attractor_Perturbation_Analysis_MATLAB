function fn_plotRaster(peaks, fs, nN, nFrames, t_P9, t_C2, recording_ID, protocol, figuresDir)
% FN_PLOTRASTER  Spike raster (top) and 1-s-binned population activity
% (bottom), with P9/C2 landmark lines.
%
%   FN_PLOTRASTER(peaks, fs, nN, nFrames, t_P9, t_C2, recording_ID, protocol, figuresDir)

fig = figure('Name','Spike Raster','Position',[50 50 1400 600]);
subplot(2,1,1);
[row, col] = find(peaks);
scatter(col/fs, row, 1, 'k', 'filled');
xlim([0 nFrames/fs]); ylim([0.5 nN+0.5]); axis tight; box on;
xlabel('Time (s)'); ylabel('Neuron #');
title(sprintf('Spike raster (n=%d) | %s | %s', nN, recording_ID, protocol));
hold on;
xline(t_P9,'r-','LineWidth',2,'Label','P9 (motor onset)');
xline(t_C2,'b-','LineWidth',2,'Label','C2 (perturbation)');

subplot(2,1,2);
bin_f = round(fs); n_bins = floor(nFrames/bin_f);
pop_rate = zeros(1,n_bins);
for b = 1:n_bins
    pop_rate(b) = sum(sum(peaks(:,(b-1)*bin_f+1:b*bin_f)));
end
bar(0.5:n_bins-0.5, pop_rate, 1, 'FaceColor',[0.4 0.4 0.8],'EdgeColor','none');
hold on;
xline(t_P9,'r-','LineWidth',2,'Label','P9');
xline(t_C2,'b-','LineWidth',2,'Label','C2');
xlabel('Time (s)'); ylabel('Population events / s');
title('Population activity (1 s bins)'); grid on; xlim([0 nFrames/fs]);
exportgraphics(fig, fullfile(figuresDir,'01b_raster.png'), 'Resolution', 500);

end
