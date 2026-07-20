function fn_plotRateDistribution(rate_all, recording_ID, nNeurons, figuresDir)
% FN_PLOTRATEDISTRIBUTION  Histogram of per-neuron mean event rates, with
% the silent/active reference lines used for quality filtering.
%
%   FN_PLOTRATEDISTRIBUTION(rate_all, recording_ID, nNeurons, figuresDir)

fig = figure('Name','Rate Distribution','Position',[100 100 800 450]);
histogram(rate_all(rate_all>0), 50, 'FaceColor',[0.3 0.5 0.8], 'EdgeColor','none');
xlabel('Mean event rate (Hz)'); ylabel('Number of neurons');
title(sprintf('Firing rate distribution | %s  (n=%d)', recording_ID, nNeurons));
rl = xline(0.01,'r--','LineWidth',1.5);
gl = xline(0.10,'g--','LineWidth',1.5);
legend([rl, gl], ...
    {'Min. threshold (0.01 Hz) - neurons below this are silent and excluded', ...
     '0.10 Hz reference - neurons above this are active'}, ...
    'Location','northeast','FontSize',9);
grid on;
exportgraphics(fig, fullfile(figuresDir,'01a_rate_distribution.png'), 'Resolution', 500);
close(fig);

end
