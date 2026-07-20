function fn_plotDominantOscillation(osc, figuresDir)
% FN_PLOTDOMINANTOSCILLATION  Time-domain (drift-removal) and frequency-
% domain (FFT amplitude spectrum) panels for the recurrence-density
% oscillation analysis.
%
%   FN_PLOTDOMINANTOSCILLATION(osc, figuresDir)
%
% INPUT
%   osc - struct returned by fn_dominantOscillationPeriod

if isempty(osc.f_fft)
    return;  % too few samples - nothing to plot
end

fig = figure('Name', 'FFT Recurrence Density', 'Color', 'w', 'Position', [100, 100, 700, 600]);

subplot(2, 1, 1);
plot(osc.time_vector, osc.rr_detrended, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
hold on;
plot(osc.time_vector, osc.rr_filtered, 'k', 'LineWidth', 1.2);
title(sprintf('Time Domain: High-Pass Filter (Cutoff = %d s)', osc.cutoff_period), 'FontWeight', 'normal');
xlabel('Time (s)'); ylabel('Recurrence Density');
legend('Original (Detrended)', 'Filtered Signal', 'Location', 'best');
grid on; box on; set(gca, 'FontSize', 11, 'TickDir', 'out');
hold off;

subplot(2, 1, 2);
plot(osc.f_fft, osc.P1, 'k', 'LineWidth', 1.2);
hold on;
plot(osc.dom_freq, osc.max_amp, 'ro', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
label_str = sprintf('Peak: %.3f Hz (T = %.1f s)', osc.dom_freq, osc.dom_period);
text(osc.dom_freq, osc.max_amp, ['  ' label_str], 'VerticalAlignment', 'bottom', ...
    'FontSize', 10, 'Color', 'r');
title('Single-Sided Amplitude Spectrum of Filtered Signal', 'FontWeight', 'normal');
xlabel('Frequency (Hz)'); ylabel('|P1(f)| (Amplitude)');
xlim([0, osc.Fs_rr/2]);
ylim([0, osc.max_amp * 1.2]);
grid on; box on; set(gca, 'FontSize', 11, 'TickDir', 'out');
hold off;

exportgraphics(fig, fullfile(figuresDir, '05e_fft_recurrence_period.png'), 'Resolution', 500);
close(fig);

end
