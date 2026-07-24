function fn_plotAlignmentTimeseries(al_t, al_v, chance_lvl, t_P9, t_evoked_start, t_C2, t_recovery_start, ...
    motor_buffer_s, recovery_delay_s, nDims, slide_win_s, recording_ID, figuresDir)
% FN_PLOTALIGNMENTTIMESERIES  Alignment of the instantaneous subspace with
% the Evoked-epoch subspace, over time, with the chance level and epoch
% landmarks marked.
%
%   FN_PLOTALIGNMENTTIMESERIES(al_t, al_v, chance_lvl, t_P9, t_evoked_start, t_C2, ...
%       t_recovery_start, motor_buffer_s, recovery_delay_s, nDims, slide_win_s, recording_ID, figuresDir)

fig = figure('Name','Alignment over time','Position',[100 400 1200 450]);
plot(al_t, al_v, 'k-', 'LineWidth', 1.5); hold on;
yline(chance_lvl,'r--','LineWidth',1.5,'Label',sprintf('Chance=%.3f',chance_lvl));
xline(t_P9,'r-','LineWidth',1.5,'Label','P9');
xline(t_evoked_start,'r:','LineWidth',1.2,'Label',sprintf('+%ds',motor_buffer_s));
xline(t_C2,'b-','LineWidth',1.5,'Label','C2');
xline(t_recovery_start,'b:','LineWidth',1.2,'Label',sprintf('+%ds',recovery_delay_s));
ylim([0 1]); xlabel('Time (s)');
ylabel('Alignment with evoked subspace'); grid on;
title(sprintf('Alignment over time | %s  [K=%d, win=%ds]', recording_ID, nDims, slide_win_s));
exportgraphics(fig, fullfile(figuresDir,'04c_alignment_timeseries.png'), 'Resolution', 500);

end
