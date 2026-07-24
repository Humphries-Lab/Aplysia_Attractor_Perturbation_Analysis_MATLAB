function fn_plotPRTimeseries(pr_t, pr_v, t_P9, t_evoked_start, t_C2, t_recovery_start, ...
    motor_buffer_s, recovery_delay_s, slide_win_s, slide_step_s, recording_ID, figuresDir)
% FN_PLOTPRTIMESERIES  Sliding-window normalised participation ratio (PR/N)
% over the whole recording, with epoch-landmark lines.
%
%   FN_PLOTPRTIMESERIES(pr_t, pr_v, t_P9, t_evoked_start, t_C2, t_recovery_start, ...
%       motor_buffer_s, recovery_delay_s, slide_win_s, slide_step_s, recording_ID, figuresDir)

fig = figure('Name','PR over time','Position',[100 400 1200 400]);
plot(pr_t, pr_v, 'k-', 'LineWidth', 1.5); hold on;
xline(t_P9,'r--','LineWidth',2,'Label','P9');
xline(t_evoked_start,'r:','LineWidth',1.5,'Label',sprintf('+%ds',motor_buffer_s));
xline(t_C2,'b--','LineWidth',2,'Label','C2');
xline(t_recovery_start,'b:','LineWidth',1.5,'Label',sprintf('+%ds',recovery_delay_s));
xlabel('Time (s)'); ylabel(sprintf('PR/N  (window=%ds)',slide_win_s)); grid on;
title(sprintf('Normalised dimensionality over time | %s  [window=%ds, step=%ds]', ...
              recording_ID, slide_win_s, slide_step_s));
exportgraphics(fig, fullfile(figuresDir,'03_PR_timeseries.png'), 'Resolution', 500);

end
