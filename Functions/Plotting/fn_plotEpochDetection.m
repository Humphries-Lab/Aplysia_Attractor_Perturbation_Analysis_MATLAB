function fn_plotEpochDetection(onset, return_, t_P9, t_C2, t_end, onsetRatio, recording_ID, figuresDir)
% FN_PLOTEPOCHDETECTION  Onset-search and return-search recurrence-density
% traces, with detected epoch windows shaded and the onset threshold
% marked.
%
%   FN_PLOTEPOCHDETECTION(onset, return_, t_P9, t_C2, t_end, onsetRatio, recording_ID, figuresDir)
%
% INPUTS
%   onset, return_ - structs from fn_detectAttractorEpoch (P9->C2 and
%                    C2->end searches, respectively)

fig = figure('Name','Epoch detection','Position',[100 100 1200 420]);
hold on;
y_lo = 0; y_hi = 1.05;
if onset.detected
    patch([onset.t_lock t_C2 t_C2 onset.t_lock], ...
        [y_lo y_lo y_hi y_hi],[1 0.87 0.87],'EdgeColor','none','DisplayName','Evoked');
end
if return_.detected
    patch([return_.t_lock t_end t_end return_.t_lock], ...
        [y_lo y_lo y_hi y_hi],[0.87 0.91 1],'EdgeColor','none','DisplayName','Recovery');
end
plot([0 t_end],[onsetRatio onsetRatio],'--','Color',[0.4 0.4 0.4],'LineWidth',1.5, ...
    'DisplayName',sprintf('%.0f%% criterion', onsetRatio*100));
plot(onset.win_t,   onset.win_v,       'r-','LineWidth',1.8,'DisplayName','P9->C2 (within-epoch eps)');
plot(return_.win_t, return_.win_v,     'b-','LineWidth',1.8,'DisplayName','C2->end (within-epoch eps)');
plot(onset.win_t,   onset.win_v_fixed, 'r:','LineWidth',1.2,'DisplayName','P9->C2 (fixed evoked eps, comparison)');
plot(return_.win_t, return_.win_v_fixed,'b:','LineWidth',1.2,'DisplayName','C2->end search (fixed evoked eps, comparison)');
xline(t_P9,'r-','LineWidth',2,'Label','P9','HandleVisibility','off');
xline(t_C2,'b-','LineWidth',2,'Label','C2','HandleVisibility','off');
ylim([y_lo y_hi]); xlim([0 t_end]);
xlabel('Time (s)'); ylabel('Recurrence density');
title(sprintf('%s | Attractor epoch detection', recording_ID), 'FontWeight', 'normal');
legend('Location','northeastoutside'); grid off; box on;
exportgraphics(fig, fullfile(figuresDir, '05c_epoch_detection.png'), 'Resolution', 550);
close(fig);

end
