function fn_plotAlignmentVsK(K_range, al_vs_k, nN, nDims, recording_ID, figuresDir)
% FN_PLOTALIGNMENTVSK  Evoked-Recovery alignment swept across K, compared
% to the chance level K/N, with the auto/fixed K choice marked.
%
%   FN_PLOTALIGNMENTVSK(K_range, al_vs_k, nN, nDims, recording_ID, figuresDir)

fig = figure('Name','Alignment vs K');
plot(K_range, al_vs_k, 'k-o', 'LineWidth',2, 'MarkerFaceColor','k'); hold on;
plot(K_range, K_range/nN, 'r--', 'LineWidth',1.5);
xline(nDims,'b:','LineWidth',1.5,'Label',sprintf('K=%d (auto)',nDims));
xlabel('K'); ylabel('Alignment');
legend('Observed','Chance K/N','Location','southeast');
title(sprintf('Alignment sensitivity to K | %s', recording_ID)); grid on;
exportgraphics(fig, fullfile(figuresDir,'04b_alignment_vs_K.png'), 'Resolution', 500);
close(fig);

end
