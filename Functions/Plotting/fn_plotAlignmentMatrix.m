function fn_plotAlignmentMatrix(align_mat, win_labels, recording_ID, figuresDir)
% FN_PLOTALIGNMENTMATRIX  Heatmap of the pairwise subspace-alignment
% matrix across epochs, with numeric annotations.
%
%   FN_PLOTALIGNMENTMATRIX(align_mat, win_labels, recording_ID, figuresDir)

nW = size(align_mat,1);
fig = figure('Name','Alignment Matrix');
imagesc(align_mat); clim([0 1]); colormap(hot); colorbar; axis square;
set(gca,'XTick',1:nW,'XTickLabel',win_labels,'YTick',1:nW,'YTickLabel',win_labels);
for i=1:nW
    for j=1:nW
        text(j,i,sprintf('%.3f',align_mat(i,j)),'HorizontalAlignment','center', ...
             'FontSize',10,'Color','w','FontWeight','bold');
    end
end
title(sprintf('Subspace alignment | %s', recording_ID));
exportgraphics(fig, fullfile(figuresDir,'04a_alignment_matrix.png'), 'Resolution', 500);
close(fig);

end
