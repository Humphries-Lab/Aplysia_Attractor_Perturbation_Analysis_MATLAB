function fn_plotTrajectory3D(scores, t_pca, t_C2, var_exp, recording_ID, figuresDir)
% FN_PLOTTRAJECTORY3D  3-D population trajectory in the top-3 PC space,
% coloured by pre-/post-C2, with markers for the trajectory start and the
% C2 perturbation onset.
%
%   FN_PLOTTRAJECTORY3D(scores, t_pca, t_C2, var_exp, recording_ID, figuresDir)
%
% INPUTS
%   scores - 3 x T top-3 PC scores (from fn_projectOntoPCs)
%   t_pca  - 1 x T times (s) matching scores
%   t_C2   - C2 perturbation onset time (s)
%   var_exp- percent variance explained per PC (for axis labels)

pre  = t_pca < t_C2;
post = t_pca >= t_C2;

fig = figure('Name','3D Trajectory','Position',[200 100 900 800]);
plot3(scores(1,pre),scores(2,pre),scores(3,pre),'r-','LineWidth',1.2); hold on;
plot3(scores(1,post),scores(2,post),scores(3,post),'b-','LineWidth',1.2);
scatter3(scores(1,1),scores(2,1),scores(3,1),100,'gs','filled');
idx_c2 = find(post,1);
if ~isempty(idx_c2)
    scatter3(scores(1,idx_c2),scores(2,idx_c2),scores(3,idx_c2),100,'b^','filled');
end
xlabel(sprintf('PC1 (%.1f%%)',var_exp(1)));
ylabel(sprintf('PC2 (%.1f%%)',var_exp(2)));
zlabel(sprintf('PC3 (%.1f%%)',var_exp(3)));
legend('Pre-C2','Post-C2','P9 onset','C2 onset'); grid on; view([-35 25]);
title(sprintf('Population trajectory | %s', recording_ID));
exportgraphics(fig, fullfile(figuresDir,'02_trajectory_3D.png'), 'Resolution', 500);
close(fig);

end
