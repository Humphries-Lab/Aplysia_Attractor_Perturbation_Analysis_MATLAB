function fn_plotScree(var_exp, nN, figuresDir)
% FN_PLOTSCREE  Scree plot: individual and cumulative variance explained
% by the leading principal components.
%
%   FN_PLOTSCREE(var_exp, nN, figuresDir)

fig = figure('Name','Scree');
nShow = min(20, nN);
bar(1:nShow, var_exp(1:nShow), 'FaceColor',[0.3 0.5 0.8]); hold on;
plot(1:nShow, cumsum(var_exp(1:nShow)), 'ro-', 'MarkerFaceColor','r');
xlabel('PC'); ylabel('Variance (%)');
legend('Individual','Cumulative'); title('Scree plot'); grid on;
exportgraphics(fig, fullfile(figuresDir,'02_scree.png'), 'Resolution', 500);

end
