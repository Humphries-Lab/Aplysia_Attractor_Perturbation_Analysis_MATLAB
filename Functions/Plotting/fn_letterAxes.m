function fn_letterAxes(axArray, letterSize, fontName)
% FN_LETTERAXES  A, B, C... labels flush to the top-left of each axes.
letters = {'A','B','C','D','E','F','G','H','I','J','K','L','M'};
for i = 1:numel(axArray)
    text(axArray(i), -0.08, 1.15, letters{i}, ...
        'Units','normalized', 'FontSize', letterSize, 'FontWeight','bold', ...
        'FontName', fontName, 'Clipping','off', 'HorizontalAlignment','right');
end
end
