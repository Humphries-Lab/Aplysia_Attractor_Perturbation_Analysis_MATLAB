function fn_formatAxes(figHandle, fontName, fontSize)
% FN_FORMATAXES  Uniform publication formatting across all axes in a figure.
allAx = findobj(figHandle, 'Type', 'axes');
set(allAx, 'FontName', fontName, 'FontSize', fontSize, ...
    'LineWidth', 1.0, 'TickDir', 'out', 'TickLength', [0.01 0.01]);
end
