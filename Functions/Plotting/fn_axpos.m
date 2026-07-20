function pos = fn_axpos(row, col, rowspan, colspan, rowFracs, margin, gapX, gapY)
% FN_AXPOS  Normalized [left bottom width height] for a manual axes grid.
%   rowFracs is a vector of row-height fractions (renormalized to sum 1),
%   letting rows differ in height (e.g. a tall time-series row, tall
%   matrix row, short summary row). All layouts here use 4 columns.
nCols = 4;
rowFracs = rowFracs / sum(rowFracs);
nRows = numel(rowFracs);
 
plotW = 1 - 2*margin - (nCols-1)*gapX;
plotH = 1 - 2*margin - (nRows-1)*gapY;
 
colW = plotW / nCols;
left  = margin + (col-1)*(colW + gapX);
width = colspan*colW + (colspan-1)*gapX;
 
rowH_px = rowFracs * plotH;
top = margin + sum(rowH_px(1:row-1)) + (row-1)*gapY;
height = sum(rowH_px(row:row+rowspan-1)) + (rowspan-1)*gapY;
bottom = 1 - top - height;
 
pos = [left, bottom, width, height];
end
