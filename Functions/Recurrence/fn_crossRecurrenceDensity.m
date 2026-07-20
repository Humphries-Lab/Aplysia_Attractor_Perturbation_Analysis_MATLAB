function [pRecur, recurs, minDist] = fn_crossRecurrenceDensity(trajA, trajB, epsilon)
% FN_CROSSRECURRENCEDENSITY  Fraction of points in trajA that recur ANYWHERE
% in trajB, within epsilon. Cross-epoch analogue of fn_recurrenceDensity.
%
% INPUTS
%   trajA   - Na x D matrix (e.g. evoked epoch trajectory)
%   trajB   - Nb x D matrix (e.g. recovery epoch trajectory)
%   epsilon - recurrence threshold (should be the SAME epsilon used to
%             define trajA's own attractor, so this is a fair
%             like-for-like comparison)
%
% OUTPUTS
%   pRecur  - scalar in [0 1]: fraction of trajA points with >=1 match in trajB
%   recurs  - Na x 1 logical: which trajA points recur in trajB
%   minDist - Na x 1: distance to the closest trajB point, for each trajA point
 
Na = size(trajA, 1);
Nb = size(trajB, 1);
recurs  = false(Na, 1);
minDist = nan(Na, 1);
 
BLK = 200;
for i0 = 1:BLK:Na
    i1 = min(i0+BLK-1, Na);
    A  = trajA(i0:i1, :);
    best = inf(size(A,1), 1);
    for j0 = 1:BLK:Nb
        j1 = min(j0+BLK-1, Nb);
        B  = trajB(j0:j1, :);
        AA = sum(A.^2, 2);
        BB = sum(B.^2, 2);
        D  = sqrt(max(bsxfun(@plus, AA, BB') - 2*(A*B'), 0));
        best = min(best, min(D, [], 2));
    end
    minDist(i0:i1) = best;
    recurs(i0:i1)  = best <= epsilon;
end
 
pRecur = mean(recurs);
 
end
