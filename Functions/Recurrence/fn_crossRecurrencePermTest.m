function [p_val, null_dist, z_score] = fn_crossRecurrencePermTest(trajA, trajB, epsilon, nPerm)
% FN_CROSSRECURRENCEPERMTEST  Permutation null for cross-recurrence density.
%
%   Null hypothesis: the group labels ("epoch A" vs "epoch B") carry no
%   attractor information -- i.e. any random split of the SAME pooled
%   points into groups of sizes nA/nB would give a similar cross-density.
%   Pool trajA+trajB, repeatedly re-split into random groups of the
%   original sizes (no time info preserved), recompute cross-recurrence
%   density each time -> null distribution.
%
%   Interpretation: this test has most power to catch a FALSE "same
%   attractor" claim -- if A and B are genuinely separated in state space,
%   mixing them randomly will pull points from each half close to points
%   from the other half far more often than the true partition does, so
%   observed << null, z very negative. For a TRUE "same attractor" case,
%   z near 0 is CONSISTENT WITH shared structure, not evidence against it.
 
if nargin < 4, nPerm = 1000; end
nA = size(trajA,1); nB = size(trajB,1);
pooled = [trajA; trajB];
nTotal = nA + nB;
 
observed = fn_crossRecurrenceDensity(trajA, trajB, epsilon);
 
null_dist = zeros(nPerm,1);
for p = 1:nPerm
    permIdx = randperm(nTotal);
    null_dist(p) = fn_crossRecurrenceDensity(pooled(permIdx(1:nA),:), pooled(permIdx(nA+1:end),:), epsilon);
end
 
p_val   = mean(null_dist <= observed);   % one-sided: observed lower than chance (evidence of true separation)
z_score = (observed - mean(null_dist)) / std(null_dist);
end
