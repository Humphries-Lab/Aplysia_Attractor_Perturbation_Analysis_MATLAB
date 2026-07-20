function eps = fn_safeEpsilon(traj)
% FN_SAFEEPSILON  10th-percentile pairwise-distance epsilon, with a
% fallback if the raw 10th percentile is exactly 0 (e.g. near-duplicate
% points from heavy downsampling of a slow-moving trajectory).
d = pdist(traj, 'euclidean');
eps = prctile(d, 10);
if eps == 0
    d_nz = d(d > 0);
    if isempty(d_nz)
        eps = eps + 1e-9;  % degenerate: all points identical
    else
        eps = prctile(d_nz, 10);
    end
end
end
