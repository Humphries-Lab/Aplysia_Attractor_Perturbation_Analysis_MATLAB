function dyn = fn_dynamicalEpochs(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels)
% FN_DYNAMICALEPOCHS  Data-driven epoch detection by cycle-by-cycle subspace
% alignment (Colins Rodriguez's jPCA / recurrence / clustering pipeline).
%
%   dyn = FN_DYNAMICALEPOCHS(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels)
%
% WHAT IT DOES (window = t_range_s of spike_conv)
%   1. Bins the smoothed activity to cfg.DYN_BIN_MS (time x neurons).
%   2. jPCA isolates the rotational component, and recurrence in the jPCA
%      plane closes each cycle (detect_cycle_JPCA_v2)  -> cycle boundaries.
%   3. Subspace alignment between EVERY pair of cycles (overlap_similarity).
%   4. Clusters the cycle x cycle alignment matrix (clustering) -> groups
%      of cycles that share a subspace = dynamical epochs.
%   5. Maps each cycle onto the protocol epochs (ep_bounds/ep_labels) so you
%      can see whether the same cluster appears either side of the
%      perturbation (e.g. Evoked vs Recovery).
%
% INPUTS
%   spike_conv - neurons x frames smoothed activity (as in the main pipeline)
%   fs         - frames per second of spike_conv
%   t_range_s  - [t_start t_end] (s) of the window to analyse
%   cfg        - needs DYN_BIN_MS, DYN_MIN_SEG_S, DYN_REC_THRESH, DYN_SEED,
%                DYN_CODE_DIR (see Config_UserSettings.m)
%   ep_bounds  - nE x 2 [start end] (s), protocol epochs used only to LABEL cycles
%   ep_labels  - 1 x nE cell of epoch names
%
% OUTPUT dyn (struct; dyn.ok = false if no cycles were found)
%   .segment_s, .startbin, .ndim_cycle  - cycle durations (s), start bins, jPCA dims
%   .cycle_t [nCyc x 2], .cycle_mid     - cycle start/end and midpoint (s, recording time)
%   .similarity      - nCyc x nCyc alignment (row i = cycle i's variance captured by cycle j's subspace)
%   .grps, .nGroups  - cluster label per cycle, number of clusters
%   .epoch_of_cycle  - index into ep_labels by cycle midpoint (0 = outside all)
%   .cluster_by_epoch- nGroups x nE cycle counts
%   .epoch_sim       - nE x nE mean (symmetrised) alignment between cycles of two epochs
%   .shared_frac     - (a,b) fraction of epoch-a cycles whose cluster also occurs in epoch b
%   .dominant_cluster- most frequent cluster per epoch
%   .EvR             - Evoked vs Recovery summary (only if both labels exist and have cycles)
%   .jpca_var_top2, .is_rotational - rotational variance (first jPC plane) of the whole window
%   .runs, .change_t, .change_from/.change_to, .returns_to - contiguous stretches of one cluster
%                    (dynamical epochs), the times where the epoch changes, clusters that recur
%   .landmarks       - protocol landmarks and the signed lag (s) to the nearest detected change
%   .sim_thresh      - within/between-cluster alignment, shuffle-chance (95th pct), and the
%                    threshold that best separates same- from different-cluster pairs (+AUC)
%   .landmark_test   - permutation test: are change points closer to the landmarks than chance?
%   .boundary        - contrast (within - between block alignment) for every split; landmark percentiles
%   .BvR             - Baseline vs Recovery summary
%   .drift           - lag-matched check: Evoked<->Recovery alignment vs what within-epoch drift predicts
%
% Requires Statistics Toolbox (pca, pdist) and the Dynamical_epochs +
% jPCA_Aplysia code in cfg.DYN_CODE_DIR (added to the path at the END so it
% never shadows the toolbox's own functions).

fn_addDynamicalEpochsPath(cfg.DYN_CODE_DIR);

dyn = struct('ok', false, 'message', '', 'jpca_error', '', 't_range', t_range_s, ...
    'ep_bounds', ep_bounds, 'ep_labels', {ep_labels});

% keep the user's global state untouched
ws = warning('off', 'stats:pca:ColRankDefX');
c1 = onCleanup(@() warning(ws));
rngState = rng; %#ok<NASGU>
c2 = onCleanup(@() rng(rngState));
rng(cfg.DYN_SEED);

%% 1. Bin to time x neurons at cfg.DYN_BIN_MS
[nN, nFrames] = size(spike_conv);
bf      = max(round(fs * cfg.DYN_BIN_MS / 1000), 1);   % frames per bin
bin_eff = 1000 * bf / fs;                              % actual bin (ms)
f1 = max(round(t_range_s(1) * fs) + 1, 1);
f2 = min(round(t_range_s(2) * fs), nFrames);
nb = floor((f2 - f1 + 1) / bf);
t0 = (f1 - 1) / fs;

if nb * bin_eff / 1000 < 2 * cfg.DYN_MIN_SEG_S
    dyn.message = 'Analysis window too short for cycle detection.';
    warning('fn_dynamicalEpochs:ShortWindow', '%s', dyn.message);
    return
end

X = zeros(nb, nN);
chunk = 500;   % bins per chunk, keeps memory low
for b0 = 1:chunk:nb
    b1  = min(b0 + chunk - 1, nb);
    seg = double(spike_conv(:, f1 + (b0-1)*bf : f1 + b1*bf - 1));
    X(b0:b1, :) = reshape(mean(reshape(seg, nN, bf, b1-b0+1), 2), nN, b1-b0+1)';
end
if isfield(cfg, 'DYN_SURR_SHIFT') && cfg.DYN_SURR_SHIFT
    % surrogate control: circularly shift each neuron independently (keeps rates and
    % autocorrelation, destroys the population coordination)
    for nn = 1:nN, X(:, nn) = circshift(X(:, nn), randi(nb)); end
end
fprintf('  Dynamical epochs: %d neurons, %d bins of %.1f ms (%.0f-%.0f s)\n', ...
    nN, nb, bin_eff, t0, t0 + nb*bin_eff/1000);

% jPCA calls keyboard (drops into the debugger) when its input-size checks fail, e.g. when
% no low-dimensional rotational structure exists and it keeps adding PCs. Turn that into an
% error, only while this function runs.
c3 = fn_installKeyboardGuard(); %#ok<NASGU>

%% 2. Rotational component + cycles (jPCA + recurrence)
% dimensionality: PCs needed for 80% of the variance (is_rotational uses this as its jPCA dimension)
[~, ~, ~, ~, expl80] = pca(X);
dyn.ndim80 = find(cumsum(expl80) >= 80, 1, 'first');
fprintf('  PCs needed for 80%% of variance: %d\n', dyn.ndim80);

% rotation strength in a FIXED 6-PC subspace, so real data and surrogates are directly comparable
dyn.rot6 = NaN;
try
    p6.params = false; p6.numPCs = 6; p6.suppressBWrosettes = true;
    p6.suppressHistograms = true; p6.suppressText = true;
    D6(1).A = X; D6(1).times = ((1:nb) * 50)';
    [~, S6] = jPCA(D6, (1:nb) * 50, p6);
    dyn.rot6 = sum(S6.varCaptEachJPC(1:2));
    fprintf('  jPCA (fixed 6 PCs): first jPC plane captures %.1f%% of variance\n', 100 * dyn.rot6);
catch ME6
    warning('fn_dynamicalEpochs:Rot6Failed', 'Fixed 6-PC jPCA failed: %s', ME6.message);
end

try
    if dyn.ndim80 > 20   % Andrea's jPCA code refuses more than 20 PCs (would stop in the debugger)
        error('needs %d PCs for 80%% of the variance (>20): no low-dimensional structure', dyn.ndim80);
    end
    [~, Summ, ~] = is_rotational(X, 0);
    dyn.jpca_var_each = Summ.varCaptEachJPC;
    dyn.jpca_var_top2 = sum(Summ.varCaptEachJPC(1:2));
    dyn.is_rotational = dyn.jpca_var_top2 >= 0.1;
    fprintf('  jPCA: first jPC plane captures %.1f%% of variance\n', 100*dyn.jpca_var_top2);
catch ME
    warning('fn_dynamicalEpochs:jPCAFailed', 'is_rotational failed: %s', ME.message);
    dyn.jpca_var_top2 = NaN; dyn.is_rotational = NaN; dyn.jpca_error = ME.message;
end

try
    [segment_s, startbin, ndim_cycle] = detect_cycle_JPCA_v2(X, cfg.DYN_MIN_SEG_S, bin_eff, cfg.DYN_REC_THRESH, 0);
catch ME
    dyn.message = sprintf('Cycle detection stopped (%s) - typically no low-dimensional rotational structure.', ME.message);
    warning('fn_dynamicalEpochs:CycleDetectionFailed', '%s', dyn.message);
    return
end
nSeg = numel(segment_s);
fprintf('  Cycles detected: %d (mean %.1f s)\n', nSeg, mean(segment_s));
if nSeg < 1
    dyn.message = 'No recurring cycles found (no rotational structure at this recurrence threshold).';
    warning('fn_dynamicalEpochs:NoCycles', '%s', dyn.message);
    return
end

% Cycles are only closed where the trajectory recurs. A stretch without rotational
% structure (e.g. quiet baseline) therefore shows up as one abnormally long "cycle".
long_cycle = segment_s > 2 * median(segment_s);
if any(long_cycle)
    fprintf('  Note: %d cycle(s) are >2x the median length (%s s) - likely non-rotational stretches, not true cycles.\n', ...
        sum(long_cycle), mat2str(round(segment_s(long_cycle)', 0)));
end

%% 3. Cycle-by-cycle subspace alignment
[similarity, coeff_1, ndim_exp] = overlap_similarity(X, startbin, 0, 0);

%% 4. Cluster the alignment matrix
grps = ones(nSeg, 1);   % default: one cluster (too few cycles to cluster)
if nSeg >= 3
    try
        g = clustering(similarity);
        grps = g(:, 1);
        if max(grps) == 0, grps = ones(nSeg, 1); end
        grps(grps == 0) = max(grps) + 1;
    catch ME
        warning('fn_dynamicalEpochs:ClusterFailed', 'Clustering failed (%s) -- using a single cluster.', ME.message);
        grps = ones(nSeg, 1);
    end
end
nGroups = max(grps);
fprintf('  Clusters (dynamical epochs) found: %d\n', nGroups);

%% 5. Map cycles onto protocol epochs (recording time)
cycle_t   = t0 + ([startbin(1:end-1), startbin(2:end)] - 1) * bin_eff / 1000;
cycle_mid = mean(cycle_t, 2);
nE = numel(ep_labels);

epoch_of_cycle = zeros(nSeg, 1);
for e = 1:nE
    epoch_of_cycle(cycle_mid >= ep_bounds(e,1) & cycle_mid < ep_bounds(e,2)) = e;
end

cluster_by_epoch = zeros(nGroups, nE);
for c = 1:nSeg
    if epoch_of_cycle(c) > 0
        cluster_by_epoch(grps(c), epoch_of_cycle(c)) = cluster_by_epoch(grps(c), epoch_of_cycle(c)) + 1;
    end
end

Ssym = (similarity + similarity') / 2;
epoch_sim   = nan(nE, nE);
shared_frac = nan(nE, nE);
for a = 1:nE
    ia = find(epoch_of_cycle == a);
    for b = 1:nE
        ib = find(epoch_of_cycle == b);
        if isempty(ia) || isempty(ib), continue; end
        M = Ssym(ia, ib);
        if a == b
            if numel(ia) < 2, continue; end
            M(logical(eye(numel(ia)))) = NaN;   % drop self-alignment
        end
        epoch_sim(a, b)   = mean(M(:), 'omitnan');
        shared_frac(a, b) = sum(cluster_by_epoch(:, a) .* (cluster_by_epoch(:, b) > 0)) / sum(cluster_by_epoch(:, a));
    end
end

dominant_cluster = nan(1, nE);
[~, dmx] = max(cluster_by_epoch, [], 1);
has = sum(cluster_by_epoch, 1) > 0;
dominant_cluster(has) = dmx(has);

%% Evoked vs Recovery: is the attractor the "same" either side of the perturbation?
iE = find(strcmp(ep_labels, 'Evoked'), 1);
iR = find(strcmp(ep_labels, 'Recovery'), 1);
if ~isempty(iE) && ~isempty(iR) && has(iE) && has(iR)
    EvR.n_cycles_evoked       = sum(epoch_of_cycle == iE);
    EvR.n_cycles_recovery     = sum(epoch_of_cycle == iR);
    EvR.shared_frac_E_in_R    = shared_frac(iE, iR);
    EvR.shared_frac_R_in_E    = shared_frac(iR, iE);
    EvR.sim_between           = epoch_sim(iE, iR);
    EvR.sim_within_evoked     = epoch_sim(iE, iE);
    EvR.sim_within_recovery   = epoch_sim(iR, iR);
    EvR.same_dominant_cluster = dominant_cluster(iE) == dominant_cluster(iR);
    fprintf('  Evoked vs Recovery: %d vs %d cycles | shared cluster (E in R) %.2f, (R in E) %.2f\n', ...
        EvR.n_cycles_evoked, EvR.n_cycles_recovery, EvR.shared_frac_E_in_R, EvR.shared_frac_R_in_E);
    fprintf('    alignment E<->R %.2f | within E %.2f | within R %.2f | same dominant cluster: %d\n', ...
        EvR.sim_between, EvR.sim_within_evoked, EvR.sim_within_recovery, EvR.same_dominant_cluster);
    dyn.EvR = EvR;
end


%% 6. Change points: where does the dynamical epoch change along time?
%  (unsupervised: a change point is a cycle whose cluster differs from the previous cycle)
chg         = find(diff(grps) ~= 0) + 1;          % first cycle of each new run
run_first   = [1; chg(:)];
run_last    = [chg(:) - 1; nSeg];
run_cluster = grps(run_first);
run_t       = [cycle_t(run_first, 1), cycle_t(run_last, 2)];
run_ncycles = run_last - run_first + 1;
change_t    = cycle_t(chg, 1);                     % time (s) of each change of dynamical epoch
change_from = grps(chg - 1);
change_to   = grps(chg);
runs_per_cluster = accumarray(grps, 1, [nGroups 1]);
returns_to  = find(runs_per_cluster > 1)';         % clusters that come back after a different one

% distance from each protocol landmark to the nearest detected change point
lm_t    = [ep_bounds(1, 2); ep_bounds(2:end, 1)];      % P9, then the start of each later epoch
lm_name = [{'P9'}, cellfun(@(x) ['start of ' x], ep_labels(2:end), 'UniformOutput', false)];
lm_short = [{'P9'}, cellfun(@(x) x(1:min(3, numel(x))), ep_labels(2:end), 'UniformOutput', false)];
lm_lag  = nan(numel(lm_t), 1);                     % signed: change point minus landmark (s)
for k = 1:numel(lm_t)
    if lm_t(k) > t_range_s(1) && lm_t(k) < t_range_s(2) && ~isempty(change_t)
        [~, ix] = min(abs(change_t - lm_t(k)));
        lm_lag(k) = change_t(ix) - lm_t(k);
    end
end
fprintf('  Change points (s): %s\n', mat2str(round(change_t', 1)));
for k = 1:numel(lm_t)
    if ~isnan(lm_lag(k))
        fprintf('    nearest change to %s (%.0f s): %+.1f s (mean cycle %.1f s)\n', ...
            lm_name{k}, lm_t(k), lm_lag(k), mean(segment_s));
    end
end


% Is the proximity of change points to the protocol landmarks better than chance?
% Null: the same number of change points placed on randomly chosen cycle boundaries.
nPerm = 10000;
if isfield(cfg, 'DYN_NPERM'), nPerm = cfg.DYN_NPERM; end
lm_in = find(~isnan(lm_lag));
landmark_test = struct('ok', false, 'observed', NaN, 'p', NaN, 'p_each', nan(numel(lm_t), 1), 'null', []);
cand = cycle_t(2:end, 1);                          % change points can only fall on cycle boundaries
nChg = numel(change_t);
if ~isempty(lm_in) && nChg >= 1 && numel(cand) >= nChg
    obsAbs = abs(lm_lag(lm_in));
    M = nan(nPerm, numel(lm_in));
    for q = 1:nPerm
        ct = cand(randperm(numel(cand), nChg));
        M(q, :) = min(abs(ct(:) - lm_t(lm_in)'), [], 1);
    end
    nullMean = mean(M, 2);
    landmark_test.ok = true;
    landmark_test.observed = mean(obsAbs);
    landmark_test.null = nullMean;
    landmark_test.p = (1 + sum(nullMean <= mean(obsAbs))) / (1 + nPerm);
    landmark_test.p_each(lm_in) = (1 + sum(M <= obsAbs', 1)') / (1 + nPerm);
    landmark_test.expected_random = mean(nullMean);
    fprintf('  Change points vs landmarks: mean |lag| %.1f s (random expectation %.1f s), permutation p = %.3f\n', ...
        landmark_test.observed, landmark_test.expected_random, landmark_test.p);
end

%% 7. How similar is "the same"? within- vs between-cluster alignment
iu = triu(true(nSeg), 1);
[pi_, pj_] = find(iu);
sv   = Ssym(iu);                                   % same (column-major) order as find(iu)
same = grps(pi_) == grps(pj_);
w = sv(same); b = sv(~same);

% empirical chance: alignment after shuffling neuron identity of the target subspace
nNull = 2000;
if isfield(cfg, 'DYN_NNULL'), nNull = cfg.DYN_NNULL; end
nullv = []; chance95 = NaN;
if nSeg >= 2
    Cs = cell(nSeg, 1); sing = cell(nSeg, 1);
    for i = 1:nSeg
        Cs{i}   = cov(X(startbin(i):startbin(i+1)-1, :));
        sing{i} = flipud(eig(Cs{i}));
    end
    nullv = nan(nNull, 1);
    for sIdx = 1:nNull
        i = randi(nSeg); j = randi(nSeg - 1); if j >= i, j = j + 1; end
        k = ndim_exp(j);
        D = coeff_1(randperm(nN), 1:k, j);
        nullv(sIdx) = shared_variance_v2(Cs{i}, D, sing{i}, k);
    end
    chance95 = prctile(nullv, 95);
end

% threshold that best separates same-cluster from different-cluster pairs
thr_sep = NaN; auc = NaN; bal_acc = NaN;
if ~isempty(w) && ~isempty(b)
    best = -Inf;
    for t = unique(sv)'
        J = mean(w >= t) + mean(b < t) - 1;
        if J > best, best = J; thr_sep = t; end
    end
    bal_acc = (best + 1) / 2;
    r   = tiedrank([w; b]);
    auc = (sum(r(1:numel(w))) - numel(w)*(numel(w)+1)/2) / (numel(w)*numel(b));
end
sim_thresh = struct('within_mean', fn_nanmean(w), 'within_min', fn_nanmin(w), ...
    'between_mean', fn_nanmean(b), 'between_max', fn_nanmax(b), ...
    'chance95', chance95, 'separation', thr_sep, 'auc', auc, 'balanced_acc', bal_acc, ...
    'n_within', numel(w), 'n_between', numel(b), 'null', nullv);
fprintf('  Similarity: within-cluster %.2f | between-cluster %.2f | chance (95th pct) %.2f | separating threshold %.2f (AUC %.2f)\n', ...
    sim_thresh.within_mean, sim_thresh.between_mean, chance95, thr_sep, auc);

%% 8. Is a lower Evoked<->Recovery alignment more than slow drift? (lag-matched check)
lag_s = abs(cycle_mid(pj_) - cycle_mid(pi_));
ei = epoch_of_cycle(pi_); ej = epoch_of_cycle(pj_);
ptype = zeros(numel(sv), 1);
ptype(ei > 0 & ei == ej) = 1;                                        % both cycles in the same protocol epoch
if ~isempty(iE) && ~isempty(iR)
    ptype((ei == iE & ej == iR) | (ei == iR & ej == iE)) = 2;        % Evoked <-> Recovery
end
drift = struct('ok', false);
m1 = ptype == 1; m2 = ptype == 2;
if sum(m1) >= 3 && any(m2) && numel(unique(lag_s(m1))) >= 2
    p = polyfit(lag_s(m1), sv(m1), 1);
    pred = polyval(p, lag_s(m2));
    drift.ok = true;
    drift.slope_per_min = p(1) * 60;
    drift.intercept = p(2);
    drift.obs_ER  = mean(sv(m2));
    drift.pred_ER = mean(pred);
    drift.resid_ER = drift.obs_ER - drift.pred_ER;
    drift.frac_extrapolated = mean(lag_s(m2) > max(lag_s(m1)));
    drift.n_within = sum(m1); drift.n_ER = sum(m2);
    fprintf('  Drift check: within-epoch alignment changes by %+.3f per min; E<->R observed %.2f vs %.2f expected from drift alone (%.0f%% extrapolated)\n', ...
        drift.slope_per_min, drift.obs_ER, drift.pred_ER, 100*drift.frac_extrapolated);
end


%% 9. Is a protocol boundary special? Segmentation contrast for every possible split
%  contrast(k) = mean alignment within the two blocks (cycles 1..k | k+1..end) minus mean alignment between them
minBlock = 3;
bc_idx = minBlock:(nSeg - minBlock);
boundary = struct('ok', false, 'idx', bc_idx, 't', [], 'contrast', [], 'best_t', NaN, 'best_contrast', NaN, ...
    'lm_contrast', nan(numel(lm_t), 1), 'lm_pctile', nan(numel(lm_t), 1));
if numel(bc_idx) >= 2
    Sd = Ssym; Sd(logical(eye(nSeg))) = NaN;
    ctr = nan(1, numel(bc_idx));
    for q = 1:numel(bc_idx)
        k = bc_idx(q);
        A = Sd(1:k, 1:k); B = Sd(k+1:end, k+1:end); Cx = Sd(1:k, k+1:end);
        ctr(q) = mean([A(~isnan(A)); B(~isnan(B))]) - mean(Cx(:));
    end
    boundary.ok = true;
    boundary.t = cycle_t(bc_idx, 2)';
    boundary.contrast = ctr;
    [boundary.best_contrast, ib] = max(ctr);
    boundary.best_t = boundary.t(ib);
    for k = 1:numel(lm_t)
        kk = sum(cycle_mid < lm_t(k));
        q = find(bc_idx == kk, 1);
        if ~isempty(q)
            boundary.lm_contrast(k) = ctr(q);
            boundary.lm_pctile(k)   = mean(ctr <= ctr(q));
        end
    end
    fprintf('  Strongest single split of the cycle sequence: t = %.0f s (contrast %.2f)\n', boundary.best_t, boundary.best_contrast);
    for k = 1:numel(lm_t)
        if ~isnan(boundary.lm_pctile(k))
            fprintf('    %s (%.0f s): contrast %.2f, percentile among all splits %.0f\n', ...
                lm_name{k}, lm_t(k), boundary.lm_contrast(k), 100 * boundary.lm_pctile(k));
        end
    end
end

%% 10. Full epoch-by-epoch alignment (incl. Baseline) and Baseline vs Recovery
fprintf('  Mean alignment between protocol epochs:\n');
fprintf('    %-10s', ''); fprintf(' %9s', ep_labels{:}); fprintf('\n');
for a = 1:nE
    fprintf('    %-10s', ep_labels{a}); fprintf(' %9.2f', epoch_sim(a, :)); fprintf('\n');
end
iB = find(strcmp(ep_labels, 'Baseline'), 1);
BvR = struct('ok', false);
if ~isempty(iB) && ~isempty(iR) && has(iB) && has(iR)
    BvR.ok = true;
    BvR.n_cycles_baseline = sum(epoch_of_cycle == iB);
    BvR.n_cycles_recovery = sum(epoch_of_cycle == iR);
    BvR.sim_between = epoch_sim(iB, iR);
    BvR.sim_within_baseline = epoch_sim(iB, iB);
    BvR.shared_frac_B_in_R = shared_frac(iB, iR);
    BvR.shared_frac_R_in_B = shared_frac(iR, iB);
    BvR.same_dominant_cluster = dominant_cluster(iB) == dominant_cluster(iR);
    fprintf('  Baseline vs Recovery: %d vs %d cycles | alignment %.2f (within Baseline %.2f) | shared cluster (B in R) %.2f, (R in B) %.2f | same dominant cluster: %d\n', ...
        BvR.n_cycles_baseline, BvR.n_cycles_recovery, BvR.sim_between, BvR.sim_within_baseline, ...
        BvR.shared_frac_B_in_R, BvR.shared_frac_R_in_B, BvR.same_dominant_cluster);
end

%% Output
dyn.ok = true;
dyn.bin_ms = bin_eff; dyn.t0 = t0;
dyn.long_cycle = long_cycle; dyn.segment_s = segment_s; dyn.startbin = startbin; dyn.ndim_cycle = ndim_cycle;
dyn.cycle_t = cycle_t; dyn.cycle_mid = cycle_mid;
dyn.similarity = similarity; dyn.grps = grps; dyn.nGroups = nGroups;
dyn.epoch_of_cycle = epoch_of_cycle; dyn.cluster_by_epoch = cluster_by_epoch;
dyn.epoch_sim = epoch_sim; dyn.shared_frac = shared_frac; dyn.dominant_cluster = dominant_cluster;
dyn.runs = struct('cluster', run_cluster, 't', run_t, 'ncycles', run_ncycles);
dyn.change_t = change_t; dyn.change_from = change_from; dyn.change_to = change_to;
dyn.returns_to = returns_to;
dyn.landmarks = struct('t', lm_t, 'name', {lm_name}, 'short', {lm_short}, 'lag_to_nearest_change_s', lm_lag);
dyn.landmark_test = landmark_test;
dyn.boundary = boundary;
dyn.BvR = BvR;
dyn.sim_thresh = sim_thresh;
dyn.drift = drift;
dyn.pairs = struct('lag_s', lag_s, 'sim', sv, 'type', ptype, 'same_cluster', same);
dyn.params = struct('bin_ms', bin_eff, 'min_seg_s', cfg.DYN_MIN_SEG_S, ...
    'rec_thresh', cfg.DYN_REC_THRESH, 'seed', cfg.DYN_SEED);
end

function fn_addDynamicalEpochsPath(codeDir)
needed = {'detect_cycle_JPCA_v2','overlap_similarity','clustering','is_rotational','jPCA_new','jPCA'};
ok = @() cellfun(@(f) exist(f, 'file') == 2, needed);
if ~all(ok()) && isfolder(codeDir)
    addpath(genpath(codeDir), '-end');
end
if ~all(ok())
    error('fn_dynamicalEpochs:MissingCode', ...
        ['Dynamical_epochs code not found (missing: %s).\n' ...
         'Clone it, with its jPCA_Aplysia submodule, into:\n  %s\n' ...
         '  git clone --recurse-submodules https://github.com/AndreaColinsR/Dynamical_epochs'], ...
        strjoin(needed(~ok()), ', '), codeDir);
end
end

function v = fn_nanmean(x), if isempty(x), v = NaN; else, v = mean(x); end, end
function v = fn_nanmin(x),  if isempty(x), v = NaN; else, v = min(x);  end, end
function v = fn_nanmax(x),  if isempty(x), v = NaN; else, v = max(x);  end, end

function c = fn_installKeyboardGuard()
guardDir = tempname; mkdir(guardDir);
fid = fopen(fullfile(guardDir, 'keyboard.m'), 'w');
fprintf(fid, 'function keyboard(varargin)\nerror(''fn_dynamicalEpochs:jPCAdebugStop'', ''jPCA reached a debugger stop (input-size check).'');\nend\n');
fclose(fid);
w = warning('off', 'MATLAB:dispatcher:nameConflict');
addpath(guardDir);
c = onCleanup(@() fn_removeKeyboardGuard(guardDir, w));
end

function fn_removeKeyboardGuard(guardDir, w)
rmpath(guardDir);
if isfolder(guardDir), rmdir(guardDir, 's'); end
warning(w);
end
