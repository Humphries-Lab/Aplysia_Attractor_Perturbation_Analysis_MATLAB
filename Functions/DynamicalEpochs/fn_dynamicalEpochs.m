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
%
% Requires Statistics Toolbox (pca, pdist) and the Dynamical_epochs +
% jPCA_Aplysia code in cfg.DYN_CODE_DIR (added to the path at the END so it
% never shadows the toolbox's own functions).

fn_addDynamicalEpochsPath(cfg.DYN_CODE_DIR);

dyn = struct('ok', false, 'message', '', 't_range', t_range_s, ...
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
fprintf('  Dynamical epochs: %d neurons, %d bins of %.1f ms (%.0f-%.0f s)\n', ...
    nN, nb, bin_eff, t0, t0 + nb*bin_eff/1000);

%% 2. Rotational component + cycles (jPCA + recurrence)
try
    [~, Summ, ~] = is_rotational(X, 0);
    dyn.jpca_var_each = Summ.varCaptEachJPC;
    dyn.jpca_var_top2 = sum(Summ.varCaptEachJPC(1:2));
    dyn.is_rotational = dyn.jpca_var_top2 >= 0.1;
    fprintf('  jPCA: first jPC plane captures %.1f%% of variance\n', 100*dyn.jpca_var_top2);
catch ME
    warning('fn_dynamicalEpochs:jPCAFailed', 'is_rotational failed: %s', ME.message);
    dyn.jpca_var_top2 = NaN; dyn.is_rotational = NaN;
end

[segment_s, startbin, ndim_cycle] = detect_cycle_JPCA_v2(X, cfg.DYN_MIN_SEG_S, bin_eff, cfg.DYN_REC_THRESH, 0);
nSeg = numel(segment_s);
fprintf('  Cycles detected: %d (mean %.1f s)\n', nSeg, mean(segment_s));
if nSeg < 1
    dyn.message = 'No recurring cycles found (no rotational structure at this recurrence threshold).';
    warning('fn_dynamicalEpochs:NoCycles', '%s', dyn.message);
    return
end

%% 3. Cycle-by-cycle subspace alignment
similarity = overlap_similarity(X, startbin, 0, 0);

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

%% Output
dyn.ok = true;
dyn.bin_ms = bin_eff; dyn.t0 = t0;
dyn.segment_s = segment_s; dyn.startbin = startbin; dyn.ndim_cycle = ndim_cycle;
dyn.cycle_t = cycle_t; dyn.cycle_mid = cycle_mid;
dyn.similarity = similarity; dyn.grps = grps; dyn.nGroups = nGroups;
dyn.epoch_of_cycle = epoch_of_cycle; dyn.cluster_by_epoch = cluster_by_epoch;
dyn.epoch_sim = epoch_sim; dyn.shared_frac = shared_frac; dyn.dominant_cluster = dominant_cluster;
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
