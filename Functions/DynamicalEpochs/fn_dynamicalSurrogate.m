function surr = fn_dynamicalSurrogate(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels, dyn_real)
% FN_DYNAMICALSURROGATE  Null control for fn_dynamicalEpochs.
% Re-runs the whole pipeline (jPCA -> cycles -> alignment -> clustering) on surrogate data in
% which every neuron is circularly shifted by an independent random amount. Rates and
% autocorrelation are kept; population coordination is destroyed. If the real data are not
% clearly beyond these surrogates, the "epochs" cannot be attributed to population structure.
%
%   surr = FN_DYNAMICALSURROGATE(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels, dyn_real)
%
% Metrics compared with the real data:
%   rot6      rotation strength (variance in the first jPC plane) in a fixed 6-PC subspace  (real > surrogate expected)
%   ndim80    PCs needed for 80% of the variance                                             (real < surrogate expected)
%   contrast  mean within-cluster minus between-cluster alignment                            (real > surrogate expected)
%   n_clusters, n_cycles, AUC are reported for reference (AUC is circular: clusters are built from the same matrix).
% p-values are one-sided, (1 + #surrogates as or more extreme than the data)/(1 + N), so with N
% surrogates the smallest attainable p is 1/(N+1).

n = cfg.DYN_NSURR;
nanv = nan(1, n);
surr = struct('n', n, 'rot6', nanv, 'ndim80', nanv, 'n_cycles', nanv, 'n_clusters', nanv, ...
    'auc', nanv, 'contrast', nanv, 'n_changes', nanv);
surr.real = struct('rot6', dyn_real.rot6, 'ndim80', dyn_real.ndim80, ...
    'n_cycles', numel(dyn_real.segment_s), 'n_clusters', dyn_real.nGroups, ...
    'auc', dyn_real.sim_thresh.auc, ...
    'contrast', dyn_real.sim_thresh.within_mean - dyn_real.sim_thresh.between_mean, ...
    'n_changes', numel(dyn_real.change_t));

fprintf('\n  Surrogate control (%d runs, each neuron circularly shifted independently):\n', n);
for i = 1:n
    c = cfg; c.DYN_SURR_SHIFT = true; c.DYN_SEED = cfg.DYN_SEED + i;
    d = struct('ok', false, 'message', '');
    try
        evalc('d = fn_dynamicalEpochs(spike_conv, fs, t_range_s, c, ep_bounds, ep_labels);');
    catch ME
        warning('fn_dynamicalSurrogate:Failed', 'Surrogate %d failed: %s', i, ME.message);
        continue
    end
    if isfield(d, 'rot6'),   surr.rot6(i)   = d.rot6;   end
    if isfield(d, 'ndim80'), surr.ndim80(i) = d.ndim80; end
    if ~d.ok
        fprintf('    Surrogate %d: %d PCs for 80%% variance, rotation %.1f%%, no cycles found (%s)\n', ...
            i, surr.ndim80(i), 100*surr.rot6(i), d.message);
        continue
    end
    surr.n_cycles(i)   = numel(d.segment_s);
    surr.n_clusters(i) = d.nGroups;
    surr.auc(i)        = d.sim_thresh.auc;
    surr.contrast(i)   = d.sim_thresh.within_mean - d.sim_thresh.between_mean;
    surr.n_changes(i)  = numel(d.change_t);
    fprintf('    Surrogate %d: %d PCs for 80%% variance, rotation %.1f%%, %d cycles, %d clusters, within-between %.2f, AUC %.2f\n', ...
        i, surr.ndim80(i), 100*surr.rot6(i), surr.n_cycles(i), surr.n_clusters(i), surr.contrast(i), surr.auc(i));
end

surr.n_failed = sum(isnan(surr.n_cycles));
surr.p = struct('rot6', fn_pval(surr.rot6, surr.real.rot6, 'ge'), ...
                'ndim80', fn_pval(surr.ndim80, surr.real.ndim80, 'le'), ...
                'contrast', fn_pval(surr.contrast, surr.real.contrast, 'ge'));
if surr.n_failed > 0, fprintf('  %d of %d surrogates yielded no cycles.\n', surr.n_failed, n); end
fprintf('  Data vs surrogate mean: rotation (6 PCs) %.1f%% vs %.1f%% (p=%.2f) | PCs for 80%% variance %d vs %.0f (p=%.2f) | within-between %.2f vs %.2f (p=%.2f)\n', ...
    100*surr.real.rot6, 100*mean(surr.rot6, 'omitnan'), surr.p.rot6, ...
    surr.real.ndim80, mean(surr.ndim80, 'omitnan'), surr.p.ndim80, ...
    surr.real.contrast, mean(surr.contrast, 'omitnan'), surr.p.contrast);
end

function p = fn_pval(v, r, dirn)
nv = sum(~isnan(v));
if nv == 0 || isnan(r)
    p = NaN;
elseif strcmp(dirn, 'ge')
    p = (1 + sum(v >= r)) / (1 + nv);
else
    p = (1 + sum(v <= r)) / (1 + nv);
end
end