function sweep = fn_dynamicalSweep(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels, dyn_main)
% FN_DYNAMICALSWEEP  Re-run fn_dynamicalEpochs at several recurrence thresholds
% (cfg.DYN_REC_SWEEP) to test whether the conclusions depend on that one setting.
%
%   sweep = FN_DYNAMICALSWEEP(spike_conv, fs, t_range_s, cfg, ep_bounds, ep_labels, dyn_main)
%
% dyn_main is the result already computed at cfg.DYN_REC_THRESH (reused, not recomputed).
% Prints a table and returns a struct of row vectors, one entry per threshold.

vals = cfg.DYN_REC_SWEEP(:)';
n = numel(vals);
nanv = nan(1, n);
sweep = struct('rec_thresh', vals, 'n_cycles', nanv, 'n_clusters', nanv, 'n_changes', nanv, ...
    'shared_E_in_R', nanv, 'shared_R_in_E', nanv, 'sim_ER', nanv, ...
    'sim_within_E', nanv, 'sim_within_R', nanv, 'same_dominant', nanv, ...
    'chance95', nanv, 'separation', nanv);

fprintf('\n  Recurrence-threshold sweep (%s %%):\n', mat2str(vals));
for k = 1:n
    if vals(k) == cfg.DYN_REC_THRESH && isfield(dyn_main, 'ok') && dyn_main.ok
        d = dyn_main;
    else
        c = cfg; c.DYN_REC_THRESH = vals(k);
        try
            d = fn_dynamicalEpochs(spike_conv, fs, t_range_s, c, ep_bounds, ep_labels);
        catch ME
            warning('fn_dynamicalSweep:Failed', 'Threshold %g%% failed: %s', vals(k), ME.message);
            continue
        end
    end
    if ~d.ok, continue; end
    sweep.n_cycles(k)   = numel(d.segment_s);
    sweep.n_clusters(k) = d.nGroups;
    sweep.n_changes(k)  = numel(d.change_t);
    sweep.chance95(k)   = d.sim_thresh.chance95;
    sweep.separation(k) = d.sim_thresh.separation;
    if isfield(d, 'EvR')
        sweep.shared_E_in_R(k) = d.EvR.shared_frac_E_in_R;
        sweep.shared_R_in_E(k) = d.EvR.shared_frac_R_in_E;
        sweep.sim_ER(k)        = d.EvR.sim_between;
        sweep.sim_within_E(k)  = d.EvR.sim_within_evoked;
        sweep.sim_within_R(k)  = d.EvR.sim_within_recovery;
        sweep.same_dominant(k) = d.EvR.same_dominant_cluster;
    end
end

fprintf('  thresh | cycles | clusters | changes | E in R | R in E | E<->R | within E | within R | same dominant\n');
for k = 1:n
    fprintf('  %5g%% | %6g | %8g | %7g | %6.2f | %6.2f | %5.2f | %8.2f | %8.2f | %g\n', vals(k), ...
        sweep.n_cycles(k), sweep.n_clusters(k), sweep.n_changes(k), sweep.shared_E_in_R(k), ...
        sweep.shared_R_in_E(k), sweep.sim_ER(k), sweep.sim_within_E(k), sweep.sim_within_R(k), sweep.same_dominant(k));
end
end
