function cyc = fn_dynamicalCycleMetrics(spike_conv, fs, dyn, stim_t, stim_names, t_ref_start)
% FN_DYNAMICALCYCLEMETRICS  Per-cycle readouts on the cycles found by fn_dynamicalEpochs.
%
%   cyc = FN_DYNAMICALCYCLEMETRICS(spike_conv, fs, dyn, stim_t, stim_names, t_ref_start)
%
%   1. PR per cycle          participation ratio of the population covariance within each
%                            cycle (PR/N) -> time-series of dimensionality.
%   2. Alignment per cycle   (a) cycle c vs c-1 (step-to-step) and (b) cycle c vs a reference
%                            cycle -> time-series of movement of the attractor. Uses the
%                            cycle x cycle alignment already computed in dyn.similarity
%                            (symmetrised). Movement = 1 - alignment.
%   3. Stimulus phase        for each stimulus time, the cycle that contains it and the
%                            fraction (0-1) and angle (rad) of the cycle elapsed at that time.
%
% INPUTS
%   spike_conv  - neurons x frames smoothed activity (same as main pipeline)
%   fs          - frames/s
%   dyn         - output of fn_dynamicalEpochs (needs dyn.ok, cycle_t, similarity, grps)
%   stim_t      - stimulus times (s), e.g. [t_P9 t_C2]
%   stim_names  - cell array of names, e.g. {'P9','C2'}
%   t_ref_start - (optional) reference cycle = first cycle with midpoint >= t_ref_start
%                 (default: cycle 1). Suggest t_evoked_start.
%
% OUTPUT cyc (struct)
%   .cycle_mid, .cycle_dur_s, .grps
%   .PR (raw), .PR_norm (PR/N)
%   .align_prev (nCyc x 1, NaN for cycle 1), .align_ref, .ref_cycle
%   .move_prev = 1 - align_prev, .move_ref = 1 - align_ref
%   .stim: struct array with name, t, cycle, phase_frac, phase_rad, in_long_cycle, PR_norm,
%          align_prev (this cycle vs previous), align_next (next cycle vs this)
%
% NOTES
%   * Alignment is asymmetric in dyn.similarity (row i = variance of cycle i captured by the
%     subspace of cycle j); it is averaged with its transpose here.
%   * Phase is the fraction of TIME elapsed within the recurrence-defined cycle. Cycles are cut
%     end-to-start, so phase 0 is where the cycle was cut, not a fixed point on the attractor:
%     phases are only comparable across cycles to the extent the cut points sit at similar states.
%   * A cycle longer than 2x the median is flagged (likely a non-rotational stretch).

cyc = struct('ok', false);
if nargin < 6 || isempty(t_ref_start), t_ref_start = -inf; end
if ~isfield(dyn, 'ok') || ~dyn.ok
    fprintf('  Cycle metrics: skipped (no cycles).\n'); return
end
if nargin < 5 || isempty(stim_names), stim_names = arrayfun(@(k) sprintf('S%d', k), 1:numel(stim_t), 'UniformOutput', false); end

[nN, nFrames] = size(spike_conv);
nCyc = size(dyn.cycle_t, 1);

%% 1. PR per cycle
PR = nan(nCyc, 1);
for c = 1:nCyc
    f1 = max(round(dyn.cycle_t(c,1) * fs) + 1, 1);
    f2 = min(round(dyn.cycle_t(c,2) * fs), nFrames);
    if f2 - f1 < 10, continue; end
    Xc = double(spike_conv(:, f1:f2));
    Xc = Xc - mean(Xc, 2);
    ev = eig(cov(Xc'));                 % neuron x neuron covariance within the cycle
    PR(c) = fn_participationRatio(real(ev));
end

%% 2. Alignment per cycle (from dyn.similarity)
S = (dyn.similarity + dyn.similarity') / 2;
align_prev = nan(nCyc, 1);
for c = 2:nCyc, align_prev(c) = S(c, c-1); end

ref = find(dyn.cycle_mid >= t_ref_start, 1, 'first');
if isempty(ref), ref = 1; end
align_ref = S(:, ref);
align_ref(ref) = NaN;                   % self-alignment is 1 by construction

%% 3. Phase of each stimulus within its cycle
isLong = dyn.cycle_t(:,2) - dyn.cycle_t(:,1) > 2 * median(dyn.cycle_t(:,2) - dyn.cycle_t(:,1));
stim = struct('name', {}, 't', {}, 'cycle', {}, 'phase_frac', {}, 'phase_rad', {}, ...
    'in_long_cycle', {}, 'PR_norm', {}, 'align_prev', {}, 'align_next', {});
for k = 1:numel(stim_t)
    t = stim_t(k);
    c = find(dyn.cycle_t(:,1) <= t & t < dyn.cycle_t(:,2), 1, 'first');
    s.name = stim_names{k}; s.t = t; s.cycle = NaN;
    s.phase_frac = NaN; s.phase_rad = NaN; s.in_long_cycle = false;
    s.PR_norm = NaN; s.align_prev = NaN; s.align_next = NaN;
    if ~isempty(c)
        s.cycle      = c;
        s.phase_frac = (t - dyn.cycle_t(c,1)) / (dyn.cycle_t(c,2) - dyn.cycle_t(c,1));
        s.phase_rad  = 2*pi*s.phase_frac;
        s.in_long_cycle = isLong(c);
        s.PR_norm    = PR(c) / nN;
        s.align_prev = align_prev(c);
        if c < nCyc, s.align_next = align_prev(c+1); end
    end
    stim(end+1) = s; %#ok<AGROW>
end

%% Output + printout
cyc.ok = true;
cyc.cycle_mid = dyn.cycle_mid; cyc.cycle_t = dyn.cycle_t;
cyc.cycle_dur_s = dyn.cycle_t(:,2) - dyn.cycle_t(:,1); cyc.grps = dyn.grps;
cyc.PR = PR; cyc.PR_norm = PR / nN;
cyc.align_prev = align_prev; cyc.align_ref = align_ref; cyc.ref_cycle = ref;
cyc.move_prev = 1 - align_prev; cyc.move_ref = 1 - align_ref;
cyc.long_cycle = isLong; cyc.stim = stim;

fprintf('  Cycle metrics: %d cycles | PR/N %.2f-%.2f | step alignment %.2f-%.2f | reference = cycle %d\n', ...
    nCyc, min(cyc.PR_norm), max(cyc.PR_norm), min(align_prev), max(align_prev), ref);
for k = 1:numel(stim)
    s = stim(k);
    if isnan(s.cycle)
        fprintf('    %s (%.1f s): not inside any detected cycle\n', s.name, s.t);
    else
        fprintf('    %s (%.1f s): cycle %d, phase %.2f (%.0f deg)%s | PR/N %.2f | align to previous %.2f, next-to-this %.2f\n', ...
            s.name, s.t, s.cycle, s.phase_frac, rad2deg(s.phase_rad), ...
            fn_ternary(s.in_long_cycle, ' [long cycle: check]', ''), s.PR_norm, s.align_prev, s.align_next);
    end
end
end
