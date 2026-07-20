%% VALIDATE_SUBSPACEALIGNMENT
%  Synthetic validation suite for fn_subspaceAlignment: checks the metric
%  against known theoretical values across five controlled scenarios, with
%  repeated trials, error bars, and automated PASS/FAIL checks.
%
%  This is a unit test of the toolbox, not part of the analysis pipeline -
%  it uses no real recording data. Run it whenever fn_subspaceAlignment (or
%  its dependents) changes, to confirm the metric still behaves as
%  expected.
%
%  Metric definition (basis-invariant):
%     A(V1,V2) = trace(V1' V2 V2' V1) / K = trace(P1*P2) / K,  Pi = Vi*Vi'
%  This equals the mean squared cosine of the principal angles between the
%  two K-dimensional subspaces, and depends only on the SPAN of V1 and V2,
%  not the particular orthonormal basis chosen for each.

clear; close all; clc;

toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(toolboxRoot, 'Functions')));
run(fullfile(toolboxRoot, 'Config_UserSettings.m'));   % defines `cfg`
if ~exist(cfg.FIGURES_DIR,'dir'), mkdir(cfg.FIGURES_DIR); end

rng(cfg.validation_seed_alignment);
nReps = cfg.validation_nReps;
tol   = cfg.validation_tol;

syn_N = 100;   % ambient dimension (e.g., neurons)
syn_K = 10;    % subspace dimension (e.g., top PCs)

fprintf('\n=== Synthetic Validation: Subspace Alignment ===\n');
results = struct();

% -------------------------------------------------------------------
%  TEST 1: Identical subspaces -> A should be exactly 1
% ---------------------------------------------------------------------
syn_V1   = orth(randn(syn_N, syn_K));
syn_A_id = fn_subspaceAlignment(syn_V1, syn_V1);
fprintf('  Test 1 (Identical)   : A = %.4f  (Expected: 1.000)  [%s]\n', ...
    syn_A_id, fn_passFail(syn_A_id, 1.0, tol));
results.test1_value = syn_A_id;

% -------------------------------------------------------------------
%  TEST 2: Random/orthogonal subspaces, swept across K/N ratios.
%  Expected chance level ~ K/N for random subspaces in N dims.
% ---------------------------------------------------------------------
K_range  = [2 5 10 20 40 60 80];
N_fixed  = 100;
A2_mean  = zeros(size(K_range));
A2_std   = zeros(size(K_range));
for i = 1:length(K_range)
    K = K_range(i);
    vals = zeros(nReps,1);
    for r = 1:nReps
        Va = orth(randn(N_fixed, K));
        Vb = orth(randn(N_fixed, K));
        vals(r) = fn_subspaceAlignment(Va, Vb);
    end
    A2_mean(i) = mean(vals); A2_std(i) = std(vals);
end
expected2 = K_range / N_fixed;
fprintf('\n  Test 2 (Random, K/N sweep):\n');
for i = 1:length(K_range)
    fprintf('    K=%3d, N=%3d : A = %.4f +/- %.4f  (Expected K/N = %.4f)  [%s]\n', ...
        K_range(i), N_fixed, A2_mean(i), A2_std(i), expected2(i), fn_passFail(A2_mean(i), expected2(i), tol));
end
results.test2_K = K_range; results.test2_N = N_fixed;
results.test2_mean = A2_mean; results.test2_std = A2_std; results.test2_expected = expected2;

% -------------------------------------------------------------------
%  TEST 3: Exact partial-overlap subspaces, swept 0% -> 100%. Built by
%  construction (shared basis + private orthogonal complements) so the
%  true overlap fraction is known exactly. The private complements are
%  independent random K_unique-frames drawn from an (N-K_shared)-dim
%  space, so they carry their own chance-level alignment:
%     A_theory(f) = [K_shared + K_unique^2 / (N - K_shared)] / K
% ---------------------------------------------------------------------
overlap_fracs = 0:0.1:1;
K3 = syn_K;
A3_mean = zeros(size(overlap_fracs)); A3_std = zeros(size(overlap_fracs));
expected3 = zeros(size(overlap_fracs));
for i = 1:length(overlap_fracs)
    f = overlap_fracs(i);
    K_shared = round(f * K3); K_unique = K3 - K_shared;
    vals = zeros(nReps,1);
    for r = 1:nReps
        [Va, Vb] = fn_makePartialOverlapSubspaces(syn_N, K_shared, K_unique);
        vals(r) = fn_subspaceAlignment(Va, Vb);
    end
    A3_mean(i) = mean(vals); A3_std(i) = std(vals);
    if (syn_N - K_shared) > 0
        expected3(i) = (K_shared + (K_unique^2) / (syn_N - K_shared)) / K3;
    else
        expected3(i) = 1;
    end
end
fprintf('\n  Test 3 (Exact overlap fraction sweep, bias-corrected expectation):\n');
for i = 1:length(overlap_fracs)
    fprintf('    Overlap=%.1f : A = %.4f +/- %.4f  (Expected = %.4f)  [%s]\n', ...
        overlap_fracs(i), A3_mean(i), A3_std(i), expected3(i), fn_passFail(A3_mean(i), expected3(i), tol));
end
results.test3_frac = overlap_fracs;
results.test3_mean = A3_mean; results.test3_std = A3_std; results.test3_expected = expected3;
[~, idx50] = min(abs(overlap_fracs - 0.5));
fprintf('  --> Test 3 (50%% Overlap), single check: A = %.4f  (Expected: %.4f)  [%s]\n', ...
    A3_mean(idx50), expected3(idx50), fn_passFail(A3_mean(idx50), expected3(idx50), tol));

% -------------------------------------------------------------------
%  TEST 4: Dimensionality control - fix K, sweep ambient dimension N.
%  Isolates whether chance level scales as K/N, not as a function of K or
%  N individually.
% ---------------------------------------------------------------------
Ns = [20 50 100 200 500 1000]; K_fixed = 10;
A4_mean = zeros(size(Ns)); A4_std = zeros(size(Ns));
for i = 1:length(Ns)
    N = Ns(i);
    vals = zeros(nReps,1);
    for r = 1:nReps
        Va = orth(randn(N, K_fixed)); Vb = orth(randn(N, K_fixed));
        vals(r) = fn_subspaceAlignment(Va, Vb);
    end
    A4_mean(i) = mean(vals); A4_std(i) = std(vals);
end
expected4 = K_fixed ./ Ns;
fprintf('\n  Test 4 (Fixed K=%d, sweep N):\n', K_fixed);
for i = 1:length(Ns)
    fprintf('    N=%4d : A = %.4f +/- %.4f  (Expected K/N = %.4f)  [%s]\n', ...
        Ns(i), A4_mean(i), A4_std(i), expected4(i), fn_passFail(A4_mean(i), expected4(i), tol));
end
results.test4_N = Ns; results.test4_K = K_fixed;
results.test4_mean = A4_mean; results.test4_std = A4_std; results.test4_expected = expected4;

% -------------------------------------------------------------------
%  TEST 5: Noise-interpolation control. Start from an identical subspace
%  (A=1) and progressively mix in an independent random subspace;
%  alignment should decay monotonically toward the chance level K/N.
% ---------------------------------------------------------------------
noise_levels = linspace(0, 1, 11);
A5_mean = zeros(size(noise_levels)); A5_std = zeros(size(noise_levels));
for i = 1:length(noise_levels)
    eps_lvl = noise_levels(i);
    vals = zeros(nReps,1);
    for r = 1:nReps
        V1 = orth(randn(syn_N, syn_K));
        Vrand = orth(randn(syn_N, syn_K));
        Vnoisy = orth((1 - eps_lvl) * V1 + eps_lvl * Vrand);
        vals(r) = fn_subspaceAlignment(V1, Vnoisy);
    end
    A5_mean(i) = mean(vals); A5_std(i) = std(vals);
end
chance_level = syn_K / syn_N;
fprintf('\n  Test 5 (Noise interpolation 0 -> 1):\n');
for i = 1:length(noise_levels)
    fprintf('    noise=%.1f : A = %.4f +/- %.4f\n', noise_levels(i), A5_mean(i), A5_std(i));
end
is_monotonic = all(diff(A5_mean) <= 1e-6);
fprintf('  --> Monotonic decay check: [%s]\n', fn_ternary(is_monotonic, 'PASS', 'FAIL'));
results.test5_noise = noise_levels; results.test5_mean = A5_mean; results.test5_std = A5_std;
results.chance_level = chance_level;

%% =========================================================================
%  SUMMARY FIGURE
%  =========================================================================
figure('Name', 'Subspace Alignment - Synthetic Validation', 'Position', [80 80 1200 800], 'Color', 'w');

subplot(2,3,1);
bar(1, syn_A_id, 0.4, 'FaceColor', [0.2 0.5 0.8]); hold on;
h_exp = yline(1, 'k--', 'LineWidth', 1.5);
ylim([0 1.1]); xlim([0.5 1.5]); set(gca, 'XTick', []);
ylabel('Alignment A'); title('Test 1: Identical Subspaces');
legend(h_exp, 'Expected = 1.0', 'Location', 'southeast', 'Box', 'off'); box off;

subplot(2,3,2);
kn_ratio = K_range / N_fixed;
errorbar(kn_ratio, A2_mean, A2_std, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.2 0.5 0.8], 'Color', [0.2 0.5 0.8]); hold on;
plot([0 max(kn_ratio)], [0 max(kn_ratio)], 'k--', 'LineWidth', 1.5);
xlabel('K / N ratio'); ylabel('Alignment A'); title('Test 2: Random Subspaces vs K/N');
legend({'Observed', 'Expected (y = K/N)'}, 'Location', 'northwest', 'Box', 'off'); box off;

subplot(2,3,3);
errorbar(overlap_fracs, A3_mean, A3_std, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.8 0.3 0.3], 'Color', [0.8 0.3 0.3]); hold on;
plot(overlap_fracs, expected3, 'k--', 'LineWidth', 1.5);
xlabel('True overlap fraction'); ylabel('Alignment A'); title('Test 3: Exact Partial Overlap');
legend({'Observed', 'Expected (bias-corrected)'}, 'Location', 'northwest', 'Box', 'off'); box off;

subplot(2,3,4);
kn_ratio4 = K_fixed ./ Ns;
errorbar(kn_ratio4, A4_mean, A4_std, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.3 0.7 0.4], 'Color', [0.3 0.7 0.4]); hold on;
plot([0 max(kn_ratio4)], [0 max(kn_ratio4)], 'k--', 'LineWidth', 1.5);
set(gca, 'XScale', 'log');
xlabel('K / N ratio (log scale, K fixed = 10)'); ylabel('Alignment A'); title('Test 4: Dimensionality Control');
legend({'Observed', 'Expected (y = K/N)'}, 'Location', 'northwest', 'Box', 'off'); box off;

subplot(2,3,5);
errorbar(noise_levels, A5_mean, A5_std, 'o-', 'LineWidth', 1.5, 'MarkerFaceColor', [0.6 0.4 0.8], 'Color', [0.6 0.4 0.8]); hold on;
h_ident = yline(1, 'k--', 'LineWidth', 1.5);
h_chance = yline(chance_level, 'r--', 'LineWidth', 1.5);
xlabel('Noise mixing fraction'); ylabel('Alignment A'); title('Test 5: Noise-Interpolation Control'); ylim([0 1.05]);
legend([h_ident, h_chance], {'A = 1 (identical)', sprintf('Chance = K/N = %.2f', chance_level)}, ...
       'Location', 'northeast', 'Box', 'off'); box off;

subplot(2,3,6);
K_hist = 10; N_hist = 100; nHistReps = 500;
null_vals = zeros(nHistReps,1);
for r = 1:nHistReps
    Va = orth(randn(N_hist, K_hist)); Vb = orth(randn(N_hist, K_hist));
    null_vals(r) = fn_subspaceAlignment(Va, Vb);
end
histogram(null_vals, 30, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'Normalization', 'pdf'); hold on;
h_emp = xline(mean(null_vals), 'b-', 'LineWidth', 2);
h_theo = xline(K_hist/N_hist, 'k--', 'LineWidth', 2);
xlabel('Alignment A'); ylabel('Probability Density'); title(sprintf('Null Distribution (K=%d, N=%d)', K_hist, N_hist));
legend([h_emp, h_theo], {'Empirical mean', 'Theoretical K/N'}, 'Location', 'northwest', 'Box', 'off'); box off;

sgtitle('Subspace Alignment: Synthetic Validation', 'FontSize', 16, 'FontWeight', 'bold', 'FontName', 'Arial');
all_axes = findobj(gcf, 'Type', 'axes');
set(all_axes, 'FontName', 'Arial', 'FontSize', 11, 'LineWidth', 1.2, 'TickDir', 'out', ...
    'TickLength', [0.015 0.015], 'XColor', 'k', 'YColor', 'k');
for i = 1:length(all_axes)
    grid(all_axes(i), 'off'); all_axes(i).Box = 'off';
end
exportgraphics(gcf, fullfile(cfg.FIGURES_DIR, 'validation_subspace_alignment.png'), 'Resolution', 500);
close(gcf);

fprintf('\n=== Synthetic Validation complete. ===\n');
