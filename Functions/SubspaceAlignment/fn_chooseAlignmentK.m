function [nDims, chance_lvl] = fn_chooseAlignmentK(K_ALIGN, WPCA_evoked, nN, varThreshPct)
% FN_CHOOSEALIGNMENTK  Choose the subspace dimension K used for subspace
% alignment: either 'auto' (dims explaining varThreshPct of the evoked
% epoch's variance) or a fixed integer.
%
%   [nDims, chance_lvl] = FN_CHOOSEALIGNMENTK(K_ALIGN, WPCA_evoked, nN, varThreshPct)
%
% INPUTS
%   K_ALIGN      - 'auto', or a positive integer
%   WPCA_evoked  - windowed-PCA struct for the Evoked epoch (from fn_windowedPCA)
%   nN           - number of neurons (ambient dimension)
%   varThreshPct - (optional) variance threshold for 'auto' mode, in percent
%                  (default 80)
%
% OUTPUTS
%   nDims      - chosen subspace dimension K
%   chance_lvl - chance-level alignment K/N for random subspaces

if nargin < 4 || isempty(varThreshPct)
    varThreshPct = 80;
end

if ischar(K_ALIGN) && strcmpi(K_ALIGN, 'auto')
    cum_var_ev = cumsum(WPCA_evoked.var_explained);
    nDims = find(cum_var_ev >= varThreshPct, 1, 'first');
    if isempty(nDims) || nDims > nN
        nDims = min(10, nN);
    end
    fprintf('  K auto = %d  (explains %.1f%% evoked variance)\n', nDims, cum_var_ev(nDims));
else
    nDims = min(K_ALIGN, nN);
    fprintf('  K fixed = %d\n', nDims);
end

chance_lvl = nDims / nN;
fprintf('  Chance K/N = %d/%d = %.3f\n', nDims, nN, chance_lvl);

end
