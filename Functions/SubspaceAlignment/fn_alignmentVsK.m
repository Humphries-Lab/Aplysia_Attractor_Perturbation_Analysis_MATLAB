function al_vs_k = fn_alignmentVsK(WPCA_i, WPCA_j, K_range)
% FN_ALIGNMENTVSK  Subspace alignment between two epochs, swept across a
% range of subspace dimensions K, to check sensitivity of the alignment
% value to the choice of K.
%
%   al_vs_k = FN_ALIGNMENTVSK(WPCA_i, WPCA_j, K_range)
%
% INPUTS
%   WPCA_i, WPCA_j - windowed-PCA structs for the two epochs being compared
%   K_range        - vector of K values to test
%
% OUTPUT
%   al_vs_k - same size as K_range, alignment value at each K

al_vs_k = zeros(size(K_range));
for k_idx = 1:numel(K_range)
    k = K_range(k_idx);
    al_vs_k(k_idx) = fn_subspaceAlignment( ...
        WPCA_i.eigenvectors(:,1:k), WPCA_j.eigenvectors(:,1:k));
end

end
