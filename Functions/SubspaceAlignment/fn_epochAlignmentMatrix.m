function align_mat = fn_epochAlignmentMatrix(WPCA, nDims)
% FN_EPOCHALIGNMENTMATRIX  Pairwise subspace alignment between every pair
% of epochs in a windowed-PCA struct array.
%
%   align_mat = FN_EPOCHALIGNMENTMATRIX(WPCA, nDims)
%
% INPUTS
%   WPCA  - 1 x nW struct array from fn_windowedPCA
%   nDims - subspace dimension K used for alignment
%
% OUTPUT
%   align_mat - nW x nW matrix, align_mat(i,j) = fn_alignTwoWindows(WPCA(i),WPCA(j),nDims)

nW = numel(WPCA);
align_mat = zeros(nW, nW);
for i = 1:nW
    for j = 1:nW
        align_mat(i,j) = fn_alignTwoWindows(WPCA(i), WPCA(j), nDims);
    end
end

end
