function A = fn_alignTwoWindows(WPCA_w1, WPCA_w2, nDims)
%FN_ALIGNTWOWINDOWS  Takes eigendecomposition from two WPCA window structs.
%  A = (1/K) * ||Q1'Q2||_F^2   (Cunningham & Yu 2014)
K = min(nDims, min(size(WPCA_w1.eigenvectors,2), size(WPCA_w2.eigenvectors,2)));
if K == 0, A = 0; return; end
A = fn_subspaceAlignment(WPCA_w1.eigenvectors(:,1:K), WPCA_w2.eigenvectors(:,1:K));
end
