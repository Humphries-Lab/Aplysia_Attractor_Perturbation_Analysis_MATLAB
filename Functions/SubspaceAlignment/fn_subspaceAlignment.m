function A = fn_subspaceAlignment(V1, V2)
% FN_SUBSPACEALIGNMENT  Basis-invariant subspace alignment metric.
%   A = trace(V1' V2 V2' V1) / K, equivalently trace(P1*P2)/K where
%   Pi = Vi*Vi' are the orthogonal projectors onto each subspace.
%   V1, V2 must have orthonormal columns (K each), same ambient dim N.
%   A in [0,1]: 1 = identical subspaces, ~K/N = chance for random
%   subspaces in N dimensions.
    K = size(V1, 2);
    assert(size(V2,2) == K, 'V1 and V2 must have the same number of columns (K).');
    M = V1' * V2;
    A = trace(M * M') / K;
end
