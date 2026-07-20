function [Vsp1, Vsp2] = fn_makePartialOverlapSubspaces(N, Kshared, Kunique)
% MAKE_PARTIAL_OVERLAP_SUBSPACES  Construct two K-dim subspaces (K =
%   Kshared + Kunique) that share an EXACT Kshared-dimensional subspace
%   and have independent, mutually orthogonal private complements.
%   NOTE: the two K-dim subspaces are exact in their shared component,
%   but the private complements U1, U2 are independent random draws and
%   so are not orthogonal to each other in general — see Test 3 fix note
%   above for the correct theoretical alignment this implies.
    if Kshared > 0
        Vshared = orth(randn(N, Kshared));
    else
        Vshared = zeros(N, 0);
    end

    Pshared = Vshared * Vshared';
    Icomp = eye(N) - Pshared;   % projector onto orthogonal complement of shared space

    if Kunique > 0
        U1 = orth(Icomp * randn(N, Kunique));
        U2 = orth(Icomp * randn(N, Kunique));
    else
        U1 = zeros(N, 0);
        U2 = zeros(N, 0);
    end

    Vsp1 = [Vshared, U1];
    Vsp2 = [Vshared, U2];
end
