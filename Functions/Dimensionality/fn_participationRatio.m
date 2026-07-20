function PR = fn_participationRatio(eigenvalues)
eigenvalues=eigenvalues(:); eigenvalues(eigenvalues<0)=0;
if sum(eigenvalues)==0, PR=0; return; end
PR=sum(eigenvalues)^2/sum(eigenvalues.^2);
end
