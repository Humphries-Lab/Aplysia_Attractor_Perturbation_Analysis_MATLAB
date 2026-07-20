function first_idx = fn_firstConsecRun(logical_vec, n_consec)
%FN_FIRST_CONSEC_RUN  Returns the index of the FIRST element in the first
%  run of n_consec consecutive TRUE values in logical_vec.
%  Returns [] if no such run exists.
%
%  Used to detect attractor onset/return: the first time RR stays above
%  threshold for n_consec consecutive windows.

first_idx = [];
logical_vec = logical_vec(:)';   % ensure row vector
n = numel(logical_vec);
count = 0;
for i = 1:n
    if logical_vec(i)
        count = count + 1;
        if count >= n_consec
            first_idx = i - n_consec + 1;   % index of the RUN START
            return;
        end
    else
        count = 0;
    end
end
end
