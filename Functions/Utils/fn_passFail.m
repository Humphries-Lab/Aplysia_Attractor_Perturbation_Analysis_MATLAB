function s = fn_passFail(observed, expected, tol)
% PASSFAIL  Simple tolerance check, returns 'PASS' or 'FAIL' string.
    if abs(observed - expected) <= tol
        s = 'PASS';
    else
        s = 'FAIL';
    end
end
