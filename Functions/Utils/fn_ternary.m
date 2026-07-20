function out = fn_ternary(cond, a, b)
% TERNARY  Small helper for inline conditional strings.
    if cond
        out = a;
    else
        out = b;
    end
end
