function ok = ld_same_optional_string(a, b)
    if isempty(a) && isempty(b)
        ok = true;
    elseif isempty(a) || isempty(b)
        ok = false;
    else
        ok = strcmp(char(a), char(b));
    end
end
