function annotation_validate_declared_binary(A, is_binary, name)
% Validate declared-binary annotation columns using stored sparse values only.
    if ~any(is_binary)
        return
    end
    vals = nonzeros(A(:, is_binary));
    bad = (vals ~= 0) & (vals ~= 1);
    if any(bad)
        v = vals(find(bad, 1, 'first'));
        error('statgen:annotations', '%s declared binary must be binary (0/1); found non-binary value %g', name, v);
    end
end
