function out = annotation_coerce_metadata_vector(value, expected_len, name)
% Coerce annotation metadata to an expected-length cell-string column.
    if ischar(value)
        error('statgen:annotations', '%s must be a string vector with length %d', name, expected_len);
    elseif isstring(value)
        value = cellstr(value(:));
    end
    if ~iscell(value)
        error('statgen:annotations', '%s must be a string vector with length %d', name, expected_len);
    end
    out = statgen.internal.ensure_cell_col(value);
    if numel(out) ~= expected_len
        error('statgen:annotations', '%s length mismatch: expected %d, got %d', name, expected_len, numel(out));
    end
end
