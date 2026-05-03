function ld_validate_sparse_matrix_values(mat, path)
    if size(mat, 1) ~= size(mat, 2)
        error('statgen:ld', '%s: ld_r must be square', path);
    end
    d = diag(mat);
    if any(abs(d - 1) > 1e-7)
        error('statgen:ld', '%s: ld_r diagonal must be explicit unit', path);
    end
    diff = mat - mat';
    if nnz(diff) > 0 && max(abs(nonzeros(diff))) > 1e-6
        error('statgen:ld', '%s: ld_r must be symmetric', path);
    end
end
