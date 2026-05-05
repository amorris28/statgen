function ld_validate_mat_metadata(meta, path)
    if ~isstruct(meta)
        error('statgen:ld', '%s: metadata must be a struct', path);
    end
    required = {'object_type', 'schema_version', 'format', 'chr', 'sex', ...
        'num_snp', 'nnz', 'matrix', 'diagonal', 'value', ...
        'reference_checksum', 'reference_bim', 'num_monomorphic_snps'};
    for i = 1:numel(required)
        if ~isfield(meta, required{i})
            error('statgen:ld', '%s: metadata missing required field %s', path, required{i});
        end
    end
    if ~strcmp(meta.object_type, 'ld_shard')
        error('statgen:ld', '%s: metadata object_type must be ld_shard', path);
    end
    if ~strcmp(meta.schema_version, '1.0')
        error('statgen:ld', '%s: unsupported metadata schema_version', path);
    end
    if ~strcmp(meta.format, 'statgen_ld_mat_sparse_double')
        error('statgen:ld', '%s: metadata format must be statgen_ld_mat_sparse_double', path);
    end
    statgen.internal.ld_validate_chr_sex(meta.chr, meta.sex, [path ': metadata']);
    if double(meta.num_snp) <= 0 || floor(double(meta.num_snp)) ~= double(meta.num_snp)
        error('statgen:ld', '%s: metadata num_snp must be a positive integer', path);
    end
    if double(meta.nnz) < 0 || floor(double(meta.nnz)) ~= double(meta.nnz)
        error('statgen:ld', '%s: metadata nnz must be a non-negative integer', path);
    end
    if double(meta.num_monomorphic_snps) < 0 || ...
            floor(double(meta.num_monomorphic_snps)) ~= double(meta.num_monomorphic_snps) || ...
            double(meta.num_monomorphic_snps) > double(meta.num_snp)
        error('statgen:ld', '%s: metadata num_monomorphic_snps must be an integer in [0, num_snp]', path);
    end
    if ~strcmp(meta.matrix, 'symmetric') || ~strcmp(meta.diagonal, 'explicit_unit') || ~strcmp(meta.value, 'r')
        error('statgen:ld', '%s: invalid LD matrix metadata', path);
    end
    validate_reference_bim_(meta.reference_bim, path);
end

function validate_reference_bim_(value, path)
    if isempty(value) || ~(ischar(value) || isstring(value))
        error('statgen:ld', '%s: metadata reference_bim must be a non-empty filename', path);
    end
    value = char(value);
    if any(value == '/') || any(value == '\') || strcmp(value, '.') || ~isempty(strfind(value, '..'))
        error('statgen:ld', '%s: metadata reference_bim must be a plain relative filename', path);
    end
    if ~isempty(regexp(value, '^[A-Za-z]:', 'once'))
        error('statgen:ld', '%s: metadata reference_bim must be a plain relative filename', path);
    end
end
