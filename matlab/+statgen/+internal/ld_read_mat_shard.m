function [shard, meta] = ld_read_mat_shard(path, check_payload_structure)
    if nargin < 2
        check_payload_structure = false;
    end
    if exist(path, 'file') ~= 2
        error('statgen:io', 'LD file not found: %s', path);
    end
    payload = load(path);
    if ~isfield(payload, 'ld_r') || ~isfield(payload, 'a1freq') || ~isfield(payload, 'metadata')
        error('statgen:ld', '%s: missing required MAT variables', path);
    end
    meta = payload.metadata;
    validate_shard_metadata_(meta, path);

    num_snp = double(meta.num_snp);
    if ~issparse(payload.ld_r)
        error('statgen:ld', '%s: ld_r must be sparse', path);
    end
    if ~isequal(size(payload.ld_r), [num_snp, num_snp])
        error('statgen:ld', '%s: ld_r shape does not match metadata num_snp', path);
    end
    if nnz(payload.ld_r) ~= double(meta.nnz)
        error('statgen:ld', '%s: ld_r nnz does not match metadata nnz', path);
    end
    if numel(payload.a1freq) ~= num_snp
        error('statgen:ld', '%s: a1freq length does not match metadata num_snp', path);
    end
    if check_payload_structure
        statgen.internal.ld_validate_sparse_matrix_values(payload.ld_r, path);
    end

    shard = statgen.LDShard(meta.chr, meta.sex, num_snp, payload.ld_r, ...
        payload.a1freq, meta.reference_checksum);
end

function validate_shard_metadata_(meta, path)
    if ~isstruct(meta)
        error('statgen:ld', '%s: metadata must be a struct', path);
    end
    required = {'object_type', 'schema_version', 'format', 'chr', 'sex', ...
        'num_snp', 'nnz', 'matrix', 'diagonal', 'value', 'reference_checksum'};
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
    if ~strcmp(meta.matrix, 'symmetric') || ~strcmp(meta.diagonal, 'explicit_unit') || ~strcmp(meta.value, 'r')
        error('statgen:ld', '%s: invalid LD matrix metadata', path);
    end
end
