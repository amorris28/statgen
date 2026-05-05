function [shard, meta] = ld_read_mat_shard(path, check_payload_structure, retain_ld_r)
    if nargin < 2
        check_payload_structure = false;
    end
    if nargin < 3 || isempty(retain_ld_r)
        retain_ld_r = true;
    end
    if exist(path, 'file') ~= 2
        error('statgen:io', 'LD file not found: %s', path);
    end
    payload = load(path);
    if ~isfield(payload, 'ld_r') || ~isfield(payload, 'a1freq') || ~isfield(payload, 'metadata')
        error('statgen:ld', '%s: missing required MAT variables', path);
    end
    meta = payload.metadata;
    statgen.internal.ld_validate_mat_metadata(meta, path);

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
        payload.a1freq, meta.reference_checksum, retain_ld_r);
end
