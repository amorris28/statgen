function [labels, checksums, start0, stop0, n_shards] = validate_cache_shard_metadata(meta, schema, cache_name, num_snp, require_nonempty_shards)
%VALIDATE_CACHE_SHARD_METADATA Validate common cache shard metadata fields.
    if nargin < 5
        require_nonempty_shards = false;
    end
    if ~isfield(meta, 'schema') || ~strcmp(meta.schema, schema)
        if isfield(meta, 'schema')
            error('statgen:cache', 'Unsupported %s cache schema: %s', cache_name, meta.schema);
        end
        error('statgen:cache', 'Unsupported %s cache schema', cache_name);
    end
    required = {'n_shards', 'shard_labels', 'shard_checksums', 'shard_start0', 'shard_stop0'};
    for i = 1:numel(required)
        if ~isfield(meta, required{i})
            error('statgen:cache', 'Invalid %s cache: missing metadata field %s', cache_name, required{i});
        end
    end

    labels = statgen.internal.ensure_cell_col(meta.shard_labels);
    checksums = statgen.internal.ensure_cell_col(meta.shard_checksums);
    start0 = double(meta.shard_start0(:));
    stop0 = double(meta.shard_stop0(:));

    n_shards = numel(labels);
    if ~isscalar(meta.n_shards) || double(meta.n_shards) ~= n_shards
        error('statgen:cache', 'Invalid %s cache: n_shards mismatch', cache_name);
    end
    if n_shards <= 0
        error('statgen:cache', 'Invalid %s cache: n_shards must be positive', cache_name);
    end
    if numel(checksums) ~= n_shards || numel(start0) ~= n_shards || numel(stop0) ~= n_shards
        error('statgen:cache', 'Invalid %s cache: shard metadata length mismatch', cache_name);
    end
    if start0(1) ~= 0 || stop0(end) ~= num_snp || any(stop0 < start0) || any(start0(2:end) ~= stop0(1:end-1))
        error('statgen:cache', 'Invalid %s cache: shard offsets are not contiguous', cache_name);
    end
    if require_nonempty_shards && any(stop0 <= start0)
        error('statgen:cache', 'Invalid %s cache: shard offsets are not contiguous', cache_name);
    end
end
