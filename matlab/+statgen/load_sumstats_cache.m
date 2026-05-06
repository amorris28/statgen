function sumstats = load_sumstats_cache(path, shards)
% Load Sumstats from MATLAB .mat cache file with optional shard subsetting.
    if nargin < 2
        shards = [];
    end
    path = char(path);
    loaded = load(path, 'metadata', 'zvec', 'nvec', 'logpvec', ...
        'beta_vec', 'se_vec', 'eaf_vec', 'info_vec');
    meta = loaded.metadata;

    zvec = loaded.zvec(:);
    nvec = loaded.nvec(:);
    logpvec = loaded.logpvec(:);
    num_snp = numel(zvec);
    if numel(nvec) ~= num_snp || numel(logpvec) ~= num_snp
        error('statgen:cache', 'Invalid sumstats cache: required vector lengths mismatch');
    end

    [labels, checksums, start0, stop0] = validate_metadata_(meta, num_snp);
    beta_vec = validate_optional_(loaded.beta_vec, meta.has_beta, num_snp, 'beta_vec');
    se_vec = validate_optional_(loaded.se_vec, meta.has_se, num_snp, 'se_vec');
    eaf_vec = validate_optional_(loaded.eaf_vec, meta.has_eaf, num_snp, 'eaf_vec');
    info_vec = validate_optional_(loaded.info_vec, meta.has_info, num_snp, 'info_vec');

    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_sumstats_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);
        out_shards{i} = statgen.SumstatsShard( ...
            labels{idx}, checksums{idx}, ...
            zvec(ix), nvec(ix), logpvec(ix), ...
            slice_optional_(beta_vec, ix), ...
            slice_optional_(se_vec, ix), ...
            slice_optional_(eaf_vec, ix), ...
            slice_optional_(info_vec, ix));
    end
    sumstats = statgen.Sumstats(out_shards);
end

function [labels, checksums, start0, stop0] = validate_metadata_(meta, num_snp)
    if ~strcmp(meta.schema, 'sumstats_cache/0.1')
        error('statgen:cache', 'Unsupported sumstats cache schema: %s', meta.schema);
    end
    labels = statgen.internal.ensure_cell_col(meta.shard_labels);
    checksums = statgen.internal.ensure_cell_col(meta.shard_checksums);
    start0 = double(meta.shard_start0(:));
    stop0 = double(meta.shard_stop0(:));

    n_shards = numel(labels);
    if ~isfield(meta, 'n_shards') || ~isscalar(meta.n_shards) || meta.n_shards ~= n_shards
        error('statgen:cache', 'Invalid sumstats cache: n_shards mismatch');
    end
    if numel(checksums) ~= n_shards || numel(start0) ~= n_shards || numel(stop0) ~= n_shards
        error('statgen:cache', 'Invalid sumstats cache: shard metadata length mismatch');
    end
    validate_offsets_(start0, stop0, num_snp);
end

function validate_offsets_(start0, stop0, num_snp)
    if isempty(start0)
        if num_snp ~= 0
            error('statgen:cache', 'Invalid sumstats cache: empty shard offsets for non-empty payload');
        end
        return
    end
    if start0(1) ~= 0 || stop0(end) ~= num_snp || any(stop0 < start0) || any(start0(2:end) ~= stop0(1:end-1))
        error('statgen:cache', 'Invalid sumstats cache: shard offsets are not contiguous');
    end
end

function out = validate_optional_(vec, has_field, num_snp, name)
    if has_field
        out = vec(:);
        if numel(out) ~= num_snp
            error('statgen:cache', 'Invalid sumstats cache: %s length mismatch', name);
        end
    else
        if ~isempty(vec)
            error('statgen:cache', 'Invalid sumstats cache: absent %s must be empty', name);
        end
        out = [];
    end
end

function out = slice_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end
