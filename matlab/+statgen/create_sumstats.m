function sumstats = create_sumstats(reference, pvec, zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)
% Create Sumstats from full-panel aligned vectors.
    if nargin < 2
        error('statgen:arg', 'create_sumstats requires reference and pvec');
    end
    if nargin < 3, zvec = []; end
    if nargin < 4, nvec = []; end
    if nargin < 5, beta_vec = []; end
    if nargin < 6, se_vec = []; end
    if nargin < 7, eaf_vec = []; end
    if nargin < 8, info_vec = []; end

    n = double(reference.num_snp);
    p_aligned = coerce_pvec_(pvec, n, 'pvec');
    logp_aligned = statgen.internal.sumstats_derive_logp(p_aligned);

    z_aligned = coerce_optional_vec_(zvec, n, 'zvec');
    n_aligned = coerce_optional_vec_(nvec, n, 'nvec');
    beta_aligned = coerce_optional_vec_(beta_vec, n, 'beta_vec');
    se_aligned = coerce_optional_vec_(se_vec, n, 'se_vec');
    eaf_aligned = coerce_optional_vec_(eaf_vec, n, 'eaf_vec');
    info_aligned = coerce_optional_vec_(info_vec, n, 'info_vec');

    n_shards = numel(reference.shards);
    out_shards = cell(n_shards, 1);
    for i = 1:n_shards
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        out_shards{i} = statgen.SumstatsShard( ...
            s_ref.label, s_ref.checksum, ...
            logp_aligned(ix), ...
            pick_optional_(z_aligned, ix), ...
            pick_optional_(n_aligned, ix), ...
            pick_optional_(beta_aligned, ix), ...
            pick_optional_(se_aligned, ix), ...
            pick_optional_(eaf_aligned, ix), ...
            pick_optional_(info_aligned, ix));
    end
    sumstats = statgen.Sumstats(out_shards);
end

function out = coerce_pvec_(x, n, name)
    out = coerce_optional_vec_(x, n, name);
    if isempty(out)
        error('statgen:sumstats', '%s is required', name);
    end
    bad = ~isfinite(out) | out < 0 | out > 1;
    if any(bad)
        idx = find(bad, 1, 'first');
        error('statgen:sumstats', '%s(%d) must be finite numeric in [0, 1]', name, idx);
    end
end

function out = coerce_optional_vec_(x, n, name)
    if isempty(x)
        out = [];
        return
    end
    out = coerce_vec_(x, n, name);
    bad = ~(isfinite(out) | isnan(out));
    if any(bad)
        idx = find(bad, 1, 'first');
        error('statgen:sumstats', '%s(%d) must be finite numeric or NaN', name, idx);
    end
end

function out = coerce_vec_(x, n, name)
    try
        out = double(x);
    catch
        error('statgen:sumstats', '%s must be numeric', name);
    end
    if ~(isvector(out) || isempty(out))
        error('statgen:sumstats', '%s must be a vector with length %d', name, n);
    end
    out = out(:);
    if numel(out) ~= n
        error('statgen:sumstats', '%s length mismatch: expected %d, got %d', name, n, numel(out));
    end
end

function out = pick_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end
