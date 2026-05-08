function save_sumstats_cache(sumstats, path, varargin)
% Save Sumstats to MATLAB binary .mat cache file.
    path = char(path);
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_sumstats_cache', 'v7', {'v7', 'v7.3', 'v5'}, varargin{:});

    n_shards = numel(sumstats.shards);
    metadata.schema = 'sumstats_cache/0.1';
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);
    metadata.has_z = ~isempty(sumstats.zvec);
    metadata.has_n = ~isempty(sumstats.nvec);
    metadata.has_beta = ~isempty(sumstats.beta_vec);
    metadata.has_se = ~isempty(sumstats.se_vec);
    metadata.has_eaf = ~isempty(sumstats.eaf_vec);
    metadata.has_info = ~isempty(sumstats.info_vec);

    for i = 1:n_shards
        s = sumstats.shards{i};
        off = sumstats.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.reference_checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
    end

    logpvec = sumstats.logpvec;
    zvec = sumstats.zvec;
    nvec = sumstats.nvec;
    beta_vec = sumstats.beta_vec;
    se_vec = sumstats.se_vec;
    eaf_vec = sumstats.eaf_vec;
    info_vec = sumstats.info_vec;

    save(path, 'metadata', 'zvec', 'nvec', 'logpvec', ...
        'beta_vec', 'se_vec', 'eaf_vec', 'info_vec', save_arg);
end
