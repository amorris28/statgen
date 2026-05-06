function save_reference_cache(panel, path, varargin)
% Save a ReferencePanel to a MATLAB binary .mat cache file.
    path = char(path);
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_reference_cache', 'v7', {'v7', 'v7.3', 'v5'}, varargin{:});

    n_shards = numel(panel.shards);
    metadata.schema = 'reference_cache/0.1';
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);

    for i = 1:n_shards
        s = panel.shards{i};
        off = panel.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
    end

    chr = panel.chr;
    snp = panel.snp;
    bp = panel.bp;
    a1 = panel.a1;
    a2 = panel.a2;

    save(path, 'metadata', 'chr', 'snp', 'bp', 'a1', 'a2', save_arg);
end
