function panel = load_reference_cache(path, varargin)
% Load a ReferencePanel from a MATLAB binary .mat cache file.
%
%   panel = statgen.load_reference_cache(path)
%   panel = statgen.load_reference_cache(path, shards)
    shards = parse_args_(varargin{:});

    path = char(path);
    loaded = load(path, 'metadata', 'chr', 'snp', 'bp', 'a1', 'a2');

    meta = loaded.metadata;
    chr = statgen.internal.ensure_cell_col(loaded.chr);
    snp = statgen.internal.ensure_cell_col(loaded.snp);
    bp = loaded.bp(:);
    a1 = statgen.internal.ensure_cell_col(loaded.a1);
    a2 = statgen.internal.ensure_cell_col(loaded.a2);

    n = numel(bp);
    if numel(chr) ~= n || numel(snp) ~= n || numel(a1) ~= n || numel(a2) ~= n
        error('statgen:cache', 'Invalid reference cache: panel-wide vector lengths mismatch');
    end

    [labels, checksums, start0, stop0] = validate_metadata_(meta, n);
    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_reference_cache');
    shard_objs = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);

        shard_obj = statgen.ReferenceShard( ...
            labels{idx}, chr(ix), snp(ix), bp(ix), a1(ix), a2(ix), ...
            checksums{idx});
        shard_objs{i} = shard_obj;
    end

    panel = statgen.ReferencePanel(shard_objs);
end

function shards = parse_args_(varargin)
    shards = [];
    if numel(varargin) > 1
        error('statgen:arg', 'load_reference_cache accepts at most path and shards');
    end
    if numel(varargin) == 1
        shards = varargin{1};
    end
end

function [labels, checksums, start0, stop0] = validate_metadata_(meta, num_snp)
    if ~strcmp(meta.schema, 'reference_cache/0.1')
        error('statgen:cache', 'Unsupported reference cache schema: %s', meta.schema);
    end
    labels = statgen.internal.ensure_cell_col(meta.shard_labels);
    checksums = statgen.internal.ensure_cell_col(meta.shard_checksums);
    start0 = double(meta.shard_start0(:));
    stop0 = double(meta.shard_stop0(:));

    n_shards = numel(labels);
    if ~isfield(meta, 'n_shards') || ~isscalar(meta.n_shards) || meta.n_shards ~= n_shards
        error('statgen:cache', 'Invalid reference cache: n_shards mismatch');
    end
    if numel(checksums) ~= n_shards || numel(start0) ~= n_shards || numel(stop0) ~= n_shards
        error('statgen:cache', 'Invalid reference cache: shard metadata length mismatch');
    end
    validate_offsets_(start0, stop0, num_snp, 'reference');
end

function validate_offsets_(start0, stop0, num_snp, label)
    if isempty(start0)
        if num_snp ~= 0
            error('statgen:cache', 'Invalid %s cache: empty shard offsets for non-empty payload', label);
        end
        return
    end
    if start0(1) ~= 0 || stop0(end) ~= num_snp || any(stop0 < start0) || any(start0(2:end) ~= stop0(1:end-1))
        error('statgen:cache', 'Invalid %s cache: shard offsets are not contiguous', label);
    end
end
