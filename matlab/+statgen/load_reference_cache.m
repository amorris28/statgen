function panel = load_reference_cache(path, varargin)
%LOAD_REFERENCE_CACHE Load a cached ReferencePanel.
%
%   panel = statgen.load_reference_cache(path)
%   panel = statgen.load_reference_cache(path, shards)
%
% Loads one MATLAB .mat reference cache file produced by
% statgen.save_reference_cache. Optional shards return a logical subset of the
% cached shard labels. The loaded ReferencePanel state is read from cache
% metadata: full references expose snp, a1, and a2; thin references support
% alignment checks and summary-statistics matching but do not expose those full
% reference fields.
%
% See also statgen.ReferencePanel, statgen.load_reference,
% statgen.save_reference_cache.
    shards = parse_args_(varargin{:});

    path = char(path);
    meta_loaded = load(path, 'metadata');
    require_field_(meta_loaded, 'metadata');
    meta = meta_loaded.metadata;
    mode = validate_cache_mode_(meta);
    if strcmp(mode, 'full')
        loaded = load(path, 'chr', 'snp', 'bp', 'a1', 'a2', 'a1_hash64', 'a2_hash64');
    else
        loaded = load(path, 'bp', 'a1_hash64', 'a2_hash64');
    end

    require_field_(loaded, 'bp');
    require_field_(loaded, 'a1_hash64');
    require_field_(loaded, 'a2_hash64');
    bp = loaded.bp(:);
    a1_hash64 = uint64(loaded.a1_hash64(:));
    a2_hash64 = uint64(loaded.a2_hash64(:));

    n = numel(bp);
    if numel(a1_hash64) ~= n || numel(a2_hash64) ~= n
        error('statgen:cache', 'Invalid reference cache: allele hash vector lengths mismatch');
    end

    if strcmp(mode, 'full')
        require_field_(loaded, 'chr');
        require_field_(loaded, 'snp');
        require_field_(loaded, 'a1');
        require_field_(loaded, 'a2');
        chr = statgen.internal.ensure_cell_col(loaded.chr);
        snp = statgen.internal.ensure_cell_col(loaded.snp);
        a1 = statgen.internal.ensure_cell_col(loaded.a1);
        a2 = statgen.internal.ensure_cell_col(loaded.a2);
        if numel(chr) ~= n || numel(snp) ~= n || numel(a1) ~= n || numel(a2) ~= n
            error('statgen:cache', 'Invalid reference cache: panel-wide vector lengths mismatch');
        end
    end

    [labels, checksums, start0, stop0] = validate_metadata_(meta, n);
    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_reference_cache');
    shard_objs = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);

        if strcmp(mode, 'full')
            shard_obj = statgen.ReferenceShard( ...
                labels{idx}, chr(ix), snp(ix), bp(ix), a1(ix), a2(ix), ...
                checksums{idx});
        else
            shard_obj = statgen.ReferenceShard.from_thin( ...
                labels{idx}, numel(ix), bp(ix), a1_hash64(ix), a2_hash64(ix), ...
                checksums{idx});
        end
        shard_objs{i} = shard_obj;
    end

    panel = statgen.ReferencePanel(shard_objs);
end

function shards = parse_args_(varargin)
    shards = [];
    if numel(varargin) > 1
        error('statgen:arg', 'load_reference_cache accepts path and optional shards; cache mode is saved in metadata');
    end
    if ~isempty(varargin)
        if is_name_(varargin{1}, 'full') || is_name_(varargin{1}, 'mode')
            error('statgen:arg', 'load_reference_cache no longer accepts full/mode options; save full or thin cache mode with save_reference_cache');
        end
        shards = varargin{1};
    end
end

function tf = is_name_(x, name)
    tf = (ischar(x) || isstring(x)) && strcmpi(char(x), name);
end

function require_field_(loaded, name)
    if ~isfield(loaded, name)
        if strcmp(name, 'a1_hash64') || strcmp(name, 'a2_hash64')
            error('statgen:cache', 'reference cache missing a1_hash64/a2_hash64; rebuild cache');
        end
        error('statgen:cache', 'Invalid reference cache: missing %s', name);
    end
end

function mode = validate_cache_mode_(meta)
    if ~isfield(meta, 'schema') || ~strcmp(meta.schema, 'reference_cache/0.1')
        if isfield(meta, 'schema')
            error('statgen:cache', 'Unsupported reference cache schema: %s', meta.schema);
        end
        error('statgen:cache', 'Unsupported reference cache schema');
    end
    if ~isfield(meta, 'mode')
        error('statgen:cache', 'reference cache missing mode; delete and rebuild old cache');
    end
    mode = lower(char(meta.mode));
    if ~any(strcmp(mode, {'full', 'thin'}))
        error('statgen:cache', 'Invalid reference cache mode: %s', mode);
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
