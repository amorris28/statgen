function panel = load_reference_cache(path, varargin)
%LOAD_REFERENCE_CACHE Load a cached ReferencePanel.
%
%   panel = statgen.load_reference_cache(path)
%   panel = statgen.load_reference_cache(path, shards)
%
% Loads one MATLAB .mat reference cache file produced by
% statgen.save_reference_cache. Optional shards return a logical subset of the
% cached shard labels. String fields are stored as per-shard text payloads and
% decoded lazily when snp, a1, or a2 is accessed.
%
% See also statgen.ReferencePanel, statgen.load_reference,
% statgen.save_reference_cache.
    shards = parse_args_(varargin{:});

    path = char(path);
    meta_loaded = load(path, 'metadata');
    require_field_(meta_loaded, 'metadata');
    meta = meta_loaded.metadata;
    validate_cache_metadata_schema_(meta);
    loaded = load(path, 'bp', 'snp_text_by_shard', 'a1_text_by_shard', ...
        'a2_text_by_shard', 'a1_hash64', 'a2_hash64');

    require_field_(loaded, 'bp');
    require_field_(loaded, 'snp_text_by_shard');
    require_field_(loaded, 'a1_text_by_shard');
    require_field_(loaded, 'a2_text_by_shard');
    require_field_(loaded, 'a1_hash64');
    require_field_(loaded, 'a2_hash64');
    bp = loaded.bp(:);
    a1_hash64 = uint64(loaded.a1_hash64(:));
    a2_hash64 = uint64(loaded.a2_hash64(:));

    n = numel(bp);
    if numel(a1_hash64) ~= n || numel(a2_hash64) ~= n
        error('statgen:cache', 'Invalid reference cache: allele hash vector lengths mismatch');
    end

    [labels, checksums, start0, stop0] = validate_metadata_(meta, n);
    snp_text_by_shard = validate_text_payloads_(loaded.snp_text_by_shard, numel(labels), 'snp_text_by_shard');
    a1_text_by_shard = validate_text_payloads_(loaded.a1_text_by_shard, numel(labels), 'a1_text_by_shard');
    a2_text_by_shard = validate_text_payloads_(loaded.a2_text_by_shard, numel(labels), 'a2_text_by_shard');
    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_reference_cache');
    shard_objs = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);

        shard_obj = statgen.ReferenceShard.from_cache_text( ...
            labels{idx}, numel(ix), bp(ix), a1_hash64(ix), a2_hash64(ix), ...
            checksums{idx}, snp_text_by_shard{idx}, a1_text_by_shard{idx}, ...
            a2_text_by_shard{idx});
        shard_objs{i} = shard_obj;
    end

    panel = statgen.ReferencePanel(shard_objs);
end

function shards = parse_args_(varargin)
    shards = [];
    if numel(varargin) > 1
        error('statgen:arg', 'load_reference_cache accepts path and optional shards');
    end
    if ~isempty(varargin)
        shards = varargin{1};
    end
end

function require_field_(loaded, name)
    if ~isfield(loaded, name)
        if strcmp(name, 'a1_hash64') || strcmp(name, 'a2_hash64')
            error('statgen:cache', 'reference cache missing a1_hash64/a2_hash64; rebuild cache');
        end
        error('statgen:cache', 'Invalid reference cache: missing %s', name);
    end
end

function validate_cache_metadata_schema_(meta)
    if ~isfield(meta, 'schema') || ~strcmp(meta.schema, 'reference_cache/0.1')
        if isfield(meta, 'schema')
            error('statgen:cache', 'Unsupported reference cache schema: %s', meta.schema);
        end
        error('statgen:cache', 'Unsupported reference cache schema');
    end
end

function [labels, checksums, start0, stop0] = validate_metadata_(meta, num_snp)
    [labels, checksums, start0, stop0] = statgen.internal.validate_cache_shard_metadata( ...
        meta, 'reference_cache/0.1', 'reference', num_snp, false);
end

function out = validate_text_payloads_(value, n_shards, name)
    out = statgen.internal.ensure_cell_col(value);
    if numel(out) ~= n_shards
        error('statgen:cache', 'Invalid reference cache: %s length mismatch', name);
    end
    for i = 1:n_shards
        if ~(ischar(out{i}) || isstring(out{i}))
            error('statgen:cache', 'Invalid reference cache: %s{%d} must be text', name, i);
        end
        out{i} = char(out{i});
    end
end
