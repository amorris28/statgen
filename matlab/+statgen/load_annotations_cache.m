function panel = load_annotations_cache(path, shards)
%LOAD_ANNOTATIONS_CACHE Load cached annotations.
%
%   annotations = statgen.load_annotations_cache(path)
%   annotations = statgen.load_annotations_cache(path, shards)
%
% Loads one MATLAB .mat AnnotationPanel cache file produced by
% statgen.save_annotations_cache. Optional shards return a logical subset of
% the cached shard labels.
%
% See also statgen.AnnotationPanel, statgen.load_annotations,
% statgen.save_annotations_cache.
    if nargin < 2
        shards = [];
    end

    path = char(path);
    loaded = load(path, 'metadata', 'annomat', 'annonames', 'is_binary', 'annotation_metadata');
    meta = loaded.metadata;
    annomat = sparse(double(loaded.annomat));
    annonames = statgen.internal.ensure_cell_col(loaded.annonames);

    [labels, checksums, start0, stop0, schema] = validate_metadata_(meta, size(annomat, 1));
    validate_annonames_(annonames);
    if strcmp(schema, 'annotations_cache/0.1')
        is_binary = true(numel(annonames), 1);
        annotation_metadata = repmat({''}, numel(annonames), 1);
    else
        if ~isfield(loaded, 'is_binary') || ~isfield(loaded, 'annotation_metadata')
            error('statgen:cache', 'Invalid annotations cache: missing is_binary or annotation_metadata');
        end
        is_binary = validate_is_binary_(loaded.is_binary, numel(annonames));
        annotation_metadata = validate_annotation_metadata_(loaded.annotation_metadata, numel(annonames));
    end

    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_annotations_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);
        out_shards{i} = statgen.AnnotationShard(labels{idx}, checksums{idx}, annomat(ix, :));
    end

    panel = statgen.AnnotationPanel(out_shards, annonames, is_binary, annotation_metadata);
end

function [labels, checksums, start0, stop0, schema] = validate_metadata_(meta, num_snp)
    if ~isfield(meta, 'schema')
        error('statgen:cache', 'Unsupported annotations cache schema');
    end
    schema = char(meta.schema);
    if ~(strcmp(schema, 'annotations_cache/0.1') || strcmp(schema, 'annotations_cache/0.2'))
        error('statgen:cache', 'Unsupported annotations cache schema: %s', schema);
    end
    [labels, checksums, start0, stop0] = statgen.internal.validate_cache_shard_metadata( ...
        meta, schema, 'annotations', num_snp, false);
end

function validate_annonames_(annonames)
    if isempty(annonames)
        error('statgen:cache', 'Invalid annotations cache: annonames must be non-empty');
    end
    if any(cellfun('isempty', annonames))
        error('statgen:cache', 'Invalid annotations cache: annonames must not contain empty names');
    end
    if numel(unique(annonames)) ~= numel(annonames)
        error('statgen:cache', 'Invalid annotations cache: annonames must be unique');
    end
end

function out = validate_is_binary_(value, expected_len)
    if ischar(value) || isstring(value)
        error('statgen:cache', 'Invalid annotations cache: is_binary must be logical');
    end
    arr = value(:);
    if numel(arr) ~= expected_len
        error('statgen:cache', 'Invalid annotations cache: is_binary length mismatch');
    end
    if ~islogical(arr)
        bad = ~(arr == 0 | arr == 1);
        if any(bad)
            error('statgen:cache', 'Invalid annotations cache: is_binary must be logical');
        end
    end
    out = logical(arr);
end

function out = validate_annotation_metadata_(value, expected_len)
    if ischar(value)
        error('statgen:cache', 'Invalid annotations cache: annotation_metadata must be a cell array or string vector');
    elseif isstring(value)
        value = cellstr(value(:));
    end
    out = statgen.internal.ensure_cell_col(value);
    if numel(out) ~= expected_len
        error('statgen:cache', 'Invalid annotations cache: annotation_metadata length mismatch');
    end
end
