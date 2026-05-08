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
    loaded = load(path, 'metadata', 'annomat', 'annonames');
    meta = loaded.metadata;
    annomat = sparse(loaded.annomat);
    annonames = statgen.internal.ensure_cell_col(loaded.annonames);

    [labels, checksums, start0, stop0] = validate_metadata_(meta, size(annomat, 1));
    validate_annonames_(annonames);

    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_annotations_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);
        out_shards{i} = statgen.AnnotationShard(labels{idx}, checksums{idx}, annomat(ix, :), false);
    end

    panel = statgen.AnnotationPanel(out_shards, annonames);
end

function [labels, checksums, start0, stop0] = validate_metadata_(meta, num_snp)
    if ~strcmp(meta.schema, 'annotations_cache/0.1')
        error('statgen:cache', 'Unsupported annotations cache schema: %s', meta.schema);
    end

    labels = statgen.internal.ensure_cell_col(meta.shard_labels);
    checksums = statgen.internal.ensure_cell_col(meta.shard_checksums);
    start0 = double(meta.shard_start0(:));
    stop0 = double(meta.shard_stop0(:));

    n_shards = numel(labels);
    if ~isfield(meta, 'n_shards') || ~isscalar(meta.n_shards) || meta.n_shards ~= n_shards
        error('statgen:cache', 'Invalid annotations cache: n_shards mismatch');
    end
    if numel(checksums) ~= n_shards || numel(start0) ~= n_shards || numel(stop0) ~= n_shards
        error('statgen:cache', 'Invalid annotations cache: shard metadata length mismatch');
    end
    validate_offsets_(start0, stop0, num_snp);
end

function validate_offsets_(start0, stop0, num_snp)
    if isempty(start0)
        if num_snp ~= 0
            error('statgen:cache', 'Invalid annotations cache: empty shard offsets for non-empty payload');
        end
        return
    end
    if start0(1) ~= 0 || stop0(end) ~= num_snp || any(stop0 < start0) || any(start0(2:end) ~= stop0(1:end-1))
        error('statgen:cache', 'Invalid annotations cache: shard offsets are not contiguous');
    end
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
