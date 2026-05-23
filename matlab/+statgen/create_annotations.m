function panel = create_annotations(reference, annotation_matrix, annotation_names, is_binary, annotation_metadata)
%CREATE_ANNOTATIONS Create annotations from an aligned numeric matrix.
%
%   annotations = statgen.create_annotations(reference, annotation_matrix, annotation_names)
%   annotations = statgen.create_annotations(..., is_binary, annotation_metadata)
%
% annotation_matrix must be a finite num_snp-by-num_annot matrix already
% aligned to reference. is_binary marks columns that must contain only 0/1
% values; when omitted it is inferred from the observed values.
%
% See also statgen.AnnotationPanel, statgen.create_annotation,
% statgen.load_annotations, statgen.load_annotation.
    if nargin < 4
        is_binary = [];
    end
    if nargin < 5
        annotation_metadata = [];
    end
    names = ensure_names_(annotation_names);

    if issparse(annotation_matrix)
        A = sparse(double(annotation_matrix));
    else
        A = sparse(double(annotation_matrix));
    end
    if ndims(A) ~= 2
        error('statgen:annotations', 'annotation_matrix must be a 2D matrix');
    end

    n = reference.num_snp;
    k = numel(names);
    if size(A, 1) ~= n || size(A, 2) ~= k
        error('statgen:annotations', ...
            'annotation_matrix shape mismatch: expected (%d, %d), got (%d, %d)', ...
            n, k, size(A, 1), size(A, 2));
    end

    vals = nonzeros(A);
    if any(~isfinite(vals))
        v = vals(find(~isfinite(vals), 1, 'first'));
        error('statgen:annotations', 'annotation_matrix contains non-finite value: %g', v);
    end

    binary = coerce_or_infer_is_binary_(is_binary, A, k);
    statgen.internal.annotation_validate_declared_binary(A, binary, 'annotation_matrix');
    metadata = coerce_metadata_(annotation_metadata, k);

    shards = cell(numel(reference.shards), 1);
    for i = 1:numel(reference.shards)
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        shards{i} = statgen.AnnotationShard(s_ref.label, s_ref.checksum, A(ix, :));
    end

    panel = statgen.AnnotationPanel(shards, names, binary, metadata);
end

function out = ensure_names_(names)
    if ischar(names) || isstring(names)
        names = cellstr(names(:));
    end
    if ~iscell(names)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end

    out = statgen.internal.ensure_cell_col(names);
    if isempty(out)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end
    if any(cellfun('isempty', out))
        error('statgen:annotations', 'annotation names must not contain empty strings');
    end
    if numel(unique(out)) ~= numel(out)
        error('statgen:annotations', 'annotation names must be unique');
    end
end

function out = coerce_or_infer_is_binary_(value, A, k)
    if isempty(value)
        [~, cols, vals] = find(A);
        out = true(k, 1);
        out(unique(cols(vals ~= 1))) = false;
        return
    end
    if ischar(value) || isstring(value)
        error('statgen:annotations', 'is_binary must be a logical vector with length %d', k);
    end
    arr = value(:);
    if numel(arr) ~= k
        error('statgen:annotations', 'is_binary length mismatch: expected %d, got %d', k, numel(arr));
    end
    if ~islogical(arr)
        bad = ~(arr == 0 | arr == 1);
        if any(bad)
            i = find(bad, 1, 'first');
            error('statgen:annotations', 'is_binary(%d) must be logical', i);
        end
    end
    out = logical(arr);
end

function out = coerce_metadata_(value, k)
    if isempty(value)
        out = repmat({''}, k, 1);
        return
    end
    out = statgen.internal.annotation_coerce_metadata_vector(value, k, 'annotation_metadata');
end
