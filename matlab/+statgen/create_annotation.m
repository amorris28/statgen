function panel = create_annotation(reference, annovec, annotation_name, is_binary, annotation_metadata)
%CREATE_ANNOTATION Create one annotation from an aligned numeric vector.
%
%   annotations = statgen.create_annotation(reference, annovec, annotation_name)
%   annotations = statgen.create_annotation(..., is_binary, annotation_metadata)
%
% annovec must have num_snp elements already aligned to reference. The output
% is an AnnotationPanel with one annotation column named annotation_name.
%
% See also statgen.AnnotationPanel, statgen.create_annotations.
    if nargin < 3
        error('statgen:arg', 'create_annotation requires reference, annovec, annotation_name');
    end
    if nargin < 4
        is_binary = [];
    end
    if nargin < 5
        annotation_metadata = '';
    end

    name = char(annotation_name);
    if isempty(name)
        error('statgen:annotations', 'annotation_name must be non-empty');
    end

    n = reference.num_snp;
    vec = annovec;
    if ndims(vec) > 2 || (ismatrix(vec) && ~any(size(vec) == 1))
        error('statgen:annotations', 'annovec must be a vector with length %d', n);
    end
    vec = double(vec(:));
    if numel(vec) ~= n
        error('statgen:annotations', 'annovec length mismatch: expected %d, got %d', n, numel(vec));
    end

    bad = ~isfinite(vec);
    if any(bad)
        i = find(bad, 1, 'first');
        error('statgen:annotations', 'annovec(%d) must be finite numeric', i);
    end

    if ~isempty(is_binary)
        if ~isscalar(is_binary) || (ischar(is_binary) || isstring(is_binary))
            error('statgen:annotations', 'is_binary must be a scalar logical for create_annotation');
        end
        binary = logical(is_binary);
    else
        binary = [];
    end

    panel = statgen.create_annotations(reference, sparse(vec), {name}, binary, {char(annotation_metadata)});
end
