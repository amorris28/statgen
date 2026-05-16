function keep = validate_variant_chr_labels(labels, path, line_offset, error_id)
%VALIDATE_VARIANT_CHR_LABELS Validate canonical source contig labels.
%
% Returns a logical mask for rows on supported canonical contigs 1-22 and X.
% Recognized non-supported Y/MT rows are valid labels but are not kept.
    labels = statgen.internal.ensure_cell_col(labels);
    if nargin < 3
        line_offset = 0;
    end
    if nargin < 4
        error_id = 'statgen:variant';
    end
    path = char(path);

    [unique_labels, ~, label_idx] = unique(labels);

    bad_chr_unique = cellfun('isempty', unique_labels);
    bad_chr = bad_chr_unique(label_idx);
    if any(bad_chr)
        idx = find(bad_chr, 1, 'first');
        error(error_id, '%s:%d: chr must be non-empty', path, idx + line_offset);
    end

    chr_style_unique = cellfun(@(c) strncmpi(c, 'chr', 3), unique_labels);
    chr_style = chr_style_unique(label_idx);
    if any(chr_style)
        idx = find(chr_style, 1, 'first');
        error(error_id, '%s:%d: chr-style labels (e.g., chr1/chrX) are not allowed', path, idx + line_offset);
    end

    canonical = statgen.internal.canonical_labels();
    keep_unique = ismember(unique_labels, canonical);
    known_unique = keep_unique | ismember(unique_labels, {'Y', 'MT'});
    keep = keep_unique(label_idx);
    known = known_unique(label_idx);
    if ~all(known)
        idx = find(~known, 1, 'first');
        error(error_id, ...
            '%s:%d: unsupported chr label %s; expected 1-22, X (Y/MT are ignored)', ...
            path, idx + line_offset, labels{idx});
    end
end
