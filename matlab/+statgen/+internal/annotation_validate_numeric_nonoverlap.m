function annotation_validate_numeric_nonoverlap(chr_col, starts, ends, path, row_base0)
% Validate numeric annotation intervals do not overlap within chromosomes.
    if nargin < 5
        row_base0 = 0;
    end
    labels = unique(chr_col, 'stable');
    for i = 1:numel(labels)
        label = labels{i};
        if strcmp(label, 'Y') || strcmp(label, 'MT')
            continue
        end
        rows = find(strcmp(chr_col, label));
        if numel(rows) <= 1
            continue
        end
        [~, order] = sortrows([double(starts(rows)), double(ends(rows))], [1, 2]);
        sorted_rows = rows(order);
        s = double(starts(sorted_rows));
        e = double(ends(sorted_rows));
        prev_max = cummax(e);
        overlap = s(2:end) < prev_max(1:end-1);
        if any(overlap)
            j = find(overlap, 1, 'first') + 1;
            error('statgen:annotations', ...
                '%s: row %d: numeric annotation intervals overlap on chromosome %s', ...
                path, row_base0 + sorted_rows(j), label);
        end
    end
end
