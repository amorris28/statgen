function out = annotation_paint_numeric_sparse(bp, intervals, values)
% Paint interval value rows onto BIM base-pair coordinates as a sparse block.
    n = numel(bp);
    k = size(values, 2);
    out = sparse(n, k);
    if isempty(intervals)
        return
    end

    starts = double(intervals(:, 1));
    ends = double(intervals(:, 2));
    pos0 = double(bp(:)) - 1;

    % PERF: histc is retained for Octave compatibility; replace with
    % discretize(pos0, [starts; inf]) once the minimum Octave version supports it.
    [~, idx] = histc(pos0, [starts; inf]);
    valid = idx > 0;
    if ~any(valid)
        return
    end

    hit = valid & (pos0 < ends(max(idx, 1)));
    hit_rows = find(hit);
    if isempty(hit_rows)
        return
    end

    hit_values = values(idx(hit_rows), :);
    [row_idx, col_idx, vals] = find(hit_values);
    if isempty(vals)
        return
    end

    out = sparse(hit_rows(row_idx), col_idx, vals, n, k);
end
