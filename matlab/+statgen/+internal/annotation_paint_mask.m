function mask = annotation_paint_mask(bp, intervals)
% Paint sorted non-overlapping intervals onto BIM base-pair coordinates.
    n = numel(bp);
    mask = zeros(n, 1);
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
    if any(valid)
        hit = valid & (pos0 < ends(max(idx, 1)));
        mask(hit) = 1;
    end
end
