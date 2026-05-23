function merged = annotation_merge_intervals(starts, ends)
% Merge overlapping or adjacent BED intervals.
    if isempty(starts)
        merged = zeros(0, 2);
        return
    end

    pairs = sortrows([double(starts(:)), double(ends(:))], [1, 2]);
    s = pairs(:, 1);
    e = pairs(:, 2);

    running_max_end = cummax(e);
    new_group = [true; s(2:end) > running_max_end(1:end-1)];

    merged_starts = s(new_group);
    group_id = cumsum(new_group);
    merged_ends = accumarray(group_id, e, [], @max);

    merged = [merged_starts, merged_ends];
end
