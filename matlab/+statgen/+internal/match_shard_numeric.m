function loc = match_shard_numeric(ref_bp, ref_a1_hash64, ref_a2_hash64, src_bp, src_a1_hash64, src_a2_hash64, label, context)
%MATCH_SHARD_NUMERIC Align source rows to reference rows by numeric variant key.
%
% Returns a length(ref_bp) vector. loc(i) is the 1-based source row matching
% reference row i, or 0 when the reference row is absent from the source.
    if nargin < 8 || isempty(context)
        context = 'source';
    end

    n_ref = numel(ref_bp);
    loc = zeros(n_ref, 1);
    if isempty(src_bp)
        return
    end

    all_bp = [ref_bp(:); src_bp(:)];
    all_a1_hash64 = [uint64(ref_a1_hash64(:)); uint64(src_a1_hash64(:))];
    all_a2_hash64 = [uint64(ref_a2_hash64(:)); uint64(src_a2_hash64(:))];
    is_ref = [true(n_ref, 1); false(numel(src_bp), 1)];
    ref_index = [(1:n_ref)'; zeros(numel(src_bp), 1)];
    src_index = [zeros(n_ref, 1); (1:numel(src_bp))'];

    order = sort_key_order_(all_bp, all_a1_hash64, all_a2_hash64);
    bp_sorted = all_bp(order);
    a1_sorted = all_a1_hash64(order);
    a2_sorted = all_a2_hash64(order);

    same_prev = [false; ...
        bp_sorted(2:end) == bp_sorted(1:end-1) & ...
        a1_sorted(2:end) == a1_sorted(1:end-1) & ...
        a2_sorted(2:end) == a2_sorted(1:end-1)];
    group_id = cumsum(~same_prev);
    n_groups = group_id(end);

    is_ref_sorted = is_ref(order);
    ref_count = accumarray(group_id, double(is_ref_sorted), [n_groups, 1], @sum, 0);
    src_count = accumarray(group_id, double(~is_ref_sorted), [n_groups, 1], @sum, 0);

    ambiguous = ref_count > 1 | src_count > 1;
    if any(ambiguous)
        g = find(ambiguous, 1, 'first');
        pos = find(group_id == g, 1, 'first');
        error('statgen:match', ...
            'Ambiguous duplicate %s/reference matching key in shard %s: bp=%d, a1_hash64=%s, a2_hash64=%s', ...
            char(context), char(label), round(bp_sorted(pos)), uint64_to_string_(a1_sorted(pos)), uint64_to_string_(a2_sorted(pos)));
    end

    ref_by_group = accumarray(group_id, ref_index(order), [n_groups, 1], @max, 0);
    src_by_group = accumarray(group_id, src_index(order), [n_groups, 1], @max, 0);
    matched = ref_count == 1 & src_count == 1;
    loc(ref_by_group(matched)) = src_by_group(matched);
end

function order = sort_key_order_(bp, a1_hash64, a2_hash64)
    order = (1:numel(bp))';
    if isempty(order)
        return
    end
    order = stable_sort_by_(order, a2_hash64);
    order = stable_sort_by_(order, a1_hash64);
    order = stable_sort_by_(order, bp);
end

function order = stable_sort_by_(order, values)
    current_rank = (1:numel(order))';
    sorted_values = values(order);
    [sorted_values, p] = sort(sorted_values, 'ascend');
    order = order(p);
    sorted_rank = current_rank(p);

    same_prev = [false; sorted_values(2:end) == sorted_values(1:end-1)];
    group_id = cumsum(~same_prev);
    [~, q] = sortrows([double(group_id(:)), double(sorted_rank(:))], [1, 2]);
    order = order(q);
end

function out = uint64_to_string_(x)
    out = dec2hex(uint64(x), 16);
end
