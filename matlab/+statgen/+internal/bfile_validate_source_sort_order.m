function bfile_validate_source_sort_order(bim, path, line_idx)
%BFILE_VALIDATE_SOURCE_SORT_ORDER Validate canonical BIM row order and keys.
    if nargin < 3 || isempty(line_idx)
        line_idx = (1:numel(bim.chr))';
    else
        line_idx = line_idx(:);
    end

    if isempty(bim.chr)
        return
    end

    canonical = statgen.internal.canonical_labels();
    [is_ok_chr, chr_rank] = ismember(bim.chr, canonical);
    bad_chr = ~is_ok_chr;
    if any(bad_chr)
        bad_i = find(bad_chr, 1, 'first');
        error('statgen:bim', '%s:%d: chr must use canonical labels 1-22 or X', ...
            path, line_idx(bad_i));
    end

    n = numel(bim.chr);
    if n <= 1
        return
    end

    bad_order = (chr_rank(2:end) < chr_rank(1:end-1)) | ...
        ((chr_rank(2:end) == chr_rank(1:end-1)) & (bim.bp(2:end) < bim.bp(1:end-1)));
    bad = find(bad_order, 1, 'first') + 1;
    if ~isempty(bad)
        bad_line = line_idx(bad);
        error('statgen:bim', ...
            '%s:%d: rows must be sorted by (chr_rank, bp) in canonical contig order', ...
            path, bad_line);
    end

    a1_hash64 = statgen.internal.allele_hash64(bim.a1);
    a2_hash64 = statgen.internal.allele_hash64(bim.a2);
    keys = [uint64(chr_rank(:)), uint64(round(bim.bp(:))), a1_hash64(:), a2_hash64(:)];
    [sorted_keys, sort_idx] = sortrows(keys);
    dup_sorted = all(sorted_keys(2:end, :) == sorted_keys(1:end-1, :), 2);
    if any(dup_sorted)
        dup_pos = find(dup_sorted, 1, 'first') + 1;
        bad_line = line_idx(sort_idx(dup_pos));
        error('statgen:bim', ...
            '%s:%d: duplicate (chr, bp, a1_hash64, a2_hash64) matching key is not allowed', ...
            path, bad_line);
    end
end
