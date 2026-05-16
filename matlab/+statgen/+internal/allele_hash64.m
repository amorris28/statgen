function out = allele_hash64(alleles)
%ALLELE_HASH64 Deterministic allele hash defined by spec/reference.md.
    alleles = statgen.internal.ensure_cell_col(alleles);
    n = numel(alleles);
    out = zeros(n, 1, 'uint64');
    if n == 0
        return
    end

    p = 2147483647;
    base1 = 257;
    base2 = 263;
    max_chars = 150;

    lens = cellfun('length', alleles);
    too_long = lens > max_chars;
    if any(too_long)
        warning('statgen:allele_hash', ...
            'Allele length exceeds 150 characters; hashing uses first 150 characters');
        idx = find(too_long);
        for k = 1:numel(idx)
            i = idx(k);
            alleles{i} = alleles{i}(1:max_chars);
        end
        lens = min(lens, max_chars);
    end

    single_base = lens == 1;
    if any(single_base)
        bases = {'A'; 'C'; 'G'; 'T'};
        base_x = double(uint8(char(bases))) + 1;
        base_x = base_x(:);
        base_h1 = mod(base1 + base_x, p);
        base_h2 = mod(base2 + base_x, p);
        base_hash = bitshift(uint64(base_h1), 32) + uint64(base_h2);

        single_idx = find(single_base);
        single_alleles = alleles(single_idx);
        assigned = false(numel(single_idx), 1);
        for k = 1:numel(bases)
            matches = strcmp(single_alleles, bases{k});
            if any(matches)
                out(single_idx(matches)) = base_hash(k);
                assigned(matches) = true;
            end
        end
        if any(~assigned)
            other_idx = single_idx(~assigned);
            x = double(uint8(char(alleles(other_idx)))) + 1;
            x = x(:);
            h1 = mod(base1 + x, p);
            h2 = mod(base2 + x, p);
            out(other_idx) = bitshift(uint64(h1), 32) + uint64(h2);
        end
    end
    if all(single_base)
        return
    end

    multi_idx = find(~single_base);
    multi_alleles = alleles(multi_idx);
    multi_lens = lens(multi_idx);

    mat = uint8(char(multi_alleles));
    L = size(mat, 2);
    cols = 1:L;
    pad = bsxfun(@gt, cols, multi_lens);
    mat(pad) = uint8(0);

    h1 = ones(numel(multi_idx), 1);
    h2 = ones(numel(multi_idx), 1);
    for j = 1:L
        active = j <= multi_lens;
        x = double(mat(active, j)) + 1;
        h1(active) = mod(h1(active) .* base1 + x, p);
        h2(active) = mod(h2(active) .* base2 + x, p);
    end

    out(multi_idx) = bitshift(uint64(h1), 32) + uint64(h2);
end
