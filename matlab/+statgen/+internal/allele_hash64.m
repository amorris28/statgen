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

    mat = uint8(char(alleles));
    L = size(mat, 2);
    cols = 1:L;
    pad = bsxfun(@gt, cols, lens);
    mat(pad) = uint8(0);

    h1 = ones(n, 1);
    h2 = ones(n, 1);
    for j = 1:L
        active = j <= lens;
        x = double(mat(active, j)) + 1;
        h1(active) = mod(h1(active) .* base1 + x, p);
        h2(active) = mod(h2(active) .* base2 + x, p);
    end

    out = bitshift(uint64(h1), 32) + uint64(h2);
end
