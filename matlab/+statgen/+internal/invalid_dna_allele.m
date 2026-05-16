function [bad, lens] = invalid_dna_allele(alleles)
%INVALID_DNA_ALLELE True for alleles outside uppercase A/C/G/T strings.
    alleles = statgen.internal.ensure_cell_col(alleles);
    lens = cellfun('length', alleles);
    if isempty(alleles)
        bad = false(0, 1);
        return
    end

    bad = lens == 0;

    single = lens == 1;
    if any(single)
        chars = char(alleles(single));
        bad(single) = chars ~= 'A' & chars ~= 'C' & chars ~= 'G' & chars ~= 'T';
    end

    multi = lens > 1;
    if any(multi)
        chars = char(alleles(multi));
        multi_lens = lens(multi);
        cols = 1:size(chars, 2);
        padding = bsxfun(@gt, cols, multi_lens);
        invalid = chars ~= 'A' & chars ~= 'C' & chars ~= 'G' & chars ~= 'T';
        bad(multi) = any(invalid & ~padding, 2);
    end
end
