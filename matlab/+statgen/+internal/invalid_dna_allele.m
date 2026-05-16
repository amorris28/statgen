function bad = invalid_dna_allele(alleles)
%INVALID_DNA_ALLELE True for alleles outside uppercase A/C/G/T strings.
    alleles = statgen.internal.ensure_cell_col(alleles);
    lens = cellfun('length', alleles);
    if isempty(alleles)
        bad = false(0, 1);
        return
    end
    chars = char(alleles);
    cols = 1:size(chars, 2);
    padding = bsxfun(@gt, cols, lens);
    invalid = chars ~= 'A' & chars ~= 'C' & chars ~= 'G' & chars ~= 'T';
    bad = any(invalid & ~padding, 2);
end
