function [is_single_nucleotide_variant, is_strand_ambiguous] = reference_hash_masks(a1_hash64, a2_hash64)
%REFERENCE_HASH_MASKS Variant-type masks from reference allele hashes.
    a1 = uint64(a1_hash64(:));
    a2 = uint64(a2_hash64(:));

    base_hashes = statgen.internal.allele_hash64({'A'; 'C'; 'G'; 'T'});
    is_single_nucleotide_variant = ismember(a1, base_hashes) & ismember(a2, base_hashes);

    a_hash = base_hashes(1);
    c_hash = base_hashes(2);
    g_hash = base_hashes(3);
    t_hash = base_hashes(4);
    is_strand_ambiguous = ...
        (a1 == a_hash & a2 == t_hash) | ...
        (a1 == t_hash & a2 == a_hash) | ...
        (a1 == c_hash & a2 == g_hash) | ...
        (a1 == g_hash & a2 == c_hash);
end
