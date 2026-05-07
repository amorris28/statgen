function n_bytes = bfile_expected_bed_size(source_num_sample, source_num_snp)
%BFILE_EXPECTED_BED_SIZE Expected byte count for SNP-major PLINK BED.
    n_bytes = 3 + ceil(double(source_num_sample) / 4) * double(source_num_snp);
end
