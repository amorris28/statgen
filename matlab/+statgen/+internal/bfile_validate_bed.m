function bed_file_size = bfile_validate_bed(path, source_num_sample, source_num_snp)
%BFILE_VALIDATE_BED Validate PLINK BED magic, SNP-major mode, and exact size.
    path = char(path);
    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open BED file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    magic = fread(fid, 3, 'uint8=>uint8');
    clear cleaner;

    if numel(magic) ~= 3 || magic(1) ~= 108 || magic(2) ~= 27
        error('statgen:bed', '%s: invalid PLINK BED magic bytes', path);
    end
    if magic(3) ~= 1
        error('statgen:bed', '%s: PLINK BED must be SNP-major mode', path);
    end

    info = dir(path);
    if isempty(info)
        error('statgen:io', 'Cannot stat BED file: %s', path);
    end
    bed_file_size = double(info.bytes);
    expected = statgen.internal.bfile_expected_bed_size(source_num_sample, source_num_snp);
    if bed_file_size ~= expected
        error('statgen:bed', '%s: BED file size mismatch: expected %.0f, got %.0f', ...
            path, expected, bed_file_size);
    end
end
