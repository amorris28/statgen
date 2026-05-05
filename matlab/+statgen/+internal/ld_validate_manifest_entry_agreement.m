function ld_validate_manifest_entry_agreement(entry, meta, path)
    if ~strcmp(entry.chr, meta.chr)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for chr', path);
    end
    if ~statgen.internal.ld_same_optional_string(entry.sex, meta.sex)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for sex', path);
    end
    if double(entry.num_snp) ~= double(meta.num_snp)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for num_snp', path);
    end
    if double(entry.nnz) ~= double(meta.nnz)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for nnz', path);
    end
    if ~strcmp(entry.reference_checksum, meta.reference_checksum)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for reference_checksum', path);
    end
    if ~strcmp(entry.reference_bim, meta.reference_bim)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for reference_bim', path);
    end
end
