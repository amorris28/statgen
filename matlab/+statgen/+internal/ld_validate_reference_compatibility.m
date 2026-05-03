function ld_validate_reference_compatibility(shard, ref_shard, path)
    if ~strcmp(shard.chr, ref_shard.label)
        error('statgen:ld', '%s: LD chr does not match reference shard', path);
    end
    if shard.num_snp ~= ref_shard.num_snp
        error('statgen:ld', '%s: LD num_snp does not match reference shard', path);
    end
    if ~strcmp(shard.reference_checksum, ref_shard.checksum)
        error('statgen:ld', '%s: LD reference_checksum does not match reference shard', path);
    end
end
