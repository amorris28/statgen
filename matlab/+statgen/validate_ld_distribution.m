function report = validate_ld_distribution(path, check_payload_structure)
%VALIDATE_LD_DISTRIBUTION Validate a MATLAB sparse LD distribution.
%
%   report = statgen.validate_ld_distribution(path)
%   report = statgen.validate_ld_distribution(path, check_payload_structure)
%
% Checks the MATLAB LD manifest, shard files, checksums, and bundled reference
% BIM files. Set check_payload_structure to true for deeper sparse-matrix
% payload checks.
% Metadata checks include num_monomorphic_snps; for forced monomorphic SNPs,
% undefined off-diagonal LD is represented by omitted sparse entries.
%
% See also statgen.load_ld, statgen.create_ld_mat_manifest.
    if nargin < 2 || isempty(check_payload_structure)
        check_payload_structure = false;
    end

    path = char(path);
    warnings_list = {};
    if ~isfolder(path)
        error('statgen:ld', 'validate_ld_distribution: path must identify a panel root directory');
    end

    manifest = statgen.internal.ld_read_manifest(fullfile(path, 'ld_manifest.json'), ...
        'matlab_mat_sparse_double');
    seen_reference_bim = {};
    for i = 1:numel(manifest.shards)
        entry = manifest.shards(i);
        shard_path = fullfile(path, entry.file);
        if exist(shard_path, 'file') ~= 2
            error('statgen:io', 'LD file not found: %s', shard_path);
        end
        if ~strcmp(statgen.internal.ld_md5_file(shard_path), entry.file_md5)
            error('statgen:ld', '%s: file_md5 does not match manifest', shard_path);
        end
        if statgen.internal.ld_is_mat_v5(shard_path)
            warnings_list{end+1} = statgen.internal.ld_v5_warning(); %#ok<AGROW>
        end
        [~, meta] = statgen.internal.ld_read_mat_shard(shard_path, check_payload_structure);
        statgen.internal.ld_validate_manifest_entry_agreement(entry, meta, shard_path);
        reference_key = sprintf('%s|%s|%d|%s', ...
            char(entry.reference_bim), char(entry.chr), double(entry.num_snp), char(entry.reference_checksum));
        if ~any(strcmp(seen_reference_bim, reference_key))
            validate_bundled_reference_(path, entry);
            seen_reference_bim{end+1} = reference_key; %#ok<AGROW>
        end
    end

    warnings_list = unique(warnings_list);
    for i = 1:numel(warnings_list)
        warning('statgen:ld:v5mat', '%s', warnings_list{i});
    end
    report = struct('ok', true, 'warnings', {warnings_list(:)});
end

function validate_bundled_reference_(root, entry)
    reference_path = fullfile(root, entry.reference_bim);
    if exist(reference_path, 'file') ~= 2
        error('statgen:io', 'LD file not found: %s', reference_path);
    end
    panel = statgen.load_reference(reference_path);
    if numel(panel.shards) ~= 1
        error('statgen:ld', '%s: bundled reference_bim must contain exactly one shard', reference_path);
    end
    ref = panel.shards{1};
    if ~strcmp(ref.label, entry.chr)
        error('statgen:ld', '%s: bundled reference_bim chr does not match manifest', reference_path);
    end
    if ref.num_snp ~= double(entry.num_snp)
        error('statgen:ld', '%s: bundled reference_bim num_snp does not match manifest', reference_path);
    end
    if ~strcmp(ref.checksum, entry.reference_checksum)
        error('statgen:ld', '%s: bundled reference_bim reference_checksum does not match manifest', reference_path);
    end
end
