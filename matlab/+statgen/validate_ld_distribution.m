function report = validate_ld_distribution(path, check_payload_structure)
% Validate a MATLAB/Octave sparse LD distribution.
    if nargin < 2 || isempty(check_payload_structure)
        check_payload_structure = false;
    end

    path = char(path);
    warnings_list = {};
    if isfolder(path)
        manifest = statgen.internal.ld_read_manifest(fullfile(path, 'ld_manifest.json'), ...
            'matlab_mat_sparse_double');
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
        end
    else
        if statgen.internal.ld_is_mat_v5(path)
            warnings_list{end+1} = statgen.internal.ld_v5_warning(); %#ok<AGROW>
        end
        statgen.internal.ld_read_mat_shard(path, check_payload_structure);
    end

    warnings_list = unique(warnings_list);
    for i = 1:numel(warnings_list)
        warning('statgen:ld:v5mat', '%s', warnings_list{i});
    end
    report = struct('ok', true, 'warnings', {warnings_list(:)});
end
