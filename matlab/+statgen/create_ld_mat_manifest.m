function manifest = create_ld_mat_manifest(input_root, output_root, shards)
%CREATE_LD_MAT_MANIFEST Create a MATLAB LD manifest for converted shards.
%
%   manifest = statgen.create_ld_mat_manifest(input_root, output_root, shards)
%
% input_root is the source Python LD distribution. output_root contains .mat
% shards created by statgen.convert_ld_npz_to_mat. shards is a non-empty list
% of canonical shard labels to include in the MATLAB LD distribution. Including
% 'X' includes all converted chrX sex-label shards present in output_root.
%
% See also statgen.convert_ld_npz_to_mat, statgen.validate_ld_distribution,
% statgen.load_ld.
    if nargin < 3 || isempty(shards)
        error('statgen:ld', 'create_ld_mat_manifest requires a non-empty shards list');
    end
    input_root = char(input_root);
    output_root = char(output_root);

    if exist(output_root, 'dir') ~= 7
        error('statgen:io', 'LD output directory not found: %s', output_root);
    end

    input_manifest = statgen.internal.ld_read_manifest( ...
        fullfile(input_root, 'ld_manifest.json'), 'python_npz_csc32');
    selected = validate_requested_manifest_shards_(input_manifest, shards);
    expected = expected_entries_(input_manifest, selected);

    records = struct('chr', {}, 'sex', {}, 'file', {}, 'file_md5', {}, ...
        'num_snp', {}, 'nnz', {}, 'reference_checksum', {}, 'reference_bim', {});
    for i = 1:numel(expected)
        entry = expected(i);
        mat_file = npz_file_to_mat_file_(entry.file);
        mat_path = fullfile(output_root, mat_file);
        if exist(mat_path, 'file') ~= 2
            error('statgen:io', 'Expected MATLAB LD shard not found: %s', mat_path);
        end

        loaded = load(mat_path, 'metadata');
        if ~isfield(loaded, 'metadata')
            error('statgen:ld', '%s: missing required MAT metadata variable', mat_path);
        end
        meta = loaded.metadata;
        statgen.internal.ld_validate_mat_metadata(meta, mat_path);
        validate_expected_metadata_(entry, meta, mat_path);

        reference_path = fullfile(output_root, meta.reference_bim);
        if exist(reference_path, 'file') ~= 2
            error('statgen:io', 'Expected bundled reference_bim not found: %s', reference_path);
        end

        records(end+1).chr = char(meta.chr); %#ok<AGROW>
        records(end).sex = meta.sex;
        records(end).file = mat_file;
        records(end).file_md5 = statgen.internal.ld_md5_file(mat_path);
        records(end).num_snp = double(meta.num_snp);
        records(end).nnz = double(meta.nnz);
        records(end).reference_checksum = char(meta.reference_checksum);
        records(end).reference_bim = char(meta.reference_bim);
    end

    manifest = struct();
    manifest.object_type = 'ld_panel_manifest';
    manifest.schema_version = '1.0';
    manifest.runtime_format = 'matlab_mat_sparse_double';
    manifest.shards = records;
    write_json_atomic_(fullfile(output_root, 'ld_manifest.json'), manifest);
end

function selected = validate_requested_manifest_shards_(manifest, shards)
    if ischar(shards) || isstring(shards)
        shards = cellstr(shards);
    elseif isnumeric(shards)
        error('statgen:shards', ...
            'create_ld_mat_manifest: shards must be a non-empty list of unique canonical contig labels');
    end
    shards = shards(:)';
    if isempty(shards)
        error('statgen:shards', ...
            'create_ld_mat_manifest: shards must be a non-empty list of unique canonical contig labels');
    end

    canonical = statgen.internal.canonical_labels();
    canonical_idx = zeros(1, numel(shards));
    selected = cell(1, numel(shards));
    seen = {};
    available = manifest_chr_labels_(manifest);
    for i = 1:numel(shards)
        label = char(shards{i});
        idx = find(strcmp(canonical, label), 1, 'first');
        if isempty(idx)
            error('statgen:shards', ...
                'create_ld_mat_manifest: unsupported shard label %s; expected canonical labels 1-22 or X', ...
                label);
        end
        if any(strcmp(seen, label))
            error('statgen:shards', 'create_ld_mat_manifest: duplicate shard label %s in shards', label);
        end
        if ~any(strcmp(available, label))
            error('statgen:ld', 'create_ld_mat_manifest: shard %s not found in Python manifest', label);
        end
        seen{end+1} = label; %#ok<AGROW>
        selected{i} = label;
        canonical_idx(i) = idx;
    end
    if any(diff(canonical_idx) <= 0)
        error('statgen:shards', 'create_ld_mat_manifest: shards must be in canonical subsequence order');
    end
end

function labels = manifest_chr_labels_(manifest)
    labels = {};
    for i = 1:numel(manifest.shards)
        label = char(manifest.shards(i).chr);
        if ~any(strcmp(labels, label))
            labels{end+1} = label; %#ok<AGROW>
        end
    end
end

function out = expected_entries_(manifest, selected)
    out = struct('chr', {}, 'sex', {}, 'file', {}, 'file_md5', {}, ...
        'num_snp', {}, 'nnz', {}, 'reference_checksum', {}, 'reference_bim', {});
    for i = 1:numel(selected)
        label = selected{i};
        idx = arrayfun(@(entry) strcmp(char(entry.chr), label), manifest.shards);
        matches = manifest.shards(idx);
        for j = 1:numel(matches)
            out(end+1) = matches(j); %#ok<AGROW>
        end
    end
end

function mat_file = npz_file_to_mat_file_(npz_file)
    npz_file = char(npz_file);
    if numel(npz_file) < 4 || ~strcmp(npz_file(end-3:end), '.npz')
        error('statgen:ld', 'Python LD manifest file must end with .npz: %s', npz_file);
    end
    mat_file = [npz_file(1:end-4) '.mat'];
end

function validate_expected_metadata_(entry, meta, mat_path)
    if ~strcmp(char(entry.chr), char(meta.chr))
        error('statgen:ld', '%s: expected .mat metadata mismatch for chr', mat_path);
    end
    if ~statgen.internal.ld_same_optional_string(entry.sex, meta.sex)
        error('statgen:ld', '%s: expected .mat metadata mismatch for sex', mat_path);
    end
    if double(entry.num_snp) ~= double(meta.num_snp)
        error('statgen:ld', '%s: expected .mat metadata mismatch for num_snp', mat_path);
    end
    if double(entry.nnz) ~= double(meta.nnz)
        error('statgen:ld', '%s: expected .mat metadata mismatch for nnz', mat_path);
    end
    if ~strcmp(char(entry.reference_checksum), char(meta.reference_checksum))
        error('statgen:ld', '%s: expected .mat metadata mismatch for reference_checksum', mat_path);
    end
    if ~strcmp(char(entry.reference_bim), char(meta.reference_bim))
        error('statgen:ld', '%s: expected .mat metadata mismatch for reference_bim', mat_path);
    end
end

function write_json_atomic_(path, obj)
    tmp = [path '.tmp'];
    fid = fopen(tmp, 'w');
    if fid < 0
        error('statgen:io', 'Cannot write manifest: %s', tmp);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', jsonencode(obj));
    clear cleaner;
    [ok, msg] = movefile(tmp, path, 'f');
    if ~ok
        error('statgen:io', 'Cannot replace manifest: %s', msg);
    end
end
