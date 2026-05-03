function panel = load_ld(path, reference, default_chrX_sex)
% Load MATLAB/Octave-native sparse LD distribution artifacts.
    if nargin < 3 || isempty(default_chrX_sex)
        default_chrX_sex = 'female';
    end
    statgen.LDPanel.validate_chrx_sex_(default_chrX_sex, 'default_chrX_sex');

    path = char(path);
    if isfolder(path)
        manifest = statgen.internal.ld_read_manifest(fullfile(path, 'ld_manifest.json'), ...
            'matlab_mat_sparse_double');
        groups = load_panel_root_(path, manifest, reference.shards);
    else
        groups = load_single_shard_(path, reference.shards);
    end

    panel = statgen.LDPanel(groups, default_chrX_sex);
end

function groups = load_panel_root_(root, manifest, ref_shards)
    entries = manifest.shards;
    groups = cell(numel(ref_shards), 1);

    for i = 1:numel(ref_shards)
        ref = ref_shards{i};
        selected = [];
        for j = 1:numel(entries)
            entry = entries(j);
            if strcmp(ref.label, 'X')
                if strcmp(entry.chr, 'X')
                    selected = [selected; j]; %#ok<AGROW>
                end
            else
                if strcmp(entry.chr, ref.label) && isempty(entry.sex)
                    selected = [selected; j]; %#ok<AGROW>
                end
            end
        end

        if strcmp(ref.label, 'X')
            if isempty(selected)
                error('statgen:ld', 'LD manifest has no chrX shard for reference shard X');
            end
        elseif numel(selected) ~= 1
            error('statgen:ld', ...
                'LD manifest must contain exactly one autosomal shard for chr %s', ref.label);
        end

        group = cell(numel(selected), 1);
        seen_sex = {};
        for k = 1:numel(selected)
            entry = entries(selected(k));
            shard_path = fullfile(root, entry.file);
            [shard, meta] = statgen.internal.ld_read_mat_shard(shard_path, false);
            statgen.internal.ld_validate_manifest_entry_agreement(entry, meta, shard_path);
            statgen.internal.ld_validate_reference_compatibility(shard, ref, shard_path);
            sex_key = shard.sex;
            if any(strcmp(seen_sex, sex_key))
                error('statgen:ld', 'LD manifest has duplicate shard for chr %s', shard.chr);
            end
            seen_sex{end+1} = sex_key; %#ok<AGROW>
            group{k} = shard;
        end
        groups{i} = group;
    end
end

function groups = load_single_shard_(path, ref_shards)
    [shard, ~] = statgen.internal.ld_read_mat_shard(path, false);
    if numel(ref_shards) ~= 1 || ~strcmp(ref_shards{1}.label, shard.chr)
        error('statgen:ld', ...
            'single-shard LD loads require a single-shard reference with the same chromosome');
    end
    statgen.internal.ld_validate_reference_compatibility(shard, ref_shards{1}, path);
    groups = {{shard}};
end
