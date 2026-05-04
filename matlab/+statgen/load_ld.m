function panel = load_ld(path, reference, default_chrX_sex, retain_ld_r)
% Load MATLAB/Octave-native sparse LD distribution artifacts.
    if nargin < 3 || isempty(default_chrX_sex)
        default_chrX_sex = 'female';
    end
    if nargin < 4 || isempty(retain_ld_r)
        retain_ld_r = true;
    end
    statgen.LDPanel.validate_chrx_sex_(default_chrX_sex, 'default_chrX_sex');
    if ~islogical(retain_ld_r) && ~(isnumeric(retain_ld_r) && isscalar(retain_ld_r))
        error('statgen:ld', 'retain_ld_r must be a scalar logical value');
    end
    retain_ld_r = logical(retain_ld_r);

    path = char(path);
    if isfolder(path)
        manifest = statgen.internal.ld_read_manifest(fullfile(path, 'ld_manifest.json'), ...
            'matlab_mat_sparse_double');
        groups = load_panel_root_(path, manifest, reference.shards, retain_ld_r);
    else
        groups = load_single_shard_(path, reference.shards, retain_ld_r);
    end

    panel = statgen.LDPanel(groups, default_chrX_sex);
end

function groups = load_panel_root_(root, manifest, ref_shards, retain_ld_r)
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
            [shard, meta] = statgen.internal.ld_read_mat_shard(shard_path, false, retain_ld_r);
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

function groups = load_single_shard_(path, ref_shards, retain_ld_r)
    [shard, ~] = statgen.internal.ld_read_mat_shard(path, false, retain_ld_r);
    if numel(ref_shards) ~= 1
        error('statgen:ld', ...
            'single-shard LD loads require a single-shard reference; reference has %d shards', ...
            numel(ref_shards));
    end
    if ~strcmp(ref_shards{1}.label, shard.chr)
        error('statgen:ld', ...
            'single-shard LD chromosome does not match single-shard reference: LD chr %s, reference chr %s', ...
            shard.chr, ref_shards{1}.label);
    end
    statgen.internal.ld_validate_reference_compatibility(shard, ref_shards{1}, path);
    groups = {{shard}};
end
