function panel = load_ld(path, varargin)
%LOAD_LD Load a MATLAB sparse LD distribution.
%
%   ld = statgen.load_ld(path)
%   ld = statgen.load_ld(path, reference)
%   ld = statgen.load_ld(path, reference, chrX_sex)
%   ld = statgen.load_ld(path, reference, shards)
%   ld = statgen.load_ld(path, reference, shards, chrX_sex)
%   ld = statgen.load_ld(..., retain_ld_r)
%
% path is an LD distribution directory containing ld_manifest.json and converted
% .mat LD shard files; load_ld does not load a single shard file or @ template.
% If reference is omitted, load_ld uses the reference BIM files bundled with the
% LD distribution. Optional shards restrict loading to selected canonical shard
% labels. chrX_sex selects the default chrX shard; use 'female', 'male', or
% 'combined'.
%
% retain_ld_r defaults to true. Set it to false to drop raw signed LD r from
% loaded shards after constructing the internal r-squared matrix used by
% multiply_r2 and fast_prune.
%
% MATLAB LD distributions are prepared by building Python .npz LD shards with
% script/statgen_build_ld.py, converting them with statgen.convert_ld_npz_to_mat,
% and finalizing the manifest with statgen.create_ld_mat_manifest.
%
% See also statgen.LDPanel, statgen.convert_ld_npz_to_mat,
% statgen.create_ld_mat_manifest, statgen.validate_ld_distribution.
    [reference, shards, default_chrX_sex, retain_ld_r] = parse_args_(varargin{:});
    statgen.LDPanel.validate_chrx_sex_(default_chrX_sex, 'default_chrX_sex');
    if ~islogical(retain_ld_r) && ~(isnumeric(retain_ld_r) && isscalar(retain_ld_r))
        error('statgen:ld', 'retain_ld_r must be a scalar logical value');
    end
    retain_ld_r = logical(retain_ld_r);

    path = char(path);
    if ~isfolder(path)
        error('statgen:ld', 'load_ld: path must identify a panel root directory');
    end

    manifest = statgen.internal.ld_read_manifest(fullfile(path, 'ld_manifest.json'), ...
        'matlab_mat_sparse_double');
    if isempty(reference)
        reference = load_bundled_reference_(path, manifest, shards);
    elseif ~isempty(shards)
        reference = reference.select_shards(shards);
    end
    groups = load_panel_root_(path, manifest, reference.shards, retain_ld_r);
    panel = statgen.LDPanel(groups, default_chrX_sex, reference);
end

function [reference, shards, default_chrX_sex, retain_ld_r] = parse_args_(varargin)
    reference = [];
    shards = [];
    default_chrX_sex = 'female';
    retain_ld_r = true;

    if numel(varargin) >= 1
        reference = varargin{1};
    end
    if numel(varargin) >= 2
        third = varargin{2};
        if is_chrx_sex_(third)
            default_chrX_sex = char(third);
            if numel(varargin) >= 3 && ~isempty(varargin{3})
                retain_ld_r = varargin{3};
            end
        else
            shards = third;
            if numel(varargin) >= 3 && ~isempty(varargin{3})
                default_chrX_sex = varargin{3};
            end
            if numel(varargin) >= 4 && ~isempty(varargin{4})
                retain_ld_r = varargin{4};
            end
        end
    end
    if numel(varargin) > 4
        error('statgen:ld', 'load_ld accepts at most five arguments');
    end
    if isempty(default_chrX_sex)
        default_chrX_sex = 'female';
    end
    if isempty(retain_ld_r)
        retain_ld_r = true;
    end
end

function tf = is_chrx_sex_(value)
    if ~(ischar(value) || (isstring(value) && isscalar(value)))
        tf = false;
        return
    end
    tf = any(strcmp(char(value), {'female', 'male', 'combined'}));
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
            warn_monomorphic_snps_(shard_path, meta);
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

function reference = load_bundled_reference_(root, manifest, shards)
    available = manifest_chr_labels_(manifest);
    selected = statgen.internal.validate_requested_shards(shards, available, 'load_ld');
    ref_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        entries = entries_for_chr_(manifest.shards, label);
        reference_path = fullfile(root, entries(1).reference_bim);
        panel = statgen.load_reference(reference_path);
        if numel(panel.shards) ~= 1
            error('statgen:ld', '%s: bundled reference_bim must contain exactly one shard', reference_path);
        end
        ref = panel.shards{1};
        for j = 1:numel(entries)
            validate_bundled_reference_entry_(ref, entries(j), reference_path);
        end
        ref_shards{i} = ref;
    end
    reference = statgen.ReferencePanel(ref_shards);
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

function entries = entries_for_chr_(entries_in, label)
    keep = false(numel(entries_in), 1);
    for i = 1:numel(entries_in)
        keep(i) = strcmp(entries_in(i).chr, label);
    end
    entries = entries_in(keep);
end

function validate_bundled_reference_entry_(ref, entry, path)
    if ~strcmp(ref.label, entry.chr)
        error('statgen:ld', '%s: bundled reference_bim chr does not match manifest', path);
    end
    if ref.num_snp ~= double(entry.num_snp)
        error('statgen:ld', '%s: bundled reference_bim num_snp does not match manifest', path);
    end
    if ~strcmp(ref.checksum, entry.reference_checksum)
        error('statgen:ld', '%s: bundled reference_bim reference_checksum does not match manifest', path);
    end
end

function warn_monomorphic_snps_(path, meta)
    if isfield(meta, 'num_monomorphic_snps') && double(meta.num_monomorphic_snps) > 0
        warning('statgen:ld:monomorphic', ...
            '%s: LD shard contains %d monomorphic SNPs; LD involving those SNPs is undefined and represented by omitted off-diagonal entries', ...
            path, double(meta.num_monomorphic_snps));
    end
end
