function sumstats = load_sumstats(path, reference)
% Load sumstats TSV(.gz) and align to reference by shard-local numeric keys.
    if nargin < 2
        error('statgen:arg', 'load_sumstats requires path and reference');
    end
    path = char(path);
    if isempty(path)
        error('statgen:arg', 'path must be non-empty');
    end

    [tbl, cleanup_fn] = parse_sumstats_table_(path);
    cleaner = onCleanup(cleanup_fn);
    sumstats = build_sumstats_(tbl, reference, path);
    clear cleaner;
end

function sumstats = build_sumstats_(tbl, reference, path)
    if is_table_like_(tbl)
        var_names = canonicalize_var_names_(tbl.Properties.VariableNames, path);
    else
        var_names = canonicalize_var_names_(tbl.statgen_var_names__(:)', path);
    end
    required = {'chr', 'bp', 'a1', 'a2', 'z', 'n'};
    for i = 1:numel(required)
        if ~any(strcmp(var_names, required{i}))
            error('statgen:sumstats', '%s: missing required column: %s', path, required{i});
        end
    end

    chr_col = ensure_cellstr_col_(get_col_(tbl, var_names, 'chr'));
    bp_num  = to_numeric_col_(get_col_(tbl, var_names, 'bp'));
    a1_col  = ensure_cellstr_col_(get_col_(tbl, var_names, 'a1'));
    a2_col  = ensure_cellstr_col_(get_col_(tbl, var_names, 'a2'));
    z_num   = to_numeric_col_(get_col_(tbl, var_names, 'z'));
    n_num   = to_numeric_col_(get_col_(tbl, var_names, 'n'));

    bad_chr = cellfun('isempty', chr_col);
    if any(bad_chr), error('statgen:sumstats', '%s: row %d: chr must be non-empty', path, find(bad_chr, 1, 'first') + 1); end
    bad_a1 = cellfun('isempty', a1_col);
    if any(bad_a1), error('statgen:sumstats', '%s: row %d: a1 must be non-empty', path, find(bad_a1, 1, 'first') + 1); end
    bad_a2 = cellfun('isempty', a2_col);
    if any(bad_a2), error('statgen:sumstats', '%s: row %d: a2 must be non-empty', path, find(bad_a2, 1, 'first') + 1); end

    bad_bp = isnan(bp_num) | (bp_num ~= floor(bp_num));
    if any(bad_bp)
        i = find(bad_bp, 1, 'first');
        error('statgen:sumstats', '%s: row %d: bp is not an integer', path, i + 1);
    end

    bad_z = ~isfinite(z_num);
    if any(bad_z)
        i = find(bad_z, 1, 'first');
        error('statgen:sumstats', '%s: row %d: z must be finite numeric', path, i + 1);
    end

    bad_n = ~isfinite(n_num);
    if any(bad_n)
        i = find(bad_n, 1, 'first');
        error('statgen:sumstats', '%s: row %d: n must be finite numeric', path, i + 1);
    end

    optional_map = struct('p', [], 'beta', [], 'se', [], 'eaf', [], 'info', []);
    optional_names = fieldnames(optional_map);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        if any(strcmp(var_names, nm))
            optional_map.(nm) = to_numeric_col_(get_col_(tbl, var_names, nm));
        end
    end

    n_ref = double(reference.num_snp);
    src_a1_hash64 = statgen.internal.allele_hash64(a1_col);
    src_a2_hash64 = statgen.internal.allele_hash64(a2_col);

    aligned_z = nan(n_ref, 1);
    aligned_n = nan(n_ref, 1);

    if isempty(optional_map.p)
        p_aligned = [];
        aligned_logp = nan(n_ref, 1);
    else
        p_aligned = nan(n_ref, 1);
        aligned_logp = [];
    end

    aligned_optional = struct('beta', [], 'se', [], 'eaf', [], 'info', []);
    for i = 1:numel(optional_names)
        nm = optional_names{i};
        vals = optional_map.(nm);
        if isempty(vals)
            aligned_optional.(nm) = [];
        else
            aligned_optional.(nm) = nan(n_ref, 1);
        end
    end

    for i = 1:numel(reference.shards)
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        src_idx = find(strcmp(chr_col, s_ref.label));
        loc = statgen.internal.match_shard_numeric( ...
            s_ref.bp, s_ref.a1_hash64, s_ref.a2_hash64, ...
            bp_num(src_idx), src_a1_hash64(src_idx), src_a2_hash64(src_idx), ...
            s_ref.label, 'sumstats');
        has_match = loc > 0;
        matched_ref_ix = ix(has_match);
        matched_src_ix = src_idx(loc(has_match));

        aligned_z(matched_ref_ix) = z_num(matched_src_ix);
        aligned_n(matched_ref_ix) = n_num(matched_src_ix);
        if ~isempty(p_aligned)
            p_aligned(matched_ref_ix) = optional_map.p(matched_src_ix);
        end
        for j = 1:numel(optional_names)
            nm = optional_names{j};
            if ~isempty(aligned_optional.(nm))
                vals = optional_map.(nm);
                tmp = aligned_optional.(nm);
                tmp(matched_ref_ix) = vals(matched_src_ix);
                aligned_optional.(nm) = tmp;
            end
        end
    end

    if ~isempty(p_aligned)
        aligned_logp = nan(n_ref, 1);
        finite_mask = isfinite(p_aligned);
        in_range = finite_mask & p_aligned >= 0 & p_aligned <= 1;
        zero_mask = in_range & p_aligned == 0;
        pos_mask = in_range & p_aligned > 0;
        aligned_logp(zero_mask) = inf;
        aligned_logp(pos_mask) = -log10(p_aligned(pos_mask));
    end

    shards = cell(numel(reference.shards), 1);
    for i = 1:numel(reference.shards)
        s_ref = reference.shards{i};
        off = reference.shard_offsets(i);
        ix = (off.start0 + 1):off.stop0;
        beta_vec = pick_optional_(aligned_optional.beta, ix);
        se_vec = pick_optional_(aligned_optional.se, ix);
        eaf_vec = pick_optional_(aligned_optional.eaf, ix);
        info_vec = pick_optional_(aligned_optional.info, ix);

        shards{i} = statgen.SumstatsShard( ...
            s_ref.label, s_ref.checksum, ...
            aligned_z(ix), aligned_n(ix), aligned_logp(ix), ...
            beta_vec, se_vec, eaf_vec, info_vec);
    end

    sumstats = statgen.Sumstats(shards);
end

function var_names = canonicalize_var_names_(raw_names, path)
    var_names = lower(raw_names);
    for i = 1:numel(var_names)
        if strcmp(var_names{i}, 'pos')
            var_names{i} = 'bp';
        elseif strcmp(var_names{i}, 'effectallele')
            var_names{i} = 'a1';
        elseif strcmp(var_names{i}, 'otherallele')
            var_names{i} = 'a2';
        end
    end
    for i = 1:numel(var_names)
        if sum(strcmp(var_names, var_names{i})) > 1
            error('statgen:sumstats', '%s: duplicate columns after column normalization: %s', path, var_names{i});
        end
    end
end

function [tbl, cleanup_fn] = parse_sumstats_table_(path)
    cleanup_fn = @() [];
    actual_path = path;
    if ends_with_(path, '.gz')
        [~, base, ~] = fileparts(path);
        scratch_root = statgen_scratch_root_(path);
        tmpdir = unique_tmpdir_(scratch_root, ['sumstats_gunzip_' base]);
        mkdir(tmpdir);
        gunzip(path, tmpdir);
        actual_path = fullfile(tmpdir, base);
        cleanup_fn = @() cleanup_tmpdir_(tmpdir);
    end

    fid = fopen(actual_path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open sumstats file: %s', path);
    end
    closer = onCleanup(@() fclose(fid));
    hdr = fgetl(fid);
    if ~ischar(hdr)
        error('statgen:sumstats', 'sumstats file is empty: %s', path);
    end
    names = strsplit(hdr, '\t');
    fmt = build_col_formats_(names);
    cols = textscan(fid, fmt, 'Delimiter', '\t', 'Whitespace', '', 'MultipleDelimsAsOne', false, 'ReturnOnError', false);
    clear closer;
    S = struct();
    S.statgen_var_names__ = names;
    for i = 1:numel(names)
        if isnumeric(cols{i})
            S.(names{i}) = double(cols{i}(:));
        else
            S.(names{i}) = ensure_cellstr_col_(cols{i});
        end
    end
    tbl = S;
end

function col = get_col_(tbl, var_names, name)
    idx = find(strcmp(var_names, name), 1, 'first');
    if isempty(idx)
        error('statgen:sumstats', 'Missing required column: %s', name);
    end
    if is_table_like_(tbl)
        col = tbl{:, idx};
    else
        field_name = tbl.statgen_var_names__{idx};
        col = tbl.(field_name);
    end
end

function tf = is_table_like_(x)
    tf = false;
    if isstruct(x) && isfield(x, 'statgen_var_names__')
        return
    end
    if exist('istable', 'file') == 2
        tf = istable(x);
    else
        tf = isa(x, 'table');
    end
end

function out = pick_optional_(vec, ix)
    if isempty(vec)
        out = [];
    else
        out = vec(ix);
    end
end

function root = statgen_scratch_root_(path)
    env_root = getenv('STATGEN_SCRATCH');
    if ~isempty(env_root)
        root = char(env_root);
    else
        [root, ~, ~] = fileparts(path);
        if isempty(root)
            root = pwd;
        end
    end
    if exist(root, 'dir') ~= 7
        [ok, msg] = mkdir(root);
        if ~ok
            error('statgen:io', 'Cannot create STATGEN scratch directory %s: %s', root, msg);
        end
    end
end

function tmpdir = unique_tmpdir_(root, prefix)
    root = char(root);
    for i = 1:100
        suffix = sprintf('%s_%06d_%06d', prefix, round(1e6 * rand()), i);
        tmpdir = fullfile(root, suffix);
        if exist(tmpdir, 'dir') ~= 7 && exist(tmpdir, 'file') ~= 2
            return
        end
    end
    error('statgen:io', 'Could not allocate unique scratch directory under %s', root);
end

function cleanup_tmpdir_(tmpdir)
    if exist(tmpdir, 'dir') == 7
        try
            rmdir(tmpdir, 's');
        catch
        end
    end
end

function tf = ends_with_(s, suffix)
    n = length(suffix);
    tf = length(s) >= n && strcmp(s(end - n + 1:end), suffix);
end

function fmt = build_col_formats_(raw_names)
    numeric_set = {'bp', 'pos', 'z', 'n', 'p', 'beta', 'se', 'eaf', 'info'};
    parts = cell(1, numel(raw_names));
    for i = 1:numel(raw_names)
        if any(strcmp(numeric_set, lower(raw_names{i})))
            parts{i} = '%f';
        else
            parts{i} = '%s';
        end
    end
    fmt = [parts{:}];
end

function out = to_numeric_col_(x)
    if isnumeric(x)
        out = double(x(:));
    else
        out = str2double(ensure_cellstr_col_(x));
    end
end

function out = ensure_cellstr_col_(x)
    if iscell(x)
        out = x(:);
    elseif isstring(x)
        out = cellstr(x(:));
    elseif ischar(x)
        out = strtrim(cellstr(x));
    elseif isnumeric(x) || islogical(x)
        out = arrayfun(@(v) num2str(v, '%.15g'), x(:), 'UniformOutput', false);
    else
        out = strtrim(cellstr(x));
    end
end
