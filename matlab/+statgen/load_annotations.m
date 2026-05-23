function panel = load_annotations(bed_paths, reference, varargin)
%LOAD_ANNOTATIONS Load BED annotations aligned to a ReferencePanel.
%
%   annotations = statgen.load_annotations(bed_paths, reference)
%   annotations = statgen.load_annotations(..., 'annotation_metadata', metadata)
%   annotations = statgen.load_annotations(..., 'annotation_metadata_paths', paths)
%
% Paints one or more BED files onto the SNPs in reference. Each BED file is one
% annotation column, not a chromosome shard. Annotation names are derived from
% BED file basenames. The returned AnnotationPanel has a sparse numeric
% num_snp-by-num_annot annomat property with is_binary true for each column.
%
% BED inputs use columns 1-3 only: chromosome, 0-based start, and 0-based
% exclusive end. Chromosome labels must match the reference exactly; statgen
% does not strip 'chr' prefixes or otherwise normalize labels.
% BED basenames become annotation names and must be unique.
%
% See also statgen.AnnotationPanel, statgen.load_annotation, statgen.create_annotations,
% statgen.load_annotations_cache.
    if nargin < 2
        error('statgen:arg', 'load_annotations requires bed_paths and reference');
    end

    opts = parse_options_(varargin{:});
    paths = coerce_bed_paths_(bed_paths);
    annonames = cell(numel(paths), 1);
    for i = 1:numel(paths)
        if exist(paths{i}, 'file') ~= 2
            error('statgen:io', 'BED file not found: %s', paths{i});
        end
        [~, stem, ~] = fileparts(paths{i});
        if isempty(stem)
            error('statgen:annotations', 'invalid annotation filename: %s', paths{i});
        end
        annonames{i} = stem;
    end
    if numel(unique(annonames)) ~= numel(annonames)
        error('statgen:annotations', 'duplicate annotation names derived from BED basenames');
    end

    n = reference.num_snp;
    k = numel(paths);
    if n == 0
        panel = statgen.create_annotations(reference, sparse([], [], [], 0, k), annonames, true(k, 1), metadata_for_paths_(paths, opts));
        return
    end

    A = sparse(n, k);
    for j = 1:k
        intervals_map = parse_bed_intervals_(paths{j});
        for s = 1:numel(reference.shards)
            ref_shard = reference.shards{s};
            off = reference.shard_offsets(s);
            ix = (off.start0 + 1):off.stop0;
            mask = zeros(numel(ix), 1);
            idx = find(strcmp(intervals_map.labels, ref_shard.label), 1, 'first');
            if ~isempty(idx)
                mask = statgen.internal.annotation_paint_mask(ref_shard.bp, intervals_map.intervals{idx});
            end
            nz = find(mask ~= 0);
            if ~isempty(nz)
                A(ix(nz), j) = 1;
            end
        end
    end

    panel = statgen.create_annotations(reference, A, annonames, true(k, 1), metadata_for_paths_(paths, opts));
end

function opts = parse_options_(varargin)
    opts.annotation_metadata = [];
    opts.annotation_metadata_paths = [];
    if mod(numel(varargin), 2) ~= 0
        error('statgen:annotations', 'load_annotations options must be name-value pairs');
    end
    for i = 1:2:numel(varargin)
        if ~(ischar(varargin{i}) || isstring(varargin{i}))
            error('statgen:annotations', 'load_annotations option names must be strings');
        end
        name = char(varargin{i});
        switch name
            case 'annotation_metadata'
                opts.annotation_metadata = varargin{i + 1};
            case 'annotation_metadata_paths'
                opts.annotation_metadata_paths = varargin{i + 1};
            otherwise
                error('statgen:annotations', 'Unknown load_annotations option: %s', name);
        end
    end
    if ~isempty(opts.annotation_metadata) && ~isempty(opts.annotation_metadata_paths)
        error('statgen:annotations', 'load_annotations accepts at most one of annotation_metadata and annotation_metadata_paths');
    end
end

function metadata = metadata_for_paths_(paths, opts)
    n = numel(paths);
    if ~isempty(opts.annotation_metadata)
        metadata = statgen.internal.annotation_coerce_metadata_vector(opts.annotation_metadata, n, 'annotation_metadata');
    elseif ~isempty(opts.annotation_metadata_paths)
        meta_paths = coerce_metadata_paths_(opts.annotation_metadata_paths, n);
        metadata = cell(n, 1);
        for i = 1:n
            if isempty(meta_paths{i})
                metadata{i} = statgen.internal.annotation_generated_metadata(paths{i}, [], []);
            else
                metadata{i} = statgen.internal.annotation_read_sidecar_exact(meta_paths{i});
            end
        end
    else
        metadata = cell(n, 1);
        for i = 1:n
            metadata{i} = statgen.internal.annotation_generated_metadata(paths{i}, [], []);
        end
    end
end

function out = coerce_metadata_paths_(value, expected_len)
    if ischar(value) || isstring(value)
        value = cellstr(value(:));
    end
    if ~iscell(value)
        error('statgen:annotations', 'annotation_metadata_paths must have one entry per BED path');
    end
    out = statgen.internal.ensure_cell_col(value);
    if numel(out) ~= expected_len
        error('statgen:annotations', 'annotation_metadata_paths length mismatch: expected %d, got %d', expected_len, numel(out));
    end
    for i = 1:numel(out)
        out{i} = char(out{i});
    end
end

function paths = coerce_bed_paths_(bed_paths)
    if ischar(bed_paths)
        bed_paths = {bed_paths};
    elseif isstring(bed_paths)
        bed_paths = cellstr(bed_paths(:));
    end
    if ~iscell(bed_paths)
        error('statgen:annotations', 'bed_paths must be a non-empty list of BED files');
    end

    paths = statgen.internal.ensure_cell_col(bed_paths);
    if isempty(paths)
        error('statgen:annotations', 'bed_paths must be a non-empty list of BED files');
    end
    for i = 1:numel(paths)
        paths{i} = char(paths{i});
    end
end

function intervals_map = parse_bed_intervals_(path)
    [cols, n_rows] = read_bed_tabular_(path);
    if n_rows == 0
        error('statgen:annotations', '%s: BED file is empty', path);
    end
    chr_col = statgen.internal.ensure_cell_col(cols{1});
    start_raw = statgen.internal.ensure_cell_col(cols{2});
    end_raw = statgen.internal.ensure_cell_col(cols{3});

    statgen.internal.validate_variant_chr_labels(chr_col, path, 0, 'statgen:annotations');

    start_num = str2double(start_raw);
    bad_start = isnan(start_num) | (start_num ~= floor(start_num)) | (start_num < 0);
    if any(bad_start)
        i = find(bad_start, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED start must be a non-negative integer', path, i);
    end

    end_num = str2double(end_raw);
    bad_end = isnan(end_num) | (end_num ~= floor(end_num)) | (end_num < 0);
    if any(bad_end)
        i = find(bad_end, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED end must be a non-negative integer', path, i);
    end

    bad_len = end_num < start_num;
    if any(bad_len)
        i = find(bad_len, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED interval end must be >= start', path, i);
    end

    intervals_map = struct('labels', {{}}, 'intervals', {{}});

    chr_unique = {};
    for i = 1:numel(chr_col)
        c = chr_col{i};
        if strcmp(c, 'Y') || strcmp(c, 'MT')
            continue
        end
        if ~any(strcmp(chr_unique, c))
            chr_unique{end + 1, 1} = c; %#ok<AGROW>
        end
    end

    for i = 1:numel(chr_unique)
        mask = strcmp(chr_col, chr_unique{i});
        merged = statgen.internal.annotation_merge_intervals(start_num(mask), end_num(mask));
        intervals_map.labels{end + 1, 1} = chr_unique{i}; %#ok<AGROW>
        intervals_map.intervals{end + 1, 1} = merged; %#ok<AGROW>
    end
end

function [cols, n_rows] = read_bed_tabular_(path)
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open BED file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));

    while true
        pos = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line)
            cols = {cell(0, 1), cell(0, 1), cell(0, 1)};
            n_rows = 0;
            return
        end
        if isempty(line) || strncmp(line, '#', 1)
            continue
        end
        first_pos = pos;
        break
    end

    n_cols = numel(strfind(line, sprintf('\t'))) + 1;
    if n_cols < 3
        error('statgen:annotations', '%s: BED must have at least 3 tab-separated columns', path);
    end
    fseek(fid, first_pos, 'bof');

    fmt_parts = repmat({'%s'}, 1, n_cols);
    fmt = strjoin(fmt_parts, '');
    data = textscan(fid, fmt, ...
        'Delimiter', '\t', ...
        'Whitespace', '', ...
        'MultipleDelimsAsOne', false, ...
        'EmptyValue', NaN, ...
        'TreatAsEmpty', {'NaN', 'nan', 'Inf', 'inf', '-Inf', '-inf'}, ...
        'ReturnOnError', false);
    n_rows = numel(data{1});
    if numel(data{2}) ~= n_rows || numel(data{3}) ~= n_rows
        error('statgen:annotations', '%s: BED must have at least 3 tab-separated columns', path);
    end
    cols = data;
end
