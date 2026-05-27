function panel = load_annotation(path, reference, varargin)
%LOAD_ANNOTATION Load one BED-like annotation source aligned to a ReferencePanel.
%
%   annotations = statgen.load_annotation(path, reference)
%   annotations = statgen.load_annotation(..., 'has_header', true, 'value_columns', columns)
%   annotations = statgen.load_annotation(..., 'group_column', column)
%   annotations = statgen.load_annotation(..., 'annotation_names', names)
%   annotations = statgen.load_annotation(..., 'annotation_metadata', metadata)
%   annotations = statgen.load_annotation(..., 'annotation_metadata_path', path)
%
% Three-column input with value_columns omitted produces one binary annotation
% named from the file stem. Four-column input with value_columns omitted selects
% column 4 automatically. Five or more columns require explicit value_columns
% or group_column. Inputs with numeric value columns produce continuous
% annotation columns. group_column produces one binary annotation per distinct
% group value. value_columns and group_column use MATLAB one-based physical
% column indices or header names.
% annotation_names overrides names inferred from the header or file stem.
% annotation_metadata is a cell or string vector with one entry per output
% column. Provide annotation_metadata or annotation_metadata_path to attach
% metadata strings; omitting both generates stable JSON provenance strings.
%
% See also statgen.load_annotations, statgen.AnnotationPanel.
    if nargin < 2
        error('statgen:arg', 'load_annotation requires path and reference');
    end
    opts = parse_options_(varargin{:});
    path = char(path);
    if exist(path, 'file') ~= 2
        error('statgen:io', 'annotation file not found: %s', path);
    end

    validate_group_options_(opts);
    [cols, header_fields, row_base0, n_cols, value_columns, value_column_metadata, group_column, group_column_metadata] = read_annotation_table_(path, opts.has_header, opts.value_columns, opts.group_column);
    if n_cols < 3
        error('statgen:annotations', '%s: annotation input must have at least 3 tab-separated columns', path);
    end

    [chr_col, starts, ends] = validate_interval_columns_(cols, path, row_base0);

    num_source_intervals = numel(cols{1});

    if ~isempty(group_column)
        group_values = statgen.internal.ensure_cell_col(cols{group_column});
        empty_group = cellfun('isempty', group_values);
        if any(empty_group)
            row = find(empty_group, 1, 'first');
            error('statgen:annotations', '%s: row %d: group_column values must be non-empty', path, row_base0 + row);
        end
        names = unique(group_values, 'stable');
        annomat = sparse(reference.num_snp, numel(names));
        metadata = cell(numel(names), 1);
        % PERF: loop retained because grouped binary interval union is defined independently per group.
        for i = 1:numel(names)
            rows = strcmp(group_values, names{i});
            annomat(:, i) = paint_binary_(chr_col(rows), starts(rows), ends(rows), reference);
            metadata{i} = statgen.internal.annotation_generated_metadata( ...
                path, opts.has_header, sum(rows), ...
                'group_column', group_column_metadata{1}, ...
                'group_value', names{i});
        end
        is_binary = true(numel(names), 1);
    elseif isempty(value_columns) && n_cols == 3
        names = default_or_explicit_names_(path, n_cols, header_fields, value_columns, opts.annotation_names);
        if numel(names) ~= 1
            error('statgen:annotations', 'annotation_names length mismatch: expected 1');
        end
        annomat = paint_binary_(chr_col, starts, ends, reference);
        is_binary = true;
        if ~isempty(opts.annotation_metadata_path)
            metadata = {statgen.internal.annotation_read_sidecar_exact(opts.annotation_metadata_path)};
        elseif ~isempty(opts.annotation_metadata)
            metadata = statgen.internal.annotation_coerce_metadata_vector(opts.annotation_metadata, 1, 'annotation_metadata');
        else
            metadata = {statgen.internal.annotation_generated_metadata(path, opts.has_header, num_source_intervals)};
        end
    else
        if isempty(value_columns)
            if n_cols == 4
                value_columns = 4;
                if isempty(header_fields)
                    value_column_metadata = {4};
                else
                    value_column_metadata = {header_fields{4}};
                end
            else
                error('statgen:annotations', '%s: input with five or more columns requires explicit value_columns or group_column', path);
            end
        end
        names = default_or_explicit_names_(path, n_cols, header_fields, value_columns, opts.annotation_names);
        if numel(names) ~= numel(value_columns)
            error('statgen:annotations', 'annotation_names length mismatch: expected %d, got %d', numel(value_columns), numel(names));
        end

        statgen.internal.annotation_validate_numeric_nonoverlap(chr_col, starts, ends, path, row_base0);
        values = numeric_value_matrix_(cols, value_columns, path, row_base0);
        annomat = paint_numeric_(chr_col, starts, ends, values, reference);
        is_binary = false(numel(value_columns), 1);

        if ~isempty(opts.annotation_metadata_path)
            sidecar_lines = read_column_metadata_sidecar_(opts.annotation_metadata_path, n_cols);
            metadata = sidecar_lines(value_columns);
        elseif ~isempty(opts.annotation_metadata)
            metadata = statgen.internal.annotation_coerce_metadata_vector(opts.annotation_metadata, numel(value_columns), 'annotation_metadata');
        else
            metadata = cell(numel(value_columns), 1);
            for i = 1:numel(value_columns)
                metadata{i} = statgen.internal.annotation_generated_metadata( ...
                    path, opts.has_header, num_source_intervals, ...
                    'value_column', value_column_metadata{i});
            end
        end
    end

    panel = statgen.create_annotations(reference, annomat, names, is_binary, metadata);
end

function opts = parse_options_(varargin)
    opts.has_header = false;
    opts.value_columns = [];
    opts.group_column = [];
    opts.annotation_names = [];
    opts.annotation_metadata = [];
    opts.annotation_metadata_path = [];
    if mod(numel(varargin), 2) ~= 0
        error('statgen:annotations', 'load_annotation options must be name-value pairs');
    end
    for i = 1:2:numel(varargin)
        if ~(ischar(varargin{i}) || isstring(varargin{i}))
            error('statgen:annotations', 'load_annotation option names must be strings');
        end
        name = char(varargin{i});
        switch name
            case 'has_header'
                opts.has_header = logical_scalar_(varargin{i + 1}, 'has_header');
            case 'value_columns'
                opts.value_columns = varargin{i + 1};
            case 'group_column'
                opts.group_column = varargin{i + 1};
            case 'annotation_names'
                opts.annotation_names = varargin{i + 1};
            case 'annotation_metadata'
                opts.annotation_metadata = varargin{i + 1};
            case 'annotation_metadata_path'
                opts.annotation_metadata_path = char(varargin{i + 1});
            otherwise
                error('statgen:annotations', 'Unknown load_annotation option: %s', name);
        end
    end
    if ~isempty(opts.annotation_metadata) && ~isempty(opts.annotation_metadata_path)
        error('statgen:annotations', 'load_annotation accepts at most one of annotation_metadata and annotation_metadata_path');
    end
end

function validate_group_options_(opts)
    if isempty(opts.group_column)
        return
    end
    if ~isempty(opts.value_columns)
        error('statgen:annotations', 'group_column and value_columns are mutually exclusive');
    end
    if ~isempty(opts.annotation_names)
        error('statgen:annotations', 'annotation_names is invalid with group_column');
    end
    if ~isempty(opts.annotation_metadata)
        error('statgen:annotations', 'annotation_metadata is invalid with group_column');
    end
    if ~isempty(opts.annotation_metadata_path)
        error('statgen:annotations', 'annotation_metadata_path is invalid with group_column');
    end
end

function out = logical_scalar_(value, name)
    if ~isscalar(value) || (ischar(value) || isstring(value))
        error('statgen:annotations', '%s must be a scalar logical', name);
    end
    out = logical(value);
end

function [cols, header_fields, row_base0, n_cols, value_columns, value_column_metadata, group_column, group_column_metadata] = read_annotation_table_(path, has_header, value_columns_raw, group_column_raw)
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open annotation file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));

    first_line = '';
    first_pos = -1;
    while true
        pos = ftell(fid);
        line = fgetl(fid);
        if ~ischar(line)
            error('statgen:annotations', '%s: BED file is empty', path);
        end
        if isempty(line) || strncmp(line, '#', 1)
            continue
        end
        first_line = line;
        first_pos = pos;
        break
    end

    n_cols = numel(strfind(first_line, sprintf('\t'))) + 1;
    if n_cols < 3
        error('statgen:annotations', '%s: BED must have at least 3 tab-separated columns', path);
    end

    if has_header
        header_fields = strsplit(first_line, '\t');
        if any(cellfun('isempty', header_fields))
            error('statgen:annotations', '%s: header names must be non-empty', path);
        end
        if numel(unique(header_fields)) ~= numel(header_fields)
            error('statgen:annotations', '%s: header names must be unique', path);
        end
        row_base0 = 1;
    else
        header_fields = [];
        row_base0 = 0;
        fseek(fid, first_pos, 'bof');
    end

    [value_columns, value_column_metadata] = normalize_column_selectors_(value_columns_raw, header_fields, n_cols, path, 'value_columns');
    [group_column, group_column_metadata] = normalize_column_selectors_(group_column_raw, header_fields, n_cols, path, 'group_column');
    if ~isempty(value_columns) && ~isempty(group_column)
        error('statgen:annotations', 'group_column and value_columns are mutually exclusive');
    end
    if numel(group_column) > 1
        error('statgen:annotations', 'group_column must identify exactly one source column');
    end
    if isempty(value_columns) && isempty(group_column) && n_cols == 4
        value_columns = 4;
        if isempty(header_fields)
            value_column_metadata = {4};
        else
            value_column_metadata = {header_fields{4}};
        end
    end
    fmt_parts = repmat({'%s'}, 1, n_cols);
    for i = 1:numel(value_columns)
        fmt_parts{value_columns(i)} = '%f';
    end
    fmt = strjoin(fmt_parts, '');
    cols = textscan(fid, fmt, ...
        'Delimiter', '\t', ...
        'Whitespace', '', ...
        'MultipleDelimsAsOne', false, ...
        'EmptyValue', NaN, ...
        'TreatAsEmpty', {'NaN', 'nan', 'Inf', 'inf', '-Inf', '-inf'}, ...
        'ReturnOnError', false);
    if isempty(cols) || numel(cols{1}) == 0
        error('statgen:annotations', '%s: BED file is empty', path);
    end
    n_rows = numel(cols{1});
    for i = 2:numel(cols)
        if numel(cols{i}) ~= n_rows
            error('statgen:annotations', '%s: malformed annotation input', path);
        end
    end
end

function [columns, metadata_values] = normalize_column_selectors_(selectors_raw, header_fields, n_cols, path, name)
    if isempty(selectors_raw)
        columns = [];
        metadata_values = {};
        return
    end
    if ischar(selectors_raw)
        selectors = {selectors_raw};
    elseif isstring(selectors_raw)
        selectors = cellstr(selectors_raw(:));
    elseif iscell(selectors_raw)
        selectors = selectors_raw(:);
    else
        selectors = num2cell(selectors_raw(:));
    end
    if isempty(selectors)
        error('statgen:annotations', '%s must be a non-empty list', name);
    end
    columns = zeros(numel(selectors), 1);
    metadata_values = cell(numel(selectors), 1);
    for i = 1:numel(selectors)
        selector = selectors{i};
        if ischar(selector) || isstring(selector)
            if isempty(header_fields)
                error('statgen:annotations', 'named %s are invalid when has_header=false', name);
            end
            idx = find(strcmp(header_fields, char(selector)), 1, 'first');
            if isempty(idx)
                error('statgen:annotations', '%s: unknown %s name %s', path, name, char(selector));
            end
            col = idx;
            metadata_value = char(selector);
        elseif isnumeric(selector) && isscalar(selector) && selector == floor(selector)
            col = double(selector);
            metadata_value = col;
        else
            error('statgen:annotations', '%s must contain column names or integer indices', name);
        end
        if col < 4 || col > n_cols
            error('statgen:annotations', ...
                '%s: %s is out of range; selected columns must be physical columns 4 or later', path, name);
        end
        columns(i) = col;
        metadata_values{i} = metadata_value;
    end
    if numel(unique(columns)) ~= numel(columns)
        error('statgen:annotations', '%s must not contain duplicates', name);
    end
end

function [chr_col, starts, ends] = validate_interval_columns_(cols, path, row_base0)
    chr_col = statgen.internal.ensure_cell_col(cols{1});
    statgen.internal.validate_variant_chr_labels(chr_col, path, row_base0, 'statgen:annotations');

    start_raw = statgen.internal.ensure_cell_col(cols{2});
    end_raw = statgen.internal.ensure_cell_col(cols{3});
    starts = str2double(start_raw);
    bad_start = isnan(starts) | (starts ~= floor(starts)) | (starts < 0);
    if any(bad_start)
        i = find(bad_start, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED start must be a non-negative integer', path, row_base0 + i);
    end
    ends = str2double(end_raw);
    bad_end = isnan(ends) | (ends ~= floor(ends)) | (ends < 0);
    if any(bad_end)
        i = find(bad_end, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED end must be a non-negative integer', path, row_base0 + i);
    end
    bad_len = ends < starts;
    if any(bad_len)
        i = find(bad_len, 1, 'first');
        error('statgen:annotations', '%s: row %d: BED interval end must be >= start', path, row_base0 + i);
    end
end

function A = paint_binary_(chr_col, starts, ends, reference)
    n = reference.num_snp;
    A = sparse(n, 1);
    labels = unique(chr_col, 'stable');
    for i = 1:numel(labels)
        label = labels{i};
        if strcmp(label, 'Y') || strcmp(label, 'MT')
            continue
        end
        rows = strcmp(chr_col, label);
        intervals = statgen.internal.annotation_merge_intervals(starts(rows), ends(rows));
        for s = 1:numel(reference.shards)
            ref_shard = reference.shards{s};
            if ~strcmp(ref_shard.label, label)
                continue
            end
            off = reference.shard_offsets(s);
            ix = (off.start0 + 1):off.stop0;
            mask = statgen.internal.annotation_paint_mask(ref_shard.bp, intervals);
            nz = find(mask ~= 0);
            if ~isempty(nz)
                A(ix(nz), 1) = 1;
            end
        end
    end
end

function values = numeric_value_matrix_(cols, value_columns, path, row_base0)
    values = zeros(numel(cols{1}), numel(value_columns));
    for i = 1:numel(value_columns)
        col = cols{value_columns(i)};
        if ~isnumeric(col)
            col = str2double(statgen.internal.ensure_cell_col(col));
        end
        col = double(col(:));
        bad = ~isfinite(col);
        if any(bad)
            row = find(bad, 1, 'first');
            error('statgen:annotations', ...
                '%s: row %d: annotation value column %d must be finite numeric', ...
                path, row_base0 + row, value_columns(i) - 1);
        end
        values(:, i) = col;
    end
end

function A = paint_numeric_(chr_col, starts, ends, values, reference)
    k = size(values, 2);
    blocks = cell(numel(reference.shards), 1);
    for s = 1:numel(reference.shards)
        blocks{s} = sparse(reference.shards{s}.num_snp, k);
    end
    labels = unique(chr_col, 'stable');
    for i = 1:numel(labels)
        label = labels{i};
        if strcmp(label, 'Y') || strcmp(label, 'MT')
            continue
        end
        rows = find(strcmp(chr_col, label));
        [~, order] = sortrows([double(starts(rows)), double(ends(rows))], [1, 2]);
        sorted_rows = rows(order);
        intervals = [starts(sorted_rows), ends(sorted_rows)];
        value_block = values(sorted_rows, :);
        for s = 1:numel(reference.shards)
            ref_shard = reference.shards{s};
            if ~strcmp(ref_shard.label, label)
                continue
            end
            blocks{s} = statgen.internal.annotation_paint_numeric_sparse(ref_shard.bp, intervals, value_block);
        end
    end
    A = vertcat(blocks{:});
end

function names = default_or_explicit_names_(path, n_cols, header_fields, value_columns, explicit)
    if ~isempty(explicit)
        names = ensure_names_(explicit);
        return
    end
    if isempty(value_columns)
        if n_cols == 3
            names = {path_stem_(path)};
        elseif n_cols == 4
            if isempty(header_fields)
                names = {path_stem_(path)};
            else
                names = {header_fields{4}};
            end
        else
            error('statgen:annotations', '%s: input with five or more columns requires explicit value_columns or group_column', path);
        end
    else
        if isempty(header_fields) && n_cols >= 5
            error('statgen:annotations', '%s: headerless input with five or more columns requires annotation_names', path);
        elseif isempty(header_fields)
            names = {path_stem_(path)};
        else
            names = header_fields(value_columns);
            names = names(:);
        end
    end
    names = ensure_names_(names);
end

function out = source_names_(header_fields, value_columns)
    out = cell(numel(value_columns), 1);
    if isempty(header_fields)
        return
    end
    for i = 1:numel(value_columns)
        out{i} = header_fields{value_columns(i)};
    end
end

function out = ensure_names_(names)
    if ischar(names) || isstring(names)
        names = cellstr(names(:));
    end
    if ~iscell(names)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end
    out = statgen.internal.ensure_cell_col(names);
    if isempty(out)
        error('statgen:annotations', 'annotation names must be a non-empty list of unique strings');
    end
    if any(cellfun('isempty', out))
        error('statgen:annotations', 'annotation names must not contain empty strings');
    end
    if numel(unique(out)) ~= numel(out)
        error('statgen:annotations', 'annotation names must be unique');
    end
end

function out = path_stem_(path)
    [~, stem, ~] = fileparts(path);
    if isempty(stem)
        error('statgen:annotations', 'invalid annotation filename: %s', path);
    end
    out = stem;
end

function out = read_column_metadata_sidecar_(path, n_cols)
    text = statgen.internal.annotation_read_sidecar_exact(path);
    if isempty(text)
        lines = {};
    else
        lines = regexp(text, '\r\n|\n|\r', 'split');
        if ~isempty(lines) && isempty(lines{end}) && ~isempty(regexp(text, '\r\n$|\n$|\r$', 'once'))
            lines(end) = [];
        end
    end
    if numel(lines) ~= n_cols
        error('statgen:annotations', ...
            '%s: annotation metadata sidecar line count mismatch: expected %d, got %d', ...
            path, n_cols, numel(lines));
    end
    if any(cellfun('isempty', lines))
        idx = find(cellfun('isempty', lines), 1, 'first');
        error('statgen:annotations', '%s: annotation metadata sidecar line %d is empty', path, idx);
    end
    out = lines(:);
end
