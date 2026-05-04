function manifest = convert_ld_npz_to_mat(input_root, output_root, production)
% Convert Python LD .npz distribution artifacts into MATLAB/Octave .mat files.
%
% MATLAB production conversions write MAT-file v7.3 output by default.
% Octave writes v5 sparse .mat files only when production=false is passed
% explicitly for fixture-scale tests and local validation; Octave output is
% not a production LD distribution artifact.
    if nargin < 3 || isempty(production)
        production = true;
    end
    input_root = char(input_root);
    output_root = char(output_root);

    if production && is_octave_()
        error('statgen:ld', ...
            'Octave cannot write production LD .mat distributions; use MATLAB for v7.3 output');
    end
    if exist(output_root, 'dir') ~= 7
        mkdir(output_root);
    end

    input_manifest = statgen.internal.ld_read_manifest( ...
        fullfile(input_root, 'ld_manifest.json'), 'python_npz_csc32');

    records = struct('chr', {}, 'sex', {}, 'file', {}, 'file_md5', {}, ...
        'num_snp', {}, 'nnz', {}, 'reference_checksum', {}, 'reference_bim', {});
    copied_reference_bim = {};
    for i = 1:numel(input_manifest.shards)
        entry = input_manifest.shards(i);
        if ~any(strcmp(copied_reference_bim, entry.reference_bim))
            copy_reference_bim_(input_root, output_root, entry.reference_bim);
            copied_reference_bim{end+1} = char(entry.reference_bim); %#ok<AGROW>
        end
        npz_path = fullfile(input_root, entry.file);
        [ld_r, a1freq, metadata] = read_npz_ld_shard_(npz_path, output_root);
        validate_npz_entry_agreement_(entry, metadata, npz_path);

        metadata = convert_metadata_(metadata);
        mat_file = strrep(entry.file, '.npz', '.mat');
        mat_path = fullfile(output_root, mat_file);
        save_mat_shard_(mat_path, ld_r, a1freq, metadata, production);

        records(end+1).chr = char(metadata.chr); %#ok<AGROW>
        records(end).sex = metadata.sex;
        records(end).file = mat_file;
        records(end).file_md5 = statgen.internal.ld_md5_file(mat_path);
        records(end).num_snp = double(metadata.num_snp);
        records(end).nnz = double(metadata.nnz);
        records(end).reference_checksum = char(metadata.reference_checksum);
        records(end).reference_bim = char(entry.reference_bim);
    end

    manifest = struct();
    manifest.object_type = 'ld_panel_manifest';
    manifest.schema_version = '1.0';
    manifest.runtime_format = 'matlab_mat_sparse_double';
    manifest.shards = records;
    write_json_(fullfile(output_root, 'ld_manifest.json'), manifest);

    statgen.validate_ld_distribution(output_root, false);
end

function copy_reference_bim_(input_root, output_root, reference_bim)
    src = fullfile(input_root, reference_bim);
    dst = fullfile(output_root, reference_bim);
    if exist(src, 'file') ~= 2
        error('statgen:io', 'LD file not found: %s', src);
    end
    [ok, msg] = copyfile(src, dst);
    if ~ok
        error('statgen:io', 'Failed to copy bundled reference_bim: %s', msg);
    end
end

function [ld_r, a1freq, metadata] = read_npz_ld_shard_(path, scratch_root)
    if exist(path, 'file') ~= 2
        error('statgen:io', 'LD file not found: %s', path);
    end
    [~, tag] = fileparts(path);
    work = fullfile(scratch_root, ['.tmp_' tag]);
    mkdir(work);
    cleanup = onCleanup(@() cleanup_dir_(work));
    unzip(path, work);

    data = read_npy_(fullfile(work, 'data.npy'));
    indices = read_npy_(fullfile(work, 'indices.npy'));
    indptr = read_npy_(fullfile(work, 'indptr.npy'));
    shape = read_npy_(fullfile(work, 'shape.npy'));
    a1freq = read_npy_(fullfile(work, 'a1freq.npy'));
    metadata_bytes = read_npy_(fullfile(work, 'metadata.npy'));

    metadata = jsondecode(char(uint8(metadata_bytes(:))'));
    if ~strcmp(metadata.format, 'statgen_ld_npz_csc32')
        error('statgen:ld', '%s: metadata format must be statgen_ld_npz_csc32', path);
    end
    if ~strcmp(metadata.sparse_layout, 'csc') || double(metadata.index_base) ~= 0
        error('statgen:ld', '%s: metadata sparse layout must be zero-based CSC', path);
    end

    num_snp = double(metadata.num_snp);
    if ~isequal(double(shape(:))', [num_snp, num_snp])
        error('statgen:ld', '%s: shape must equal metadata num_snp', path);
    end
    if numel(indptr) ~= num_snp + 1
        error('statgen:ld', '%s: indptr length must be num_snp + 1', path);
    end
    if numel(data) ~= numel(indices) || numel(data) ~= double(metadata.nnz)
        error('statgen:ld', '%s: data, indices, and metadata nnz length mismatch', path);
    end
    if double(indptr(end)) ~= double(metadata.nnz)
        error('statgen:ld', '%s: indptr(end) must equal metadata nnz', path);
    end
    if numel(a1freq) ~= num_snp
        error('statgen:ld', '%s: a1freq length must equal num_snp', path);
    end

    rows = zeros(numel(data), 1);
    cols = zeros(numel(data), 1);
    for j = 1:num_snp
        start0 = double(indptr(j));
        stop0 = double(indptr(j + 1));
        if stop0 < start0
            error('statgen:ld', '%s: indptr must be monotonic nondecreasing', path);
        end
        if stop0 > start0
            idx = (start0 + 1):stop0;
            rows(idx) = double(indices(idx)) + 1;
            cols(idx) = j;
        end
    end
    if any(rows < 1 | rows > num_snp)
        error('statgen:ld', '%s: sparse row indices out of bounds', path);
    end
    ld_r = sparse(rows, cols, double(data(:)), num_snp, num_snp);
    a1freq = double(a1freq(:));
end

function metadata = convert_metadata_(metadata)
    metadata.format = 'statgen_ld_mat_sparse_double';
    if isfield(metadata, 'sparse_layout')
        metadata = rmfield(metadata, 'sparse_layout');
    end
    if isfield(metadata, 'index_base')
        metadata = rmfield(metadata, 'index_base');
    end
end

function validate_npz_entry_agreement_(entry, metadata, path)
    if ~strcmp(entry.chr, metadata.chr)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for chr', path);
    end
    if ~statgen.internal.ld_same_optional_string(entry.sex, metadata.sex)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for sex', path);
    end
    if double(entry.num_snp) ~= double(metadata.num_snp)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for num_snp', path);
    end
    if double(entry.nnz) ~= double(metadata.nnz)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for nnz', path);
    end
    if ~strcmp(entry.reference_checksum, metadata.reference_checksum)
        error('statgen:ld', '%s: manifest/per-file metadata mismatch for reference_checksum', path);
    end
end

function save_mat_shard_(path, ld_r, a1freq, metadata, production)
    if production
        save(path, 'ld_r', 'a1freq', 'metadata', '-v7.3');
    else
        save(path, 'ld_r', 'a1freq', 'metadata', '-mat');
    end
end

function value = read_npy_(path)
    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Missing NPY member: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));

    magic = fread(fid, 6, '*uint8')';
    if numel(magic) ~= 6 || ~isequal(magic, uint8([147, double('NUMPY')]))
        error('statgen:ld', '%s: invalid NPY magic', path);
    end
    version = fread(fid, 2, '*uint8')';
    if version(1) == 1
        header_len = double(fread(fid, 1, 'uint16', 0, 'ieee-le'));
    elseif version(1) == 2
        header_len = double(fread(fid, 1, 'uint32', 0, 'ieee-le'));
    else
        error('statgen:ld', '%s: unsupported NPY version', path);
    end
    header = char(fread(fid, header_len, '*char')');
    descr = regexp(header, '''descr'':\s*''([^'']+)''', 'tokens', 'once');
    fortran = regexp(header, '''fortran_order'':\s*(True|False)', 'tokens', 'once');
    shape_text = regexp(header, '''shape'':\s*\(([^\)]*)\)', 'tokens', 'once');
    if isempty(descr) || isempty(fortran) || isempty(shape_text)
        error('statgen:ld', '%s: unsupported NPY header', path);
    end
    if strcmp(fortran{1}, 'True')
        error('statgen:ld', '%s: Fortran-order NPY arrays are not supported', path);
    end
    dims = parse_shape_(shape_text{1});
    raw = fread(fid, Inf, '*uint8');
    value = decode_npy_data_(raw, descr{1});
    if prod(dims) ~= numel(value)
        error('statgen:ld', '%s: NPY data length does not match shape', path);
    end
    value = reshape(value, dims);
end

function dims = parse_shape_(text)
    parts = strsplit(strtrim(text), ',');
    vals = [];
    for i = 1:numel(parts)
        p = strtrim(parts{i});
        if isempty(p)
            continue;
        end
        vals(end+1) = str2double(p); %#ok<AGROW>
    end
    if isempty(vals)
        dims = [1, 1];
    elseif numel(vals) == 1
        dims = [vals(1), 1];
    else
        dims = vals;
    end
end

function value = decode_npy_data_(raw, descr)
    switch descr
        case '<f4'
            value = typecast(uint8(raw), 'single')';
        case '<i4'
            value = typecast(uint8(raw), 'int32')';
        case '<i8'
            value = typecast(uint8(raw), 'int64')';
        case '|u1'
            value = uint8(raw);
        otherwise
            error('statgen:ld', 'Unsupported NPY dtype: %s', descr);
    end
end

function write_json_(path, obj)
    fid = fopen(path, 'w');
    if fid < 0
        error('statgen:io', 'Cannot write manifest: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    fprintf(fid, '%s\n', jsonencode(obj));
    clear cleaner;
end

function cleanup_dir_(path)
    if exist(path, 'dir') == 7
        rmdir(path, 's');
    end
end

function tf = is_octave_()
    persistent cached;
    if isempty(cached)
        cached = exist('OCTAVE_VERSION', 'builtin') ~= 0;
    end
    tf = cached;
end
