function panel = load_genotype_cache(path, shards)
% Load genotype metadata from MATLAB .mat cache with optional shard subsetting.
    if nargin < 2
        shards = [];
    end
    path = char(path);
    loaded = load(path);
    require_fields_(loaded, {'metadata', 'is_present', 'ploidy_male', 'ploidy_female', ...
        'source_row0', 'subject_present', 'source_subject_row0', 'fid', 'iid', ...
        'father_id', 'mother_id', 'sex', 'is_male', 'is_female'});

    meta = loaded.metadata;
    is_present = logical(loaded.is_present(:));
    ploidy_male = double(loaded.ploidy_male(:));
    ploidy_female = double(loaded.ploidy_female(:));
    source_row0 = double(loaded.source_row0(:));
    num_snp = numel(is_present);
    if numel(ploidy_male) ~= num_snp || numel(ploidy_female) ~= num_snp || numel(source_row0) ~= num_snp
        error('statgen:cache', 'Invalid genotype cache: SNP-axis vector lengths mismatch');
    end

    [labels, checksums, start0, stop0, source_layout, bed_paths, bed_file_sizes, source_num_snp, source_num_sample, num_sample] = validate_metadata_(meta, num_snp);

    fid = statgen.internal.ensure_cell_col(loaded.fid);
    iid = statgen.internal.ensure_cell_col(loaded.iid);
    father_id = statgen.internal.ensure_cell_col(loaded.father_id);
    mother_id = statgen.internal.ensure_cell_col(loaded.mother_id);
    sex = double(loaded.sex(:));
    is_male = logical(loaded.is_male(:));
    is_female = logical(loaded.is_female(:));
    if numel(fid) ~= num_sample || numel(iid) ~= num_sample || numel(father_id) ~= num_sample || numel(mother_id) ~= num_sample || numel(sex) ~= num_sample || numel(is_male) ~= num_sample || numel(is_female) ~= num_sample
        error('statgen:cache', 'Invalid genotype cache: FAM vector length mismatch');
    end
    if any(sex ~= floor(sex) | ~ismember(sex, [0; 1; 2]))
        error('statgen:cache', 'Invalid genotype cache: sex values must be 0, 1, or 2');
    end
    if ~isequal(is_male, sex == 1)
        error('statgen:cache', 'Invalid genotype cache: is_male does not match sex');
    end
    if ~isequal(is_female, sex == 2)
        error('statgen:cache', 'Invalid genotype cache: is_female does not match sex');
    end
    if ~isequal(size(loaded.subject_present), [num_sample, numel(labels)])
        error('statgen:cache', 'Invalid genotype cache: subject_present shape mismatch');
    end
    if ~isequal(size(loaded.source_subject_row0), [num_sample, numel(labels)])
        error('statgen:cache', 'Invalid genotype cache: source_subject_row0 shape mismatch');
    end
    subject_present = logical(loaded.subject_present);
    source_subject_row0 = double(loaded.source_subject_row0);

    selected = statgen.internal.validate_requested_shards(shards, labels, 'load_genotype_cache');
    out_shards = cell(numel(selected), 1);
    for i = 1:numel(selected)
        label = selected{i};
        idx = find(strcmp(labels, label), 1, 'first');
        ix = (start0(idx) + 1):stop0(idx);
        out_shards{i} = statgen.GenotypeShard( ...
            labels{idx}, bed_paths{idx}, bed_file_sizes(idx), ...
            source_num_snp(idx), source_num_sample(idx), ...
            source_row0(ix), subject_present(:, idx), source_subject_row0(:, idx), ...
            is_present(ix), ploidy_male(ix), ploidy_female(ix), checksums{idx});
    end

    panel = statgen.GenotypePanel(out_shards, fid, iid, father_id, mother_id, sex, source_layout);
end

function require_fields_(loaded, fields)
    for i = 1:numel(fields)
        if ~isfield(loaded, fields{i})
            error('statgen:cache', 'Invalid genotype cache: missing field %s', fields{i});
        end
    end
end

function [labels, checksums, start0, stop0, source_layout, bed_paths, bed_file_sizes, source_num_snp, source_num_sample, num_sample] = validate_metadata_(meta, num_snp)
    if ~isfield(meta, 'schema') || ~strcmp(meta.schema, 'genotype_cache/0.1')
        if isfield(meta, 'schema')
            schema = meta.schema;
        else
            schema = '<missing>';
        end
        error('statgen:cache', 'Unsupported genotype cache schema: %s', schema);
    end
    need = {'n_shards', 'shard_labels', 'shard_checksums', 'shard_start0', ...
        'shard_stop0', 'source_layout', 'bed_paths', 'bed_file_sizes', ...
        'source_num_snp', 'source_num_sample', 'num_sample'};
    for i = 1:numel(need)
        if ~isfield(meta, need{i})
            error('statgen:cache', 'Invalid genotype cache: missing metadata field %s', need{i});
        end
    end

    labels = statgen.internal.ensure_cell_col(meta.shard_labels);
    checksums = statgen.internal.ensure_cell_col(meta.shard_checksums);
    start0 = double(meta.shard_start0(:));
    stop0 = double(meta.shard_stop0(:));
    bed_paths = statgen.internal.ensure_cell_col(meta.bed_paths);
    bed_file_sizes = double(meta.bed_file_sizes(:));
    source_num_snp = double(meta.source_num_snp(:));
    source_num_sample = double(meta.source_num_sample(:));
    source_layout = char(meta.source_layout);
    num_sample = double(meta.num_sample);

    n_shards = numel(labels);
    if ~isscalar(meta.n_shards) || double(meta.n_shards) ~= n_shards
        error('statgen:cache', 'Invalid genotype cache: n_shards mismatch');
    end
    if numel(checksums) ~= n_shards || numel(start0) ~= n_shards || numel(stop0) ~= n_shards || ...
            numel(bed_paths) ~= n_shards || numel(bed_file_sizes) ~= n_shards || ...
            numel(source_num_snp) ~= n_shards || numel(source_num_sample) ~= n_shards
        error('statgen:cache', 'Invalid genotype cache: shard metadata length mismatch');
    end
    statgen.internal.validate_requested_shards(labels, labels, 'load_genotype_cache');
    if ~strcmp(source_layout, 'non_sharded') && ~strcmp(source_layout, 'sharded')
        error('statgen:cache', 'Invalid genotype cache: source_layout must be ''non_sharded'' or ''sharded''');
    end
    if ~isscalar(num_sample) || num_sample < 0 || num_sample ~= floor(num_sample)
        error('statgen:cache', 'Invalid genotype cache: num_sample must be a non-negative integer');
    end

    validate_offsets_(start0, stop0, num_snp);
    for i = 1:n_shards
        expected = statgen.internal.bfile_expected_bed_size(source_num_sample(i), source_num_snp(i));
        if bed_file_sizes(i) ~= expected
            error('statgen:cache', 'Invalid genotype cache: bed_file_size mismatch for shard %s', labels{i});
        end
    end
end

function validate_offsets_(start0, stop0, num_snp)
    if isempty(start0)
        if num_snp ~= 0
            error('statgen:cache', 'Invalid genotype cache: empty shard offsets for non-empty payload');
        end
        return
    end
    if start0(1) ~= 0 || stop0(end) ~= num_snp || any(stop0 <= start0) || any(start0(2:end) ~= stop0(1:end-1))
        error('statgen:cache', 'Invalid genotype cache: shard offsets are not contiguous');
    end
end
