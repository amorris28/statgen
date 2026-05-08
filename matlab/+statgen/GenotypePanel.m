classdef GenotypePanel
% Ordered collection of GenotypeShard objects with lazy SNP-axis accessors.
    properties (SetAccess = private)
        num_snp
        num_sample
        shard_offsets
        shards
        fid
        iid
        father_id
        mother_id
        sex
        is_male
        is_female
        source_layout
    end
    properties (Dependent)
        is_present
        ploidy_male
        ploidy_female
        source_row0
    end

    methods
        function obj = GenotypePanel(shards_cell, fid, iid, father_id, mother_id, sex, source_layout)
            if nargin == 0, return; end
            source_layout = char(source_layout);
            if ~strcmp(source_layout, 'non_sharded') && ~strcmp(source_layout, 'sharded')
                error('statgen:genotype', 'source_layout must be ''non_sharded'' or ''sharded''');
            end
            if isempty(shards_cell)
                error('statgen:genotype', 'GenotypePanel requires at least one shard');
            end
            obj.source_layout = source_layout;
            obj.shards = shards_cell;
            obj.fid = statgen.internal.ensure_cell_col(fid);
            obj.iid = statgen.internal.ensure_cell_col(iid);
            obj.father_id = statgen.internal.ensure_cell_col(father_id);
            obj.mother_id = statgen.internal.ensure_cell_col(mother_id);
            obj.sex = double(sex(:));
            obj.num_sample = numel(obj.fid);
            if numel(obj.iid) ~= obj.num_sample || numel(obj.father_id) ~= obj.num_sample || numel(obj.mother_id) ~= obj.num_sample || numel(obj.sex) ~= obj.num_sample
                error('statgen:genotype', 'GenotypePanel FAM vector lengths must match');
            end
            if any(obj.sex ~= floor(obj.sex) | ~ismember(obj.sex, [0; 1; 2]))
                error('statgen:genotype', 'GenotypePanel sex values must be 0, 1, or 2');
            end
            obj.is_male = obj.sex == 1;
            obj.is_female = obj.sex == 2;

            n_shards = numel(shards_cell);
            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});
            total = 0;
            for i = 1:n_shards
                s = shards_cell{i};
                if numel(s.subject_present) ~= obj.num_sample
                    error('statgen:genotype', ...
                        'Genotype shard %s: subject_present length does not match panel sample axis', ...
                        s.label);
                end
                offsets(i).shard_label = s.label;
                offsets(i).start0 = total;
                offsets(i).stop0 = total + s.num_snp;
                total = total + s.num_snp;
            end
            obj.num_snp = total;
            obj.shard_offsets = offsets;
        end

        function out = get.is_present(obj)
            out = concat_required_(obj.shards, 'is_present');
        end

        function out = get.ploidy_male(obj)
            out = concat_required_(obj.shards, 'ploidy_male');
        end

        function out = get.ploidy_female(obj)
            out = concat_required_(obj.shards, 'ploidy_female');
        end

        function out = get.source_row0(obj)
            out = concat_required_(obj.shards, 'source_row0');
        end

        function out = is_subject_present(obj, shard)
            label = char(shard);
            for i = 1:numel(obj.shards)
                if strcmp(obj.shards{i}.label, label)
                    out = obj.shards{i}.subject_present;
                    return
                end
            end
            error('statgen:genotype', ...
                'GenotypePanel.is_subject_present: requested shard %s is not loaded', label);
        end

        function out = select_shards(obj, shards)
            available = cell(numel(obj.shards), 1);
            for i = 1:numel(obj.shards)
                available{i} = obj.shards{i}.label;
            end
            selected = statgen.internal.validate_requested_shards( ...
                shards, available, 'GenotypePanel.select_shards');

            out_shards = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_shards{i} = obj.shards{idx};
            end
            out = statgen.GenotypePanel(out_shards, obj.fid, obj.iid, ...
                obj.father_id, obj.mother_id, obj.sex, obj.source_layout);
        end

        function out = fetch_genotypes_int8(obj, snp_indices, varargin)
            bed_path = parse_bed_path_(varargin{:});
            indices0 = normalize_snp_indices_(snp_indices, obj.num_snp);
            out = int8(-ones(obj.num_sample, numel(indices0)));
            if isempty(indices0)
                return
            end

            addressed = addressed_shards_(obj, indices0);
            bed_paths = resolve_bed_paths_(obj, addressed, bed_path);
            validate_requested_present_(obj, addressed);

            for i = 1:numel(addressed)
                a = addressed{i};
                shard = obj.shards{a.shard_index};
                effective_bed = bed_paths{i};
                statgen.internal.bfile_validate_bed( ...
                    effective_bed, shard.source_num_sample, shard.source_num_snp);
                source_rows0 = shard.source_row0(a.local_indices0 + 1);
                source_geno = statgen.internal.bfile_read_bed_rows_int8( ...
                    effective_bed, source_rows0, shard.source_num_sample);
                panel_rows = find(shard.subject_present);
                if ~isempty(panel_rows)
                    source_subject_rows = shard.source_subject_row0(panel_rows) + 1;
                    out(panel_rows, a.columns) = source_geno(source_subject_rows, :);
                end
            end
        end

        function out = fetch_genotypes(obj, snp_indices, varargin)
            geno_int8 = obj.fetch_genotypes_int8(snp_indices, varargin{:});
            out = double(geno_int8);
            out(geno_int8 == int8(-1)) = NaN;
        end

        function save_cache(obj, path, varargin)
            statgen.save_genotype_cache(obj, path, varargin{:});
        end
    end
end

function out = concat_required_(shards, field_name)
    if isempty(shards)
        out = [];
        return
    end
    vals = cell(numel(shards), 1);
    for i = 1:numel(shards)
        vals{i} = shards{i}.(field_name);
    end
    out = vertcat(vals{:});
end

function bed_path = parse_bed_path_(varargin)
    if nargin == 0
        bed_path = [];
    elseif nargin == 1
        bed_path = varargin{1};
    else
        error('statgen:genotype', ...
            'GenotypePanel.fetch_genotypes accepts at most one optional bed_path argument');
    end
    if isstring(bed_path)
        if ~isscalar(bed_path)
            error('statgen:genotype', 'GenotypePanel.fetch_genotypes: bed_path must be a string scalar');
        end
        bed_path = char(bed_path);
    elseif ~isempty(bed_path) && ~ischar(bed_path)
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: bed_path must be a character vector');
    end
end

function indices0 = normalize_snp_indices_(snp_indices, num_snp)
    if isempty(snp_indices)
        indices0 = zeros(0, 1);
        return
    end
    if ~isnumeric(snp_indices)
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: snp_indices must be integer indices');
    end
    indices = double(snp_indices(:));
    if any(~isfinite(indices)) || any(indices ~= floor(indices))
        error('statgen:genotype', 'GenotypePanel.fetch_genotypes: snp_indices must be integer indices');
    end
    if any(indices < 1) || any(indices > num_snp)
        bad = indices(find(indices < 1 | indices > num_snp, 1, 'first'));
        error('statgen:genotype', ...
            'GenotypePanel.fetch_genotypes: SNP index %.0f is out of bounds for num_snp=%.0f', ...
            bad, num_snp);
    end
    indices0 = indices - 1;
end

function addressed = addressed_shards_(obj, indices0)
    addressed = {};
    for i = 1:numel(obj.shards)
        offset = obj.shard_offsets(i);
        cols = find(indices0 >= offset.start0 & indices0 < offset.stop0);
        if ~isempty(cols)
            a.shard_index = i;
            a.columns = cols(:)';
            a.local_indices0 = indices0(cols) - offset.start0;
            addressed{end + 1, 1} = a; %#ok<AGROW>
        end
    end
end

function bed_paths = resolve_bed_paths_(obj, addressed, bed_path)
    bed_paths = cell(numel(addressed), 1);
    if isempty(bed_path)
        for i = 1:numel(addressed)
            bed_paths{i} = obj.shards{addressed{i}.shard_index}.bed_path;
        end
        return
    end

    if ~isempty(strfind(bed_path, '@'))
        if strcmp(obj.source_layout, 'non_sharded')
            error('statgen:genotype', ...
                'GenotypePanel.fetch_genotypes: @ override incompatible with non-sharded panel metadata');
        end
        for i = 1:numel(addressed)
            shard = obj.shards{addressed{i}.shard_index};
            bed_paths{i} = strrep(bed_path, '@', shard.label);
        end
        return
    end

    if strcmp(obj.source_layout, 'sharded') && ~isempty(addressed)
        source_num_snp = zeros(numel(addressed), 1);
        source_num_sample = zeros(numel(addressed), 1);
        for i = 1:numel(addressed)
            shard = obj.shards{addressed{i}.shard_index};
            source_num_snp(i) = shard.source_num_snp;
            source_num_sample(i) = shard.source_num_sample;
        end
        if numel(unique(source_num_snp)) ~= 1 || numel(unique(source_num_sample)) ~= 1
            error('statgen:genotype', ...
                ['GenotypePanel.fetch_genotypes: flat bed_path override requires equal ' ...
                 'source_num_snp and source_num_sample across addressed shards']);
        end
    end

    for i = 1:numel(addressed)
        bed_paths{i} = bed_path;
    end
end

function validate_requested_present_(obj, addressed)
    for i = 1:numel(addressed)
        a = addressed{i};
        shard = obj.shards{a.shard_index};
        present = shard.is_present(a.local_indices0 + 1);
        if ~all(present)
            j = find(~present, 1, 'first');
            panel_index = obj.shard_offsets(a.shard_index).start0 + a.local_indices0(j) + 1;
            error('statgen:genotype', ...
                'GenotypePanel.fetch_genotypes: requested SNP %.0f in shard %s is not present in the genotype source', ...
                panel_index, shard.label);
        end
    end
end
