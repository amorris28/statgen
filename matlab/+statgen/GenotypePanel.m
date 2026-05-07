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

        function out = fetch_genotypes_int8(obj, snp_indices, varargin) %#ok<INUSD>
            error('statgen:genotype', 'GenotypePanel.fetch_genotypes_int8 is not implemented in MATLAB/Octave yet');
        end

        function out = fetch_genotypes(obj, snp_indices, varargin) %#ok<INUSD>
            error('statgen:genotype', 'GenotypePanel.fetch_genotypes is not implemented in MATLAB/Octave yet');
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
