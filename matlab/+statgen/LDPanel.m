classdef LDPanel
% Ordered collection of LD shard groups matching a ReferencePanel.
    properties (SetAccess = private)
        num_snp
        shard_offsets
        shard_groups
        default_chrX_sex
    end
    properties (Dependent)
        shards
    end

    methods
        function obj = LDPanel(shard_groups, default_chrX_sex)
            if nargin == 0, return; end
            if nargin < 2 || isempty(default_chrX_sex)
                default_chrX_sex = 'female';
            end
            statgen.LDPanel.validate_chrx_sex_(default_chrX_sex, 'default_chrX_sex');
            obj.shard_groups = shard_groups;
            obj.default_chrX_sex = char(default_chrX_sex);

            total = 0;
            offsets = struct('shard_label', {}, 'start0', {}, 'stop0', {});
            for i = 1:numel(obj.shard_groups)
                group = obj.shard_groups{i};
                if isempty(group)
                    error('statgen:ld', 'LDPanel shard groups must be non-empty');
                end
                first = group{1};
                for j = 1:numel(group)
                    if ~strcmp(group{j}.chr, first.chr)
                        error('statgen:ld', 'all LD shards in a group must share chr');
                    end
                    if group{j}.num_snp ~= first.num_snp
                        error('statgen:ld', 'all LD shards in a group must share num_snp');
                    end
                end
                offsets(i).shard_label = first.label;
                offsets(i).start0 = total;
                offsets(i).stop0 = total + first.num_snp;
                total = total + first.num_snp;
            end
            obj.num_snp = total;
            obj.shard_offsets = offsets;
            obj.validate_default_chrX_sex_();
        end

        function out = get.shards(obj)
            out = cell(numel(obj.shard_groups), 1);
            for i = 1:numel(obj.shard_groups)
                group = obj.shard_groups{i};
                if isempty(group)
                    out{i} = [];
                    continue
                end
                first = group{1};
                if strcmp(first.chr, 'X')
                    found = [];
                    for j = 1:numel(group)
                        if strcmp(group{j}.sex, obj.default_chrX_sex)
                            found = group{j};
                            break
                        end
                    end
                    out{i} = found;
                else
                    out{i} = first;
                end
            end
        end

        function out = a1freq(obj, chrX_sex)
            if nargin < 2
                chrX_sex = [];
            end
            selected = statgen.internal.ld_selected_shards(obj, chrX_sex);
            if isempty(selected)
                out = [];
                return
            end
            vals = cell(numel(selected), 1);
            for i = 1:numel(selected)
                vals{i} = selected{i}.a1freq;
            end
            out = vertcat(vals{:});
        end

        function out = select_shards(obj, shards)
            available = cell(numel(obj.shard_groups), 1);
            for i = 1:numel(obj.shard_groups)
                available{i} = obj.shard_groups{i}{1}.label;
            end
            selected = statgen.internal.validate_requested_shards( ...
                shards, available, 'LDPanel.select_shards');
            out_groups = cell(numel(selected), 1);
            for i = 1:numel(selected)
                idx = find(strcmp(available, selected{i}), 1, 'first');
                out_groups{i} = obj.shard_groups{idx};
            end
            out = statgen.LDPanel(out_groups, obj.default_chrX_sex);
        end

        function out = multiply_r2(obj, M, chrX_sex)
            if nargin < 3
                chrX_sex = [];
            end

            original_size = size(M);
            vector_input = isvector(M);
            row_vector = isrow(M) && numel(M) == obj.num_snp;
            if vector_input
                if numel(M) ~= obj.num_snp
                    error('statgen:ld', ...
                        'multiply_r2: vector length mismatch: expected %d, got %d', ...
                        obj.num_snp, numel(M));
                end
                work = M(:);
            else
                if ndims(M) ~= 2 || size(M, 1) ~= obj.num_snp
                    error('statgen:ld', ...
                        'multiply_r2: matrix row count mismatch: expected %d', obj.num_snp);
                end
                work = M;
            end

            selected = statgen.internal.ld_selected_shards(obj, chrX_sex);
            parts = cell(numel(selected), 1);
            for i = 1:numel(selected)
                offset = obj.shard_offsets(i);
                rows = (offset.start0 + 1):offset.stop0;
                if numel(rows) ~= selected{i}.num_snp
                    error('statgen:ld', 'multiply_r2: shard offset length does not match LD shard size');
                end
                parts{i} = statgen.internal.ld_multiply_r2(selected{i}.ld_r2, work(rows, :));
            end
            out = vertcat(parts{:});

            if vector_input
                if row_vector
                    out = reshape(out, original_size);
                else
                    out = out(:);
                end
            end
        end
    end

    methods (Access = private)
        function validate_default_chrX_sex_(obj)
            for i = 1:numel(obj.shard_groups)
                group = obj.shard_groups{i};
                if isempty(group) || ~strcmp(group{1}.chr, 'X')
                    continue
                end
                present = cell(numel(group), 1);
                for j = 1:numel(group)
                    present{j} = group{j}.sex;
                end
                if ~any(strcmp(present, obj.default_chrX_sex))
                    error('statgen:ld', ...
                        'default_chrX_sex must name a loaded chrX LD shard');
                end
            end
        end

    end

    methods (Static)
        function validate_chrx_sex_(sex, where)
            valid = {'female', 'male', 'combined'};
            if isempty(sex) || ~any(strcmp(char(sex), valid))
                error('statgen:ld', '%s: chrX sex must be one of female, male, combined', where);
            end
        end
    end
end
