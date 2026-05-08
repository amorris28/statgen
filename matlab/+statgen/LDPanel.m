classdef LDPanel
%STATGEN.LDPANEL Sparse LD matrices aligned to a ReferencePanel.
%
%   ld = statgen.load_ld(path)
%   ld = statgen.load_ld(path, reference)
%
% An LDPanel stores per-shard sparse LD matrices in reference-panel order.
% Use multiply_r2 to multiply by LD r-squared and statgen.fast_prune to prune
% a score vector.
%
% chrX LD may have female, male, and/or combined shards. default_chrX_sex
% selects which chrX shard is used when a method needs one LD matrix. Method
% arguments named chrX_sex override that default for chrX only; autosomes are
% always sex-agnostic.
%
% Common properties:
%   num_snp     Number of SNPs in the LD panel.
%   reference   Paired ReferencePanel, when available.
%   default_chrX_sex
%               Default chrX shard selector: 'female', 'male', or 'combined'.
%
% Common methods:
%   a1freq        Return allele frequencies in panel order.
%   select_shards Restrict the panel to selected shards.
%   multiply_r2   Multiply a vector or matrix by LD r-squared.
%
% See also statgen.load_ld, statgen.fast_prune,
% statgen.convert_ld_npz_to_mat, statgen.create_ld_mat_manifest.
    properties (SetAccess = private)
        num_snp
        shard_offsets
        shard_groups
        default_chrX_sex
        reference
    end
    properties (Dependent)
        shards
    end

    methods
        function obj = LDPanel(shard_groups, default_chrX_sex, reference)
            if nargin == 0, return; end
            if nargin < 2 || isempty(default_chrX_sex)
                default_chrX_sex = 'female';
            end
            if nargin < 3
                reference = [];
            end
            if isempty(shard_groups)
                error('statgen:ld', 'LDPanel requires at least one shard group');
            end
            statgen.LDPanel.validate_chrx_sex_(default_chrX_sex, 'default_chrX_sex');
            obj.shard_groups = shard_groups;
            obj.default_chrX_sex = char(default_chrX_sex);
            obj.reference = reference;

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
            obj.validate_reference_shape_();
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
        %A1FREQ Return LD allele frequencies in panel order.
        %
        %   f = ld.a1freq()
        %   f = ld.a1freq(chrX_sex)
        %
        % Returns a num_snp-by-1 vector. chrX_sex selects the chrX LD shard
        % when chrX is loaded; use 'female', 'male', or 'combined'.
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
        %SELECT_SHARDS Return an LD panel restricted to selected shards.
        %
        %   out = ld.select_shards(shards)
        %
        % shards is a cell array or string array of canonical shard labels, for
        % example {'21', '22'}. The returned LDPanel preserves the requested
        % shard order.
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
            out_reference = [];
            if ~isempty(obj.reference)
                out_reference = obj.reference.select_shards(selected);
            end
            out = statgen.LDPanel(out_groups, obj.default_chrX_sex, out_reference);
        end

        function out = multiply_r2(obj, M, chrX_sex)
        %MULTIPLY_R2 Multiply a vector or matrix by LD r-squared.
        %
        %   out = ld.multiply_r2(M)
        %   out = ld.multiply_r2(M, chrX_sex)
        %
        % M must have num_snp rows, or be a vector with num_snp elements.
        % A row-vector input returns a row vector; a column-vector input returns
        % a column vector. Matrix inputs keep their matrix shape. chrX_sex
        % selects the chrX LD shard when chrX is loaded; autosomes are
        % sex-agnostic.
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

        function validate_reference_shape_(obj)
            if isempty(obj.reference)
                return
            end
            if numel(obj.reference.shards) ~= numel(obj.shard_groups)
                error('statgen:ld', 'LDPanel reference shard count must match LD shard groups');
            end
            for i = 1:numel(obj.shard_groups)
                ref = obj.reference.shards{i};
                first = obj.shard_groups{i}{1};
                if ~strcmp(ref.label, first.label)
                    error('statgen:ld', 'LDPanel reference shard labels must match LD shard groups');
                end
                if ref.num_snp ~= first.num_snp
                    error('statgen:ld', 'LDPanel reference shard sizes must match LD shard groups');
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
