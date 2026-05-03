classdef LDPanel
% Ordered collection of LD shard groups matching a ReferencePanel.
    properties (SetAccess = private)
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
