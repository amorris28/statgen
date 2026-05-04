function selected = ld_selected_shards(ld_panel, chrX_sex)
% Select one LD shard per reference shard group, applying chrX sex semantics.
    if nargin < 2
        chrX_sex = [];
    end

    groups = ld_panel.shard_groups;
    selected = cell(numel(groups), 1);
    for i = 1:numel(groups)
        group = groups{i};
        first = group{1};
        if ~strcmp(first.chr, 'X')
            selected{i} = first;
            continue
        end

        if isempty(chrX_sex)
            sex = ld_panel.default_chrX_sex;
        else
            statgen.LDPanel.validate_chrx_sex_(chrX_sex, 'chrX_sex');
            sex = char(chrX_sex);
        end

        found = [];
        present = cell(numel(group), 1);
        for j = 1:numel(group)
            present{j} = group{j}.sex;
            if strcmp(group{j}.sex, sex)
                found = group{j};
            end
        end
        if isempty(found)
            error('statgen:ld', ...
                'chrX_sex must name a loaded chrX LD shard; got %s, present %s', ...
                sex, strjoin(present, ', '));
        end
        selected{i} = found;
    end
end
