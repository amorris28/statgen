function ld_validate_chr_sex(chr_label, sex, where)
    canonical = statgen.internal.canonical_labels();
    if ~any(strcmp(char(chr_label), canonical))
        error('statgen:ld', '%s: chr must be one of 1-22 or X', where);
    end
    if strcmp(char(chr_label), 'X')
        valid = {'female', 'male', 'combined'};
        if isempty(sex) || ~any(strcmp(char(sex), valid))
            error('statgen:ld', '%s: sex: chrX sex must be one of female, male, combined', where);
        end
    elseif ~isempty(sex)
        error('statgen:ld', '%s: autosomal LD shards must have sex empty', where);
    end
end
