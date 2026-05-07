function fam = bfile_parse_fam(path)
%BFILE_PARSE_FAM Parse and validate a PLINK FAM sidecar.
    path = char(path);
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open FAM file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    raw = textscan(fid, '%s%s%s%s%f%f', ...
        'Delimiter', sprintf(' \t'), ...
        'MultipleDelimsAsOne', true, ...
        'ReturnOnError', false);
    clear cleaner;

    n = numel(raw{1});
    for c = 2:6
        if numel(raw{c}) ~= n
            error('statgen:fam', '%s: expected 6 whitespace-delimited columns', path);
        end
    end

    fam.fid = statgen.internal.ensure_cell_col(raw{1});
    fam.iid = statgen.internal.ensure_cell_col(raw{2});
    fam.father_id = statgen.internal.ensure_cell_col(raw{3});
    fam.mother_id = statgen.internal.ensure_cell_col(raw{4});
    fam.sex = double(raw{5}(:));

    empty = cellfun('isempty', fam.fid) | cellfun('isempty', fam.iid) | ...
        cellfun('isempty', fam.father_id) | cellfun('isempty', fam.mother_id) | isnan(fam.sex);
    if any(empty)
        row = find(empty, 1, 'first');
        error('statgen:fam', '%s:%d: FAM must contain exactly 6 non-empty columns', path, row);
    end

    bad_sex = isnan(fam.sex) | fam.sex ~= floor(fam.sex) | ~ismember(fam.sex, [0; 1; 2]);
    if any(bad_sex)
        row = find(bad_sex, 1, 'first');
        error('statgen:fam', '%s:%d: FAM sex must be one of PLINK values 1, 2, or 0', path, row);
    end

    keys = strcat(fam.fid, char(9), fam.iid);
    [~, first_idx] = unique(keys, 'stable');
    dup = true(n, 1);
    dup(first_idx) = false;
    if any(dup)
        row = find(dup, 1, 'first');
        error('statgen:fam', '%s:%d: duplicate FAM subject pair (fid, iid)', path, row);
    end
end
