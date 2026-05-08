function [ploidy_male, ploidy_female] = bfile_parse_ploidy(path, source_num_snp, source_chr)
%BFILE_PARSE_PLOIDY Parse optional genomatch BFILE .ploidy sidecar.
    if nargin < 3
        source_chr = [];
    end
    if isempty(path)
        ploidy_male = 2 * ones(source_num_snp, 1);
        ploidy_female = 2 * ones(source_num_snp, 1);
        if ~isempty(source_chr)
            if numel(source_chr) ~= source_num_snp
                error('statgen:ploidy', 'source_chr length must match source_num_snp');
            end
            ploidy_male(strcmp(cellstr(source_chr(:)), 'X')) = 1;
        end
        return
    end

    path = char(path);
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open PLOIDY file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    raw = textscan(fid, '%f%f', ...
        'Delimiter', sprintf('\t'), ...
        'Whitespace', '', ...
        'MultipleDelimsAsOne', false, ...
        'ReturnOnError', false);
    clear cleaner;

    n = numel(raw{1});
    if numel(raw{2}) ~= n
        error('statgen:ploidy', '%s: expected 2 tab-delimited columns', path);
    end
    if n ~= source_num_snp
        error('statgen:ploidy', ...
            '%s: PLOIDY row count mismatch: expected %d, got %d', path, source_num_snp, n);
    end

    ploidy_male = double(raw{1}(:));
    ploidy_female = double(raw{2}(:));
    bad = isnan(ploidy_male) | isnan(ploidy_female) | ...
        ploidy_male ~= floor(ploidy_male) | ploidy_female ~= floor(ploidy_female) | ...
        ~ismember(ploidy_male, [0; 1; 2]) | ~ismember(ploidy_female, [0; 1; 2]);
    if any(bad)
        row = find(bad, 1, 'first');
        error('statgen:ploidy', '%s:%d: ploidy values must be 0, 1, or 2', path, row);
    end
end
