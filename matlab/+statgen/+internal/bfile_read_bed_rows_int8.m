function decoded = bfile_read_bed_rows_int8(path, source_rows0, source_num_sample)
%BFILE_READ_BED_ROWS_INT8 Read selected SNP-major PLINK BED rows as int8 calls.
    path = char(path);
    source_rows0 = double(source_rows0(:));
    source_num_sample = double(source_num_sample);
    bytes_per_snp = ceil(source_num_sample / 4);
    decoded = int8(zeros(source_num_sample, numel(source_rows0)));
    lookup = bed_lookup_int8_();

    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open BED file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));

    % PERF: loop over requested SNPs retained; each SNP is a separate packed BED
    %       row requiring an individual seek. Inner decode is vectorized across subjects.
    for j = 1:numel(source_rows0)
        offset = 3 + source_rows0(j) * bytes_per_snp;
        status = fseek(fid, offset, 'bof');
        if status ~= 0
            error('statgen:bed', '%s: failed to seek to source SNP row %.0f', path, source_rows0(j));
        end
        packed = fread(fid, bytes_per_snp, 'uint8=>uint8');
        if numel(packed) ~= bytes_per_snp
            error('statgen:bed', ...
                '%s: short BED read for source SNP row %.0f: expected %.0f bytes, got %.0f', ...
                path, source_rows0(j), bytes_per_snp, numel(packed));
        end
        chunk = lookup(double(packed) + 1, :);
        values = reshape(chunk.', [], 1);
        decoded(:, j) = values(1:source_num_sample);
    end
    clear cleaner;
end

function lookup = bed_lookup_int8_()
    persistent cached_lookup
    if isempty(cached_lookup)
        values = uint16((0:255)');
        shifts = uint16([0, 2, 4, 6]);
        two_bit = bitand( ...
            bitshift(repmat(values, 1, 4), -double(repmat(shifts, 256, 1))), ...
            uint16(3));
        lookup = int8(zeros(256, 4));
        lookup(two_bit == uint16(0)) = int8(2);
        lookup(two_bit == uint16(1)) = int8(-1);
        lookup(two_bit == uint16(2)) = int8(1);
        lookup(two_bit == uint16(3)) = int8(0);
        cached_lookup = lookup;
    end
    lookup = cached_lookup;
end
