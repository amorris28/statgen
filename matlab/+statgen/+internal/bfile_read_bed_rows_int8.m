function decoded = bfile_read_bed_rows_int8(path, source_rows0, source_num_sample)
%BFILE_READ_BED_ROWS_INT8 Read selected SNP-major PLINK BED rows as int8 calls.
    path = char(path);
    source_rows0 = double(source_rows0(:));
    source_num_sample = double(source_num_sample);
    bytes_per_snp = ceil(source_num_sample / 4);
    lookup = bed_lookup_int8_();
    if isempty(source_rows0)
        decoded = int8(zeros(source_num_sample, 0));
        return;
    end
    [unique_rows0, ~, inverse] = unique(source_rows0, 'sorted');
    unique_decoded = int8(zeros(source_num_sample, numel(unique_rows0)));
    max_read_bytes = 64 * 1024^2;
    if bytes_per_snp == 0
        rows_per_read = numel(unique_rows0);
    else
        rows_per_read = max(1, floor(max_read_bytes / bytes_per_snp));
    end

    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open BED file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));

    group_starts = [1; find(diff(unique_rows0) ~= 1) + 1; numel(unique_rows0) + 1];
    % PERF: one seek/read per contiguous BED row block; requested order and
    %       repeats are restored after block decoding. Blocks are capped to
    %       bound transient packed/decode allocations for large requests.
    for g = 1:(numel(group_starts) - 1)
        group_cols = group_starts(g):(group_starts(g + 1) - 1);
        for chunk_start = 1:rows_per_read:numel(group_cols)
            chunk_stop = min(chunk_start + rows_per_read - 1, numel(group_cols));
            block_cols = group_cols(chunk_start:chunk_stop);
            block_rows0 = unique_rows0(block_cols);
            offset = 3 + block_rows0(1) * bytes_per_snp;
            expected_bytes = bytes_per_snp * numel(block_rows0);
            status = fseek(fid, offset, 'bof');
            if status ~= 0
                error('statgen:bed', '%s: failed to seek to source SNP row %.0f', path, block_rows0(1));
            end
            packed = fread(fid, expected_bytes, 'uint8=>uint8');
            if numel(packed) ~= expected_bytes
                error('statgen:bed', ...
                    '%s: short BED read for source SNP rows %.0f-%.0f: expected %.0f bytes, got %.0f', ...
                    path, block_rows0(1), block_rows0(end), expected_bytes, numel(packed));
            end
            chunk = lookup(double(packed) + 1, :);
            values = reshape(chunk.', [], numel(block_rows0));
            unique_decoded(:, block_cols) = values(1:source_num_sample, :);
        end
    end
    decoded = unique_decoded(:, inverse);
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
