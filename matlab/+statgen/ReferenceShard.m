classdef ReferenceShard
% Immutable in-memory representation of one .bim shard.
    properties (SetAccess = private)
        label      % chromosome label string
        num_snp    % scalar integer
        bp         % n×1 double vector
        a1_hash64  % n×1 uint64 vector
        a2_hash64  % n×1 uint64 vector
        checksum   % lowercase MD5 hex string
    end
    properties (Dependent)
        chr        % n×1 cell array of strings
        snp        % n×1 cell array of strings
        a1         % n×1 cell array of strings
        a2         % n×1 cell array of strings
        is_single_nucleotide_variant  % n×1 logical vector
        is_strand_ambiguous     % n×1 logical vector
    end
    properties (Access = private)
        chr_data
        snp_data
        a1_data
        a2_data
        is_thin = false
    end

    methods
        function obj = ReferenceShard(label, chr_vec, snp_vec, bp_vec, a1_vec, a2_vec, checksum)
            if nargin == 0, return; end
            if nargin < 7, checksum = ''; end
            obj.label   = char(label);
            obj.num_snp = numel(chr_vec);
            obj.chr_data = statgen.internal.ensure_cell_col(chr_vec);
            obj.snp_data = statgen.internal.ensure_cell_col(snp_vec);
            obj.bp      = double(bp_vec(:));
            obj.a1_data = statgen.internal.ensure_cell_col(a1_vec);
            obj.a2_data = statgen.internal.ensure_cell_col(a2_vec);
            if numel(obj.chr_data) ~= obj.num_snp || numel(obj.snp_data) ~= obj.num_snp || numel(obj.bp) ~= obj.num_snp || numel(obj.a1_data) ~= obj.num_snp || numel(obj.a2_data) ~= obj.num_snp
                error('statgen:reference', 'ReferenceShard vector lengths must match');
            end
            bad = ~strcmp(obj.chr_data, obj.label);
            if any(bad)
                error('statgen:reference', 'Reference shard %s contains multiple chr labels', obj.label);
            end
            validate_distinct_alleles_(obj.label, obj.a1_data, obj.a2_data);
            obj.a1_hash64 = statgen.internal.allele_hash64(obj.a1_data);
            obj.a2_hash64 = statgen.internal.allele_hash64(obj.a2_data);
            if ~isempty(checksum)
                obj.checksum = char(checksum);
            else
                obj.checksum = reference_checksum_(obj.chr_data, obj.bp, obj.a1_data, obj.a2_data);
            end
        end

        function out = get.chr(obj)
            if obj.is_thin
                out = repmat({obj.label}, obj.num_snp, 1);
            else
                out = obj.chr_data;
            end
        end

        function out = get.snp(obj)
            if obj.is_thin
                error('statgen:cache', 'Reference field snp is unavailable in thin reference cache; load a full reference cache');
            end
            out = obj.snp_data;
        end

        function out = get.a1(obj)
            if obj.is_thin
                error('statgen:cache', 'Reference field a1 is unavailable in thin reference cache; load a full reference cache');
            end
            out = obj.a1_data;
        end

        function out = get.a2(obj)
            if obj.is_thin
                error('statgen:cache', 'Reference field a2 is unavailable in thin reference cache; load a full reference cache');
            end
            out = obj.a2_data;
        end

        function out = get.is_single_nucleotide_variant(obj)
            [out, ~] = statgen.internal.reference_hash_masks(obj.a1_hash64, obj.a2_hash64);
        end

        function out = get.is_strand_ambiguous(obj)
            [~, out] = statgen.internal.reference_hash_masks(obj.a1_hash64, obj.a2_hash64);
        end
    end

    methods (Static)
        function obj = from_thin(label, num_snp, bp_vec, a1_hash64, a2_hash64, checksum)
            obj = statgen.ReferenceShard();
            obj.label = char(label);
            obj.num_snp = double(num_snp);
            obj.bp = double(bp_vec(:));
            obj.a1_hash64 = uint64(a1_hash64(:));
            obj.a2_hash64 = uint64(a2_hash64(:));
            obj.checksum = char(checksum);
            obj.is_thin = true;
            if numel(obj.bp) ~= obj.num_snp || numel(obj.a1_hash64) ~= obj.num_snp || numel(obj.a2_hash64) ~= obj.num_snp
                error('statgen:cache', 'Invalid reference cache: thin payload vector lengths mismatch');
            end
            same_hash = obj.a1_hash64 == obj.a2_hash64;
            if any(same_hash)
                idx = find(same_hash, 1, 'first');
                error('statgen:cache', ...
                    'Invalid reference cache: shard %s variant %d has equal allele hashes for a1 and a2', ...
                    obj.label, idx);
            end
        end
    end
end

function validate_distinct_alleles_(label, a1, a2)
    same = strcmp(a1, a2);
    if any(same)
        idx = find(same, 1, 'first');
        error('statgen:reference', 'Reference shard %s variant %d: a1 and a2 must differ', label, idx);
    end
end

function out = reference_checksum_(chr, bp, a1, a2)
    % MD5 over 'chr:bp:a1:a2\n' lines in row order
    bp_str = cellstr(num2str(round(bp), '%d'));
    parts = strcat(chr, {':'}, bp_str, {':'}, a1, {':'}, a2, {sprintf('\n')});
    text_payload = [parts{:}];
    out = statgen.internal.md5_hex(text_payload);
end
