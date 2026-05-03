classdef LDShard
% Immutable raw LD shard loaded from a runtime-native distribution artifact.
    properties (SetAccess = private)
        chr
        label
        sex
        num_snp
        ld_r
        a1freq
        reference_checksum
        checksum
    end

    methods
        function obj = LDShard(chr_label, sex, num_snp, ld_r, a1freq, reference_checksum)
            if nargin == 0, return; end
            obj.chr = char(chr_label);
            obj.label = char(chr_label);
            if isempty(sex)
                obj.sex = [];
            else
                obj.sex = char(sex);
            end
            obj.num_snp = double(num_snp);
            obj.ld_r = sparse(double(ld_r));
            obj.a1freq = double(a1freq(:));
            obj.reference_checksum = char(reference_checksum);
            obj.checksum = obj.reference_checksum;
        end
    end
end
