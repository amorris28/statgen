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
    end
    properties (SetAccess = private, Hidden)
        ld_r2
        retain_ld_r
    end

    methods
        function obj = LDShard(chr_label, sex, num_snp, ld_r, a1freq, reference_checksum, retain_ld_r)
            if nargin == 0, return; end
            if nargin < 7 || isempty(retain_ld_r)
                retain_ld_r = true;
            end
            obj.chr = char(chr_label);
            obj.label = char(chr_label);
            if isempty(sex)
                obj.sex = [];
            else
                obj.sex = char(sex);
            end
            obj.num_snp = double(num_snp);
            ld_r = sparse(double(ld_r));
            if retain_ld_r
                obj.ld_r = ld_r;
            else
                obj.ld_r = [];
            end
            obj.retain_ld_r = logical(retain_ld_r);
            obj.ld_r2 = spfun(@(x) x .^ 2, ld_r);
            obj.a1freq = double(a1freq(:));
            obj.reference_checksum = char(reference_checksum);
        end
    end
end
