classdef SumstatsShard
% Immutable in-memory aligned sumstats vectors for one reference shard.
    properties (SetAccess = private)
        label
        reference_checksum
        num_snp
        zvec
        nvec
        logpvec
        beta_vec
        se_vec
        eaf_vec
        info_vec
        is_present
    end

    methods
        function obj = SumstatsShard(label, reference_checksum, logpvec, zvec, nvec, beta_vec, se_vec, eaf_vec, info_vec)
            if nargin == 0, return; end
            obj.label = char(label);
            obj.reference_checksum = char(reference_checksum);
            obj.logpvec = double(logpvec(:));
            if nargin < 4 || isempty(zvec), obj.zvec = []; else, obj.zvec = double(zvec(:)); end
            if nargin < 5 || isempty(nvec), obj.nvec = []; else, obj.nvec = double(nvec(:)); end
            if nargin < 6 || isempty(beta_vec), obj.beta_vec = []; else, obj.beta_vec = double(beta_vec(:)); end
            if nargin < 7 || isempty(se_vec),   obj.se_vec = [];   else, obj.se_vec = double(se_vec(:));   end
            if nargin < 8 || isempty(eaf_vec),  obj.eaf_vec = [];  else, obj.eaf_vec = double(eaf_vec(:)); end
            if nargin < 9 || isempty(info_vec), obj.info_vec = []; else, obj.info_vec = double(info_vec(:)); end
            obj.num_snp = numel(obj.logpvec);
            validate_optional_length_(obj.zvec, obj.num_snp, 'zvec');
            validate_optional_length_(obj.nvec, obj.num_snp, 'nvec');
            validate_optional_length_(obj.beta_vec, obj.num_snp, 'beta_vec');
            validate_optional_length_(obj.se_vec, obj.num_snp, 'se_vec');
            validate_optional_length_(obj.eaf_vec, obj.num_snp, 'eaf_vec');
            validate_optional_length_(obj.info_vec, obj.num_snp, 'info_vec');
            obj.is_present = ~isnan(obj.logpvec);
        end
    end
end

function validate_optional_length_(vec, n, name)
    if ~isempty(vec) && numel(vec) ~= n
        error('statgen:sumstats', '%s length mismatch: expected %d, got %d', name, n, numel(vec));
    end
end
