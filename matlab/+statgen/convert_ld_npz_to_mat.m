function convert_ld_npz_to_mat(input_root, output_root, shard)
% Convert one Python LD .npz reference shard into production MATLAB .mat files.
    if nargin < 3 || isempty(shard) || ~(ischar(shard) || isstring(shard))
        error('statgen:ld', 'convert_ld_npz_to_mat requires a shard label');
    end
    statgen.internal.convert_ld_npz_to_mat(input_root, output_root, shard, true);
end
