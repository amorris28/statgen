function convert_ld_npz_to_mat(input_root, output_root, shard, varargin)
%CONVERT_LD_NPZ_TO_MAT Convert one Python LD shard to MATLAB .mat format.
%
%   statgen.convert_ld_npz_to_mat(input_root, output_root, shard)
%   statgen.convert_ld_npz_to_mat(input_root, output_root, shard, 'format', format)
%
% input_root is a Python LD distribution created by script/statgen_build_ld.py.
% output_root receives the converted .mat LD shard and bundled reference BIM.
% After converting all needed shards, run statgen.create_ld_mat_manifest.
% For shard 'X', all chrX sex-label shards present in input_root are converted
% together: female, male, and/or combined.
%
% See also statgen.create_ld_mat_manifest, statgen.load_ld,
% statgen.validate_ld_distribution.
    if nargin < 3 || isempty(shard) || ~(ischar(shard) || isstring(shard))
        error('statgen:ld', 'convert_ld_npz_to_mat requires a shard label');
    end
    [format, save_arg] = statgen.internal.parse_mat_format( ...
        'convert_ld_npz_to_mat', 'v7.3', {'v7.3', 'v5'}, varargin{:});
    statgen.internal.convert_ld_npz_to_mat(input_root, output_root, shard, format, save_arg);
end
