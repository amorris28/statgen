function convert_ld_npz_to_mat(input_root, output_root, shard, varargin)
% Convert one Python LD .npz reference shard into MATLAB/Octave .mat files.
    if nargin < 3 || isempty(shard) || ~(ischar(shard) || isstring(shard))
        error('statgen:ld', 'convert_ld_npz_to_mat requires a shard label');
    end
    [format, save_arg] = statgen.internal.parse_mat_format( ...
        'convert_ld_npz_to_mat', 'v7.3', {'v7.3', 'v5'}, varargin{:});
    statgen.internal.convert_ld_npz_to_mat(input_root, output_root, shard, format, save_arg);
end
