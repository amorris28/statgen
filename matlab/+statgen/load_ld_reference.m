function reference = load_ld_reference(path, varargin)
%LOAD_LD_REFERENCE Load the ReferencePanel bundled with an LD distribution.
%
%   reference = statgen.load_ld_reference(path)
%   reference = statgen.load_ld_reference(path, shards)
%
% path is an LD distribution directory containing ld_manifest.json and the
% manifest-declared reference cache. Optional shards restrict loading to
% selected canonical shard labels. LD matrices are not loaded.
%
% See also statgen.load_ld, statgen.load_reference_cache.
    shards = parse_args_(varargin{:});
    path = char(path);
    if ~isfolder(path)
        error('statgen:ld', 'load_ld_reference: path must identify a panel root directory');
    end

    manifest = statgen.internal.ld_read_manifest(fullfile(path, 'ld_manifest.json'), ...
        'matlab_mat_sparse_double');
    reference = statgen.load_reference_cache(fullfile(path, manifest.reference_cache), shards);
end

function shards = parse_args_(varargin)
    shards = [];
    if numel(varargin) > 1
        error('statgen:arg', 'load_ld_reference accepts path and optional shards');
    end
    if ~isempty(varargin)
        shards = varargin{1};
    end
end
