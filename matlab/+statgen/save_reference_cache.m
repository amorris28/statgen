function save_reference_cache(panel, path, varargin)
% Save a ReferencePanel to a MATLAB binary .mat cache file.
%
%   statgen.save_reference_cache(panel, path)
%   statgen.save_reference_cache(panel, path, 'mode', 'full')
%   statgen.save_reference_cache(panel, path, 'mode', 'thin')
%   statgen.save_reference_cache(panel, path, 'format', 'v7')
    path = char(path);
    [mode, save_arg] = parse_args_(varargin{:});

    n_shards = numel(panel.shards);
    metadata.schema = 'reference_cache/0.1';
    metadata.mode = mode;
    metadata.n_shards = n_shards;
    metadata.shard_labels = cell(n_shards, 1);
    metadata.shard_checksums = cell(n_shards, 1);
    metadata.shard_start0 = zeros(n_shards, 1);
    metadata.shard_stop0 = zeros(n_shards, 1);

    for i = 1:n_shards
        s = panel.shards{i};
        off = panel.shard_offsets(i);
        metadata.shard_labels{i} = s.label;
        metadata.shard_checksums{i} = s.checksum;
        metadata.shard_start0(i) = off.start0;
        metadata.shard_stop0(i) = off.stop0;
    end

    bp = panel.bp;
    a1_hash64 = panel.a1_hash64;
    a2_hash64 = panel.a2_hash64;

    if strcmp(mode, 'full')
        chr = panel.chr;
        snp = panel.snp;
        a1 = panel.a1;
        a2 = panel.a2;
        save(path, 'metadata', 'chr', 'snp', 'bp', 'a1', 'a2', 'a1_hash64', 'a2_hash64', save_arg);
    else
        save(path, 'metadata', 'bp', 'a1_hash64', 'a2_hash64', save_arg);
    end
end

function [mode, save_arg] = parse_args_(varargin)
    [~, save_arg] = statgen.internal.parse_mat_format( ...
        'save_reference_cache', 'v7', {'v7', 'v7.3', 'v5'});
    mode = 'full';

    if mod(numel(varargin), 2) ~= 0
        error('statgen:io', 'save_reference_cache options must be name-value pairs');
    end
    seen_format = false;
    seen_mode = false;
    for i = 1:2:numel(varargin)
        if ~(ischar(varargin{i}) || isstring(varargin{i}))
            error('statgen:io', 'save_reference_cache option names must be strings');
        end
        name = char(varargin{i});
        if strcmpi(name, 'format')
            if seen_format
                error('statgen:io', 'save_reference_cache format option specified more than once');
            end
            seen_format = true;
            [~, save_arg] = statgen.internal.parse_mat_format( ...
                'save_reference_cache', varargin{i + 1}, {'v7', 'v7.3', 'v5'});
        elseif strcmpi(name, 'mode')
            if seen_mode
                error('statgen:io', 'save_reference_cache mode option specified more than once');
            end
            seen_mode = true;
            mode = normalize_mode_(varargin{i + 1});
        else
            error('statgen:io', 'save_reference_cache unknown option: %s', name);
        end
    end
end

function mode = normalize_mode_(value)
    if ~(ischar(value) || isstring(value)) || isempty(value)
        error('statgen:io', 'save_reference_cache mode must be full or thin');
    end
    mode = lower(char(value));
    if ~any(strcmp(mode, {'full', 'thin'}))
        error('statgen:io', 'save_reference_cache mode must be full or thin');
    end
end
