function manifest = ld_read_manifest(path, expected_runtime_format)
    if nargin < 2
        expected_runtime_format = [];
    end
    if exist(path, 'file') ~= 2
        error('statgen:io', 'LD manifest not found: %s', path);
    end
    manifest = jsondecode(fileread(path));
    if ~isfield(manifest, 'object_type') || ~strcmp(manifest.object_type, 'ld_panel_manifest')
        error('statgen:ld', '%s: object_type must be ld_panel_manifest', path);
    end
    if ~isfield(manifest, 'schema_version') || ~strcmp(manifest.schema_version, '1.0')
        error('statgen:ld', '%s: unsupported schema_version', path);
    end
    if ~isfield(manifest, 'runtime_format')
        error('statgen:ld', '%s: missing runtime_format', path);
    end
    if ~isempty(expected_runtime_format) && ~strcmp(manifest.runtime_format, expected_runtime_format)
        error('statgen:ld', 'Expected MATLAB/Octave LD manifest runtime_format');
    end
    if ~isfield(manifest, 'shards') || isempty(manifest.shards)
        error('statgen:ld', '%s: shards must be a non-empty list', path);
    end
    for i = 1:numel(manifest.shards)
        validate_manifest_entry_(manifest.shards(i), sprintf('%s:shards(%d)', path, i));
    end
end

function validate_manifest_entry_(entry, where)
    required = {'chr', 'sex', 'file', 'file_md5', 'num_snp', 'nnz', 'reference_checksum'};
    for i = 1:numel(required)
        if ~isfield(entry, required{i})
            error('statgen:ld', '%s: missing required field %s', where, required{i});
        end
    end
    statgen.internal.ld_validate_chr_sex(entry.chr, entry.sex, where);
    if isempty(entry.file)
        error('statgen:ld', '%s: file must be non-empty', where);
    end
end
