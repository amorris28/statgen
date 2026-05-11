function out = ld_md5_file(path)
    path = char(path);
    if exist(path, 'file') ~= 2
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    out = md5_file_system_(path);
    if ~isempty(out)
        return
    end

    out = md5_file_stream_(path);
end

function out = md5_file_stream_(path)
    try
        md = java.security.MessageDigest.getInstance('MD5');
    catch ME
        out = md5_file_hash_full_(path, ME);
        return
    end

    warning('statgen:ld:md5Fallback', ...
        '%s: system MD5 command unavailable; computing MD5 in MATLAB', path);

    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    chunk_size = 16 * 1024 * 1024;
    while true
        bytes = fread(fid, chunk_size, '*uint8');
        if isempty(bytes)
            break
        end
        md.update(uint8(bytes(:)));
    end
    clear cleaner;

    d = typecast(md.digest(), 'uint8');
    out = lower(reshape(dec2hex(d)', 1, []));
end

function out = md5_file_hash_full_(path, java_error)
    if exist('hash', 'builtin') ~= 5 && exist('hash', 'file') ~= 2
        error('statgen:io', ...
            'Cannot compute MD5 for %s without md5sum/md5 command, Java MessageDigest, or hash(): %s', ...
            path, java_error.message);
    end

    warning('statgen:ld:md5Fallback', ...
        '%s: system MD5 command unavailable; computing MD5 in MATLAB', path);

    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
    clear cleaner;
    out = lower(hash('md5', char(bytes')));
end

function out = md5_file_system_(path)
    out = '';
    if ~isunix()
        return
    end

    quoted = shell_quote_(path);
    out = run_md5_command_(sprintf('md5sum %s 2>/dev/null', quoted), true);
    if ~isempty(out)
        return
    end
    out = run_md5_command_(sprintf('md5 -q %s 2>/dev/null', quoted), false);
end

function out = run_md5_command_(cmd, first_field)
    out = '';
    [status, text] = system(cmd);
    if status ~= 0
        return
    end
    text = strtrim(text);
    if first_field
        parts = regexp(text, '\s+', 'split');
        if isempty(parts)
            return
        end
        candidate = parts{1};
    else
        lines = regexp(text, '\r?\n', 'split');
        candidate = strtrim(lines{1});
    end
    candidate = lower(candidate);
    if ~isempty(regexp(candidate, '^[0-9a-f]{32}$', 'once'))
        out = candidate;
    end
end

function out = shell_quote_(path)
    replacement = ['''' '\' '''' ''''];
    out = ['''' strrep(path, '''', replacement) ''''];
end
