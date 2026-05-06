function out = ld_md5_file(path)
    path = char(path);
    if exist(path, 'file') ~= 2
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    out = md5_file_system_(path);
    if ~isempty(out)
        return
    end

    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
    clear cleaner;
    out = md5_hex_(bytes);
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

function out = md5_hex_(bytes)
    try
        out = lower(hash('md5', char(bytes')));
        return
    catch
        % MATLAB path: use Java MessageDigest when hash(...) is unavailable.
    end

    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(bytes(:)));
    d = typecast(md.digest(), 'uint8');
    out = lower(reshape(dec2hex(d)', 1, []));
end
