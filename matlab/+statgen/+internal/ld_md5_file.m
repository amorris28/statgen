function out = ld_md5_file(path)
    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
    clear cleaner;
    out = md5_hex_(bytes);
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
