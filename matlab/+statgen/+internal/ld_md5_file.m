function out = ld_md5_file(path)
    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*uint8');
    clear cleaner;
    out = lower(hash('md5', char(bytes')));
end
