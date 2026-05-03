function tf = ld_is_mat_v5(path)
    fid = fopen(path, 'rb');
    if fid < 0
        error('statgen:io', 'Cannot open LD file: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    header = fread(fid, 116, '*char')';
    clear cleaner;
    tf = strncmp(header, 'MATLAB 5.0 MAT-file', 19);
end
