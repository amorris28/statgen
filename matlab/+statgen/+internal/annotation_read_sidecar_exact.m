function out = annotation_read_sidecar_exact(path)
% Read an annotation metadata sidecar as exact text.
    path = char(path);
    fid = fopen(path, 'r');
    if fid < 0
        error('statgen:io', 'Cannot open annotation metadata sidecar: %s', path);
    end
    cleaner = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, '*char')';
    out = char(bytes);
end
