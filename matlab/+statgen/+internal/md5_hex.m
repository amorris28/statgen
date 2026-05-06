function out = md5_hex(text_payload)
% Return lowercase MD5 hex digest for a char/uint8 payload.
    try
        out = lower(hash('md5', text_payload));
        return
    catch
        % MATLAB path: use Java MessageDigest when hash(...) is unavailable.
    end

    md = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(text_payload));
    d = typecast(md.digest(), 'uint8');
    out = lower(reshape(dec2hex(d)', 1, []));
end
