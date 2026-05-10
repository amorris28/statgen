function out = display_size_string(sz)
%DISPLAY_SIZE_STRING Format an array size vector for object display.
    parts = cell(1, numel(sz));
    for i = 1:numel(sz)
        parts{i} = sprintf('%d', sz(i));
    end
    out = strjoin(parts, 'x');
end
