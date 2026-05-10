function out = display_join_strings(values, max_items)
%DISPLAY_JOIN_STRINGS Format a short, truncated cell/string list.
    if nargin < 2 || isempty(max_items)
        max_items = 8;
    end
    if ischar(values) || isstring(values)
        values = cellstr(values(:));
    end
    values = statgen.internal.ensure_cell_col(values);
    if isempty(values)
        out = '{}';
        return;
    end
    n = min(numel(values), max_items);
    shown = values(1:n);
    if numel(values) > max_items
        shown{end + 1} = '...';
    end
    out = ['{' strjoin(shown(:)', ', ') '}'];
end
