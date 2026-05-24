function out = annotation_generated_metadata(path, source_column0, source_column_name)
% Generate stable annotation provenance metadata.
    if isempty(source_column0)
        col = 'null';
    else
        col = sprintf('%d', source_column0);
    end
    if isempty(source_column_name)
        col_name = 'null';
    else
        col_name = ['"' json_escape_(source_column_name) '"'];
    end
    out = sprintf('{"source_column0":%s,"source_column_name":%s,"source_file":"%s"}', ...
        col, col_name, json_escape_(path));
end

function out = json_escape_(value)
    s = char(value);
    s = strrep(s, '\', '\\');
    s = strrep(s, '"', '\"');
    s = strrep(s, sprintf('\n'), '\n');
    s = strrep(s, sprintf('\r'), '\r');
    s = strrep(s, sprintf('\t'), '\t');
    out = s;
end
