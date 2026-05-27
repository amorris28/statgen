function out = annotation_generated_metadata(path, has_header, num_source_intervals, varargin)
% Generate stable annotation provenance metadata.
    fields = {
        'source_file', ['"' json_escape_(path) '"']
        'source_file_has_header', sprintf('%d', logical(has_header))
        'num_source_intervals', sprintf('%d', num_source_intervals)
    };
    if mod(numel(varargin), 2) ~= 0
        error('statgen:annotations', 'annotation_generated_metadata options must be name-value pairs');
    end
    for i = 1:2:numel(varargin)
        name = char(varargin{i});
        value = varargin{i + 1};
        if isempty(value)
            continue
        end
        fields(end + 1, :) = {name, json_value_(value)}; %#ok<AGROW>
    end
    parts = cell(size(fields, 1), 1);
    for i = 1:size(fields, 1)
        parts{i} = sprintf('"%s":%s', fields{i, 1}, fields{i, 2});
    end
    out = ['{' strjoin(parts, ',') '}'];
end

function out = json_value_(value)
    if isnumeric(value) || islogical(value)
        if ~isscalar(value)
            error('statgen:annotations', 'metadata values must be scalar');
        end
        value = double(value);
        if floor(value) == value
            out = sprintf('%d', int64(value));
        else
            out = sprintf('%.17g', value);
        end
    else
        out = ['"' json_escape_(value) '"'];
    end
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
