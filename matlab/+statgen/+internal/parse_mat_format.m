function [format, save_arg] = parse_mat_format(where, default_format, allowed_formats, varargin)
% Parse shared MATLAB .mat format name-value options.
    if nargin < 3 || isempty(allowed_formats)
        allowed_formats = {'v7', 'v7.3', 'v5'};
    end
    format = normalize_format_(default_format, where, allowed_formats);

    if mod(numel(varargin), 2) ~= 0
        error('statgen:io', '%s options must be name-value pairs', where);
    end
    seen_format = false;
    for i = 1:2:numel(varargin)
        if ~(ischar(varargin{i}) || isstring(varargin{i}))
            error('statgen:io', '%s option names must be strings', where);
        end
        name = char(varargin{i});
        if strcmpi(name, 'format')
            if seen_format
                error('statgen:io', '%s format option specified more than once', where);
            end
            seen_format = true;
            format = normalize_format_(varargin{i + 1}, where, allowed_formats);
        else
            error('statgen:io', '%s unknown option: %s', where, name);
        end
    end
    save_arg = format_save_arg_(format);
end

function format = normalize_format_(value, where, allowed_formats)
    if ~(ischar(value) || isstring(value)) || isempty(value)
        error('statgen:io', '%s format must be one of: %s', where, strjoin(allowed_formats, ', '));
    end
    format = lower(char(value));
    if ~any(strcmp(format, allowed_formats))
        error('statgen:io', '%s format must be one of: %s', where, strjoin(allowed_formats, ', '));
    end
end

function save_arg = format_save_arg_(format)
    switch format
        case 'v7.3'
            save_arg = '-v7.3';
        case 'v7'
            save_arg = '-v7';
        case 'v5'
            save_arg = '-mat';
    end
end
