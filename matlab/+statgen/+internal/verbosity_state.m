function out = verbosity_state(action, arg)
%VERBOSITY_STATE Store package runtime verbosity for the current session.
    persistent current_level
    if isempty(current_level)
        current_level = 'info';
    end

    switch char(action)
        case 'get'
            out = current_level;
        case 'set'
            validate_level_(arg);
            current_level = char(arg);
            out = current_level;
        case 'info'
            if strcmp(current_level, 'info')
                fprintf(2, '%s\n', char(arg));
            end
            out = current_level;
        otherwise
            error('statgen:verbosity', 'Unknown verbosity action: %s', char(action));
    end
end

function validate_level_(level)
    if ~(ischar(level) || (isstring(level) && isscalar(level)))
        error('statgen:verbosity', 'verbosity must be one of: quiet, info');
    end
    if ~any(strcmp(char(level), {'quiet', 'info'}))
        error('statgen:verbosity', 'verbosity must be one of: quiet, info');
    end
end
