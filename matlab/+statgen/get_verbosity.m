function level = get_verbosity()
%GET_VERBOSITY Return statgen runtime verbosity.
%
%   level = statgen.get_verbosity()
%
% The returned level is one of 'quiet' or 'info'.
%
% See also statgen.set_verbosity.
    level = statgen.internal.verbosity_state('get', '');
end
