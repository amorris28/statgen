function varargout = statgen_version(varargin)
%STATGEN_VERSION Compatibility wrapper for statgen.version.
%
%   See help statgen.version.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.version(varargin{:});
end
