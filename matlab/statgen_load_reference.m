function varargout = statgen_load_reference(varargin)
%STATGEN_LOAD_REFERENCE Compatibility wrapper for statgen.load_reference.
%
%   See help statgen.load_reference.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_reference(varargin{:});
end
