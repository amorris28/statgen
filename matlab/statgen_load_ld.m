function varargout = statgen_load_ld(varargin)
%STATGEN_LOAD_LD Compatibility wrapper for statgen.load_ld.
%
%   See help statgen.load_ld.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_ld(varargin{:});
end
