function varargout = statgen_load_ld_reference(varargin)
%STATGEN_LOAD_LD_REFERENCE Compatibility wrapper for statgen.load_ld_reference.
%
%   See help statgen.load_ld_reference.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_ld_reference(varargin{:});
end
