function varargout = statgen_load_annotation(varargin)
%STATGEN_LOAD_ANNOTATION Compatibility wrapper for statgen.load_annotation.
%
%   See help statgen.load_annotation.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_annotation(varargin{:});
end
