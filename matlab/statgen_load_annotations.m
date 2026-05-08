function varargout = statgen_load_annotations(varargin)
%STATGEN_LOAD_ANNOTATIONS Compatibility wrapper for statgen.load_annotations.
%
%   See help statgen.load_annotations.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_annotations(varargin{:});
end
