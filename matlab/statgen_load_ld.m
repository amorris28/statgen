function varargout = statgen_load_ld(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_ld(varargin{:});
end
