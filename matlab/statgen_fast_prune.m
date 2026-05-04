function varargout = statgen_fast_prune(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.fast_prune(varargin{:});
end
