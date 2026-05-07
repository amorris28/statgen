function varargout = statgen_load_genotype(varargin)
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_genotype(varargin{:});
end
