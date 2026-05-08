function varargout = statgen_load_genotype(varargin)
%STATGEN_LOAD_GENOTYPE Compatibility wrapper for statgen.load_genotype.
%
%   See help statgen.load_genotype.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_genotype(varargin{:});
end
