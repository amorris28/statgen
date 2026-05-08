function varargout = statgen_load_sumstats(varargin)
%STATGEN_LOAD_SUMSTATS Compatibility wrapper for statgen.load_sumstats.
%
%   See help statgen.load_sumstats.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.load_sumstats(varargin{:});
end
