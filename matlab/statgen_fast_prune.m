function varargout = statgen_fast_prune(varargin)
%STATGEN_FAST_PRUNE Compatibility wrapper for statgen.fast_prune.
%
%   See help statgen.fast_prune.
    varargout = cell(1, max(nargout, 1));
    [varargout{:}] = statgen.fast_prune(varargin{:});
end
