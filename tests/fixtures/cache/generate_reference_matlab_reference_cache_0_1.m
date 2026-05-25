% Generate tests/fixtures/cache/reference_matlab_reference_cache_0_1.mat.
% Run from the repository root:
%   octave --no-gui --quiet tests/fixtures/cache/generate_reference_matlab_reference_cache_0_1.m

addpath(fileparts(mfilename('fullpath')));
cache_fixture_helpers('reference');
