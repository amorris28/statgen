% Generate tests/fixtures/cache/sumstats_matlab_sumstats_cache_0_1.mat.
% Run from the repository root:
%   octave --no-gui --quiet tests/fixtures/cache/generate_sumstats_matlab_sumstats_cache_0_1.m

addpath(fileparts(mfilename('fullpath')));
cache_fixture_helpers('sumstats');
