% Generate tests/fixtures/cache/genotype_matlab_genotype_cache_0_1.mat.
% Run from the repository root:
%   octave --no-gui --quiet tests/fixtures/cache/generate_genotype_matlab_genotype_cache_0_1.m

addpath(fileparts(mfilename('fullpath')));
cache_fixture_helpers('genotype');
