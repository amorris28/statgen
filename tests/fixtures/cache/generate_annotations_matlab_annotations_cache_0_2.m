% Generate tests/fixtures/cache/annotations_matlab_annotations_cache_0_2.mat.
% Run from the repository root:
%   octave --no-gui --quiet tests/fixtures/cache/generate_annotations_matlab_annotations_cache_0_2.m

addpath(fileparts(mfilename('fullpath')));
cache_fixture_helpers('annotations');
