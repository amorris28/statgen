% Generate tests/fixtures/cache/annotations_matlab_annotations_cache_0_1.mat.
%
% Run from a v0.3.2 checkout at the repository root:
%   octave --no-gui --quiet tests/fixtures/cache/generate_annotations_matlab_annotations_cache_0_1.m
%
% This uses the official v0.3.2 MATLAB/Octave implementation, whose annotation
% writer produces annotations_cache/0.1.

cache_dir = fileparts(mfilename('fullpath'));
repo_root = pwd;
fixture_dir = fullfile('tests', 'fixtures');
addpath(fullfile(repo_root, 'matlab'));

ref = statgen.load_reference(fullfile(fixture_dir, 'reference', 'sharded', '@.bim'));
annotations = statgen.load_annotations( ...
    {fullfile(fixture_dir, 'annotations', 'anno1.bed'), ...
     fullfile(fixture_dir, 'annotations', 'anno2.bed')}, ...
    ref);

if exist(cache_dir, 'dir') ~= 7
    mkdir(cache_dir);
end
statgen.save_annotations_cache( ...
    annotations, ...
    fullfile(cache_dir, 'annotations_matlab_annotations_cache_0_1.mat'), ...
    'format', 'v7');
