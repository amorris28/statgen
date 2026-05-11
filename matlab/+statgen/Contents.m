%STATGEN Statistical genetics data structures for MATLAB and Octave.
%
% statgen provides reference-aligned objects for SNP metadata, LD matrices,
% GWAS summary statistics, SNP annotations, and PLINK-backed genotype access.
%
% Loaders:
%   load_reference          Load a ReferencePanel from BIM metadata.
%   load_ld                 Load an LDPanel from a MATLAB LD distribution.
%   load_sumstats           Load summary statistics aligned to a reference.
%   load_annotations        Load BED annotations aligned to a reference.
%   load_genotype           Load PLINK genotype metadata aligned to a reference.
%
% Cache loaders:
%   load_reference_cache    Load a cached ReferencePanel.
%   load_sumstats_cache     Load cached summary statistics.
%   load_annotations_cache  Load cached annotations.
%   load_genotype_cache     Load cached genotype metadata.
%
% Cache writers:
%   save_reference_cache    Save a ReferencePanel cache.
%   save_sumstats_cache     Save a Sumstats cache.
%   save_annotations_cache  Save an AnnotationPanel cache.
%   save_genotype_cache     Save a GenotypePanel metadata cache.
%
% Core operations:
%   fast_prune              LD-prune a score vector.
%   create_sumstats         Create Sumstats from aligned vectors.
%   create_annotation       Create a one-column AnnotationPanel.
%   create_annotations      Create annotations from an aligned matrix.
%   validate_ld_distribution Validate a MATLAB LD distribution.
%   convert_ld_npz_to_mat   Convert Python LD shards to MATLAB .mat shards.
%   create_ld_mat_manifest  Create a manifest for converted MATLAB LD shards.
%   set_verbosity          Set runtime verbosity: quiet or info.
%   get_verbosity          Return current runtime verbosity.
%
% LD preparation:
%   Build Python .npz LD shards with script/statgen_build_ld.py, convert them
%   with statgen.convert_ld_npz_to_mat, then create the MATLAB manifest with
%   statgen.create_ld_mat_manifest. See docs/TUTORIAL_1_PREPARE_DATA.md.
%
% Public objects:
%   ReferencePanel          Reference SNP coordinate system.
%   LDPanel                 Sparse LD panel aligned to a ReferencePanel.
%   Sumstats                Aligned summary statistics for one trait.
%   AnnotationPanel         Aligned annotation matrix.
%   GenotypePanel           Reference-aligned PLINK genotype metadata.
%
% Use help statgen.<name> for function and object-specific help.
