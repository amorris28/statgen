# Docstrings and help text

## MATLAB/Octave

Public MATLAB/Octave entry points must provide concise `help` text for
happy-path use. The help should tell users what the function or object is for,
show the common call forms, mention the main outputs or properties, and point
to nearby public APIs with `See also`.

Document user-facing package functions and these public object classes:

- `ReferencePanel`
- `LDPanel`
- `Sumstats`
- `AnnotationPanel`
- `GenotypePanel`

Do not write user-facing help for shard classes or `+statgen/+internal/`
helpers. Users are not expected to construct or operate on shards directly.

The package should provide a short overview in `matlab/+statgen/Contents.m`.
MATLAB uses this file for `help statgen`. Octave does not resolve bare package
help the same way, but the same text is available with `help statgen.Contents`.
Do not add a top-level `matlab/statgen.m` help shim, because it shadows the
`statgen` package and can break calls such as `statgen.load_reference(...)`.

`Contents.m` should briefly state what `statgen` is, list the main public
loaders/cache loaders and core operations, and point users to the public panel
objects for object-specific help. It should remain an overview, not a tutorial
or a duplicate of every function docstring.

Class help belongs in the first contiguous comment block after `classdef`.
Public method and package-function help belongs in the first contiguous comment
block after the `function` line. Keep implementation comments such as `%#ok`
after the help block so both MATLAB and Octave can find the intended text.

Docstrings should be short enough to scan. Prefer one or two common call forms,
the few properties or methods users normally inspect, and any array orientation
that is easy to get wrong. Do not copy full spec text, list every field, explain
cache internals, or catalog every warning and error.

## Python

Python docstrings are out of scope for now. Python documentation requirements
may be added later if the Python API needs interactive help conventions beyond
the existing specs and tests.
