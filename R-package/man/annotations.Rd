\name{annotations}
\alias{load_annotations}
\alias{load_annotation}
\alias{create_annotations}
\alias{create_annotation}
\alias{save_annotations_cache}
\alias{load_annotations_cache}
\alias{annomat}
\alias{annonames}
\alias{num_annot}
\alias{is_binary}
\alias{annotation_metadata}
\alias{select_annotations}
\alias{union_annotations}
\title{AnnotationPanel objects}
\description{
Paint binary BED and continuous BED-like interval annotations onto a reference
panel, inspect sparse annotation matrices, and save or load R-native annotation
caches.
}
\usage{
load_annotations(bed_paths, reference, annotation_metadata = NULL,
  annotation_metadata_paths = NULL)
load_annotation(path, reference, header = FALSE, value_columns = NULL,
  annotation_names = NULL, annotation_metadata = NULL,
  annotation_metadata_path = NULL)
create_annotations(reference, annotation_matrix, annotation_names,
  is_binary = NULL, annotation_metadata = NULL)
create_annotation(reference, annovec, annotation_name,
  is_binary = NULL, annotation_metadata = NULL)
save_annotations_cache(panel, path)
load_annotations_cache(path, shards = NULL)

annomat(x, ...)
annonames(x, ...)
num_annot(x, ...)
is_binary(x, ...)
annotation_metadata(x, ...)
select_annotations(x, names, ...)
union_annotations(x, other, mode = "by_name", ...)
}
\arguments{
\item{bed_paths}{Character path or character vector of 3-column BED files.}
\item{path}{Character scalar annotation source or cache path.}
\item{reference}{A \code{ReferencePanel}.}
\item{header}{Logical scalar; whether \code{path} has a header row.}
\item{value_columns}{Optional value columns for \code{load_annotation()}; R uses
one-based physical column indices or header names.}
\item{annotation_matrix}{A finite numeric matrix with SNPs as rows and annotations as columns.}
\item{annotation_names, names}{Character annotation names.}
\item{annovec}{Finite numeric annotation vector aligned to \code{reference}.}
\item{annotation_name}{Character scalar annotation name.}
\item{is_binary}{Optional logical vector marking columns that must contain only
\code{0} and \code{1}; inferred when omitted for \code{create_annotations()}.}
\item{annotation_metadata}{Optional character vector with one metadata string per
annotation column.}
\item{annotation_metadata_paths}{Optional character vector of sidecar paths for
\code{load_annotations()}; use \code{NA_character_} or \code{""} for generated
metadata.}
\item{annotation_metadata_path}{Optional sidecar path for \code{load_annotation()}.}
\item{panel, x, other}{Annotation panel objects.}
\item{shards}{Optional character vector of shard labels in canonical order.}
\item{mode}{Union mode. Currently only \code{"by_name"} is supported.}
\item{...}{Reserved for S3 methods.}
}
\value{
\code{load_annotations()}, \code{load_annotation()}, \code{create_annotations()},
\code{create_annotation()}, and \code{load_annotations_cache()} return
\code{AnnotationPanel} S3 objects. \code{annomat(x)} returns a sparse numeric
\code{Matrix::dgCMatrix} with column names equal to \code{annonames(x)}.
}
\details{
BED intervals are interpreted as 0-based half-open intervals and are matched to
reference chromosome labels exactly. \code{load_annotations()} treats each BED
file as one binary annotation. \code{load_annotation()} loads one source file
that may be a 3-column binary BED file or a BED-like file with numeric value
columns.

For \code{load_annotation()}, a 3-column input with \code{value_columns = NULL}
produces one binary annotation. A 4-column input with
\code{value_columns = NULL} uses column 4 as one continuous annotation.
Inputs with five or more columns require explicit \code{value_columns}; for
headerless files they also require \code{annotation_names}. Integer
\code{value_columns} are one-based physical file columns and must select
columns 4 or later. Named selectors require \code{header = TRUE}.

\code{is_binary} marks columns that should be treated as binary membership
annotations. When omitted in \code{create_annotations()}, it is inferred by
checking each sparse matrix column: a column is binary when all stored non-zero
values are exactly \code{1}; implicit zeros are allowed. Declared binary
columns must contain only \code{0} and \code{1}.

\code{annotation_metadata} contains opaque per-column strings. When metadata is
omitted, loaders generate compact JSON provenance strings containing the source
file and selected source column. For \code{load_annotations()},
\code{annotation_metadata_paths} entries are file-level sidecars whose full
contents are preserved. For column-valued \code{load_annotation()} inputs,
\code{annotation_metadata_path} must contain one non-empty physical line per
physical input column; metadata for selected value columns is taken from the
matching sidecar lines.
}
\examples{
ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen")
bed_path <- system.file("extdata", "anno1.bed", package = "statgen")
ref <- load_reference(ref_path)
ann <- load_annotations(bed_path, ref)
dim(annomat(ann))
is_binary(ann)

continuous_path <- tempfile(fileext = ".annot")
writeLines(c(
  "chrom\tstart0\tend0\tscore",
  "1\t99\t200\t0.5",
  "1\t299\t400\t1.5"
), continuous_path)
continuous <- load_annotation(
  continuous_path,
  ref,
  header = TRUE,
  value_columns = "score"
)
is_binary(continuous)
annotation_metadata(continuous)

cache <- tempfile(fileext = ".rds")
save_annotations_cache(ann, cache)
ann_cached <- load_annotations_cache(cache)
dim(annomat(ann_cached))
}
