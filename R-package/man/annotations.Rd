\name{annotations}
\alias{load_annotations}
\alias{create_annotations}
\alias{create_annotation}
\alias{save_annotations_cache}
\alias{load_annotations_cache}
\alias{annomat}
\alias{annonames}
\alias{num_annot}
\alias{select_annotations}
\alias{union_annotations}
\title{AnnotationPanel objects}
\description{
Paint BED interval annotations onto a reference panel, inspect sparse
annotation matrices, and save or load R-native annotation caches.
}
\usage{
load_annotations(bed_paths, reference)
create_annotations(reference, annotation_matrix, annotation_names)
create_annotation(reference, annovec, annoname)
save_annotations_cache(panel, path)
load_annotations_cache(path, shards = NULL)

annomat(x, ...)
annonames(x, ...)
num_annot(x, ...)
select_annotations(x, names, ...)
union_annotations(x, other, mode = "by_name", ...)
}
\arguments{
\item{bed_paths}{Character path or character vector of BED files.}
\item{reference}{A \code{ReferencePanel}.}
\item{annotation_matrix}{A binary matrix with SNPs as rows and annotations as columns.}
\item{annotation_names, names}{Character annotation names.}
\item{annovec}{Binary annotation vector aligned to \code{reference}.}
\item{annoname}{Character scalar annotation name.}
\item{panel, x, other}{Annotation panel objects.}
\item{path}{Character scalar cache path.}
\item{shards}{Optional character vector of shard labels in canonical order.}
\item{mode}{Union mode. Currently only \code{"by_name"} is supported.}
\item{...}{Reserved for S3 methods.}
}
\value{
\code{load_annotations()}, \code{create_annotations()},
\code{create_annotation()}, and \code{load_annotations_cache()} return
\code{AnnotationPanel} S3 objects. \code{annomat(x)} returns a sparse
\code{Matrix::lgCMatrix} with column names equal to \code{annonames(x)}.
}
\details{
BED intervals are interpreted as 0-based half-open intervals and are matched to
reference chromosome labels exactly. The annotation matrix is binary and
reference-aligned.
}
\examples{
ref_path <- system.file("extdata", "reference_chr1.bim", package = "statgen")
bed_path <- system.file("extdata", "anno1.bed", package = "statgen")
if (nzchar(ref_path) && nzchar(bed_path)) {
  ref <- load_reference(ref_path)
  ann <- load_annotations(bed_path, ref)
  dim(annomat(ann))
}
}
