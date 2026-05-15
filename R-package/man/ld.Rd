\name{ld}
\alias{load_ld}
\alias{load_ld_reference}
\alias{validate_ld_distribution}
\alias{prepare_ld_npz_for_r}
\alias{a1freq}
\alias{multiply_r2}
\alias{fast_prune}
\alias{reference}
\alias{default_chrX_sex}
\title{LD panel loading, preparation, and operations}
\usage{
load_ld(path, shards = NULL, default_chrX_sex = NULL, retain_ld_r = TRUE)
load_ld_reference(path, shards = NULL)
validate_ld_distribution(path, check_payload_structure = FALSE)
prepare_ld_npz_for_r(npz_root, extract_npz = FALSE)
a1freq(x, chrX_sex = NULL, ...)
multiply_r2(ld_panel, M, chrX_sex = NULL, ...)
fast_prune(logpvec, ld_panel, r2_threshold = 0.2, chrX_sex = NULL)
reference(x, ...)
default_chrX_sex(x, ...)
}
\arguments{
\item{path}{Path to a Python \code{.npz} LD distribution root prepared for R loading.}
\item{shards}{Optional character vector of canonical shard labels to load.}
\item{default_chrX_sex}{Default chrX LD shard selector, one of \code{"female"}, \code{"male"}, or \code{"combined"}.}
\item{retain_ld_r}{Logical; retain raw signed LD \code{r} matrices in loaded shards. Operations keep working when this is \code{FALSE}.}
\item{check_payload_structure}{Logical; perform expensive sparse payload checks during distribution validation.}
\item{npz_root}{Path to a Python \code{.npz} LD distribution root.}
\item{extract_npz}{Logical; when \code{TRUE}, create or refresh sibling
\code{*.npz.d} directories containing extracted \code{.npy} arrays for faster
R loading.}
\item{x}{An LD panel or shard object.}
\item{ld_panel}{An \code{LDPanel}.}
\item{M}{A numeric vector or matrix with rows aligned to the LD panel reference.}
\item{logpvec}{Numeric vector of signed or unsigned log-p values aligned to the LD panel reference.}
\item{r2_threshold}{Finite non-negative pruning threshold; defaults to \code{0.2}.}
\item{chrX_sex}{Optional chrX LD shard selector for operations. Ignored when chrX is absent.}
\item{...}{Reserved for S3 method dispatch.}
}
\description{
Prepare Python LD \code{.npz} distributions for R loading, validate their
manifests and sparse CSC payloads, load LD panels, and run the core LD
operations.
}
\details{
R reads Python LD \code{.npz} shard files directly. R-loadable distributions use
\code{ld_manifest.json}, bundled reference BIM files, Python \code{.npz} LD
shards, and an R reference cache sidecar recorded in manifest fields
\code{r_reference_cache} and \code{r_reference_cache_md5}. \code{load_ld()}
constructs \code{Matrix::dgCMatrix} objects in memory and derives the
elementwise squared LD matrix used by \code{multiply_r2()} and
\code{fast_prune()}.

\code{prepare_ld_npz_for_r()} builds or refreshes the R reference cache sidecar
from the bundled BIM files for all shards described in the existing Python
manifest and updates that manifest in place. With \code{extract_npz = TRUE}, it
also creates or refreshes extracted \code{*.npz.d} cache directories next to the
manifest-declared shard files. The operation is idempotent.
}
\value{
\code{load_ld()} returns an \code{LDPanel}. \code{load_ld_reference()} returns
the paired \code{ReferencePanel}. \code{validate_ld_distribution()} returns a
list with \code{ok = TRUE} on success. \code{prepare_ld_npz_for_r()} returns the
updated manifest invisibly.
\code{a1freq()} returns aligned allele frequencies, \code{multiply_r2()} returns
a vector or matrix with the same shape as \code{M}, and \code{fast_prune()}
returns a numeric vector or one-dimensional matrix with pruned entries set to
\code{NaN}.
}
\examples{
ld_src <- system.file("extdata", "ld", "python", package = "statgen")
if (nzchar(ld_src)) {
  ld_root <- tempfile("ld_npz_")
  dir.create(ld_root)
  file.copy(list.files(ld_src, full.names = TRUE), ld_root, recursive = TRUE)

  prepare_ld_npz_for_r(ld_root)
  ld <- load_ld(ld_root)

  num_snp(ld)
  head(a1freq(ld))

  scores <- seq_len(num_snp(ld))
  multiply_r2(ld, scores)

  fast_prune(c(9, 9, 7, 1, 2, 6, 5, 4), ld, r2_threshold = 0.35)
}
}
\keyword{datasets}
