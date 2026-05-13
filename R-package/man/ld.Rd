\name{ld}
\alias{load_ld}
\alias{load_ld_reference}
\alias{validate_ld_distribution}
\alias{convert_ld_npz_to_rds}
\alias{create_ld_rds_manifest}
\alias{a1freq}
\alias{multiply_r2}
\alias{fast_prune}
\alias{reference}
\alias{default_chrX_sex}
\title{LD panel loading, conversion, and operations}
\usage{
load_ld(path, shards = NULL, default_chrX_sex = NULL, retain_ld_r = TRUE)
load_ld_reference(path, shards = NULL)
validate_ld_distribution(path, check_payload_structure = FALSE)
convert_ld_npz_to_rds(npz_root, rds_root, shard)
create_ld_rds_manifest(npz_root, rds_root, shards)
a1freq(x, chrX_sex = NULL, ...)
multiply_r2(ld_panel, M, chrX_sex = NULL, ...)
fast_prune(logpvec, ld_panel, r2_threshold = 0.2, chrX_sex = NULL)
reference(x, ...)
default_chrX_sex(x, ...)
}
\arguments{
\item{path}{Path to an R LD distribution root directory.}
\item{shards}{Optional character vector of canonical shard labels to load or finalize.}
\item{default_chrX_sex}{Default chrX LD shard selector, one of \code{"female"}, \code{"male"}, or \code{"combined"}.}
\item{retain_ld_r}{Logical; retain raw signed LD \code{r} matrices in loaded shards. Operations keep working when this is \code{FALSE}.}
\item{check_payload_structure}{Logical; perform expensive sparse payload checks during distribution validation.}
\item{npz_root}{Path to a Python \code{.npz} LD handoff distribution root.}
\item{rds_root}{Path to the output R \code{.rds} LD distribution root.}
\item{shard}{Single canonical shard label to convert. For \code{"X"}, all chrX sex-label shards in the Python manifest are converted.}
\item{x}{An LD panel or shard object.}
\item{ld_panel}{An \code{LDPanel}.}
\item{M}{A numeric vector or matrix with rows aligned to the LD panel reference.}
\item{logpvec}{Numeric vector of signed or unsigned log-p values aligned to the LD panel reference.}
\item{r2_threshold}{Finite non-negative pruning threshold; defaults to \code{0.2}.}
\item{chrX_sex}{Optional chrX LD shard selector for operations. Ignored when chrX is absent.}
\item{...}{Reserved for S3 method dispatch.}
}
\description{
Load R-native LD distributions, validate their manifests and RDS sparse CSC
payloads, convert documented Python LD \code{.npz} handoff shards to RDS, and
run the core LD operations.
}
\details{
R LD distributions use \code{ld_manifest.json}, bundled reference BIM files, an
RDS reference cache, and one RDS file per LD shard. RDS shard files store signed
LD \code{r} in CSC components; \code{load_ld()} constructs
\code{Matrix::dgCMatrix} objects in memory and derives the elementwise squared
LD matrix used by \code{multiply_r2()} and \code{fast_prune()}.

\code{convert_ld_npz_to_rds()} is a limited handoff reader for the documented
statgen LD \code{.npz} format. It is not a general NumPy archive API.
}
\value{
\code{load_ld()} returns an \code{LDPanel}. \code{load_ld_reference()} returns
the paired \code{ReferencePanel}. \code{validate_ld_distribution()} returns a
list with \code{ok = TRUE} on success. Conversion helpers return invisibly.
\code{a1freq()} returns aligned allele frequencies, \code{multiply_r2()} returns
a vector or matrix with the same shape as \code{M}, and \code{fast_prune()}
returns a numeric vector or one-dimensional matrix with pruned entries set to
\code{NaN}.
}
\keyword{datasets}
