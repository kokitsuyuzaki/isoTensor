\name{isoTensor}
\docType{methods}
\alias{isoTensor}

\title{
Cell-Type-Resolved Isoform Abundance Estimation
}
\description{
isoTensor estimates cell-type-resolved isoform abundance from a bulk
isoform abundance matrix X (individual x isoform) and an externally
estimated cell type fraction matrix A (individual x cell type).
A is treated as fixed and is never updated; only the
cell-type-by-isoform matrix B is estimated by multiplicative updates
so that X is approximated by A B.
}
\usage{
isoTensor(X, A, model=c("matrix", "tensor"),
  rank_individual=NULL, rank_isoform=NULL,
  initB=NULL, pseudocount=.Machine$double.eps,
  thr=1e-10, num.iter=100, verbose=FALSE)
}

\arguments{
\item{X}{
The input matrix which has N-rows and I-columns, where N is the number
of individuals and I is the number of isoforms. Each element must be
non-negative (e.g., TPM or estimated counts of isoforms).
}
\item{A}{
The cell type fraction matrix which has N-rows and C-columns, where C
is the number of cell types. A is estimated in advance by external
reference-based deconvolution tools (e.g., Bisque) and is fixed
throughout the optimization. Each element must be non-negative.
}
\item{model}{
"matrix" performs the population-level isoform deconvolution
(X ~ A B). "tensor" (individual-specific cell-type-resolved isoform
abundance) is not implemented yet (Default: "matrix").
}
\item{rank_individual}{
Reserved for model="tensor". The rank of the individual mode
(Default: NULL).
}
\item{rank_isoform}{
Reserved for model="tensor". The rank of the isoform mode
(Default: NULL).
}
\item{initB}{
The initial values of the cell-type-by-isoform matrix B, which has
C-rows and I-columns (Default: NULL, i.e., random initialization).
}
\item{pseudocount}{
The pseudo count added to the denominator of the multiplicative update
to avoid division by zero (Default: Machine Epsilon).
}
\item{thr}{
When the relative change of the reconstruction error is lower than
thr, the iteration is terminated (Default: 1E-10).
}
\item{num.iter}{
The number of iterations (Default: 100).
}
\item{verbose}{
If verbose == TRUE, the relative change of the reconstruction error is
shown in each iteration (Default: FALSE).
}
}

\value{
  B : A matrix which has C-rows and I-columns, the population-level
  cell-type-resolved isoform abundance.
  model : The model actually used ("matrix").
  algorithm : The objective function ("Frobenius").
  RecError : The history of the reconstruction error
  ||X - AB||_F (0-th iteration is named "offset").
  RelChange : The history of the relative change of the reconstruction
  error (0-th iteration is named "offset").
}

\examples{
toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
set.seed(1234)
fit <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000, thr=1e-12)
max(abs(fit$B - toy$B))
}

\author{Koki Tsuyuzaki}

\seealso{\code{\link{isoToyModel}}}
\keyword{methods}
