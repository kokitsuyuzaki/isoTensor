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
A is treated as fixed and is never updated.

With model="matrix" (Phase 1), the population-level
cell-type-by-isoform matrix B is estimated by multiplicative updates
so that X is approximated by A B.

With model="tensor" (Phase 2), an individual-specific latent tensor
B (individual x cell type x isoform) with a plain partial non-negative
Tucker structure B_nci = sum_pq U_np G_pcq W_iq (the cell type mode is
not compressed) is estimated by block multiplicative updates so that
X_ni is approximated by sum_c A_nc B_nci. Each block update is a
Lee-Seung non-negative least squares step, so the objective
||X - Xhat||_F is non-increasing (see notes/tensor-MU.md in the
GitHub repository for the derivation).
}
\usage{
isoTensor(X, A, model=c("matrix", "tensor"),
  rank_individual=NULL, rank_isoform=NULL,
  initB=NULL, initU=NULL, initG=NULL, initW=NULL,
  pseudocount=.Machine$double.eps,
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
(X ~ A B). "tensor" estimates the individual-specific latent tensor
under the partial Tucker structure (Default: "matrix").
}
\item{rank_individual}{
model="tensor" only (required there). The rank R_N of the individual
mode (Default: NULL).
}
\item{rank_isoform}{
model="tensor" only (required there). The rank R_I of the isoform
mode (Default: NULL).
}
\item{initB}{
model="matrix" only. The initial values of the cell-type-by-isoform
matrix B, which has C-rows and I-columns (Default: NULL, i.e., random
initialization).
}
\item{initU}{
model="tensor" only. Initial N x rank_individual non-negative matrix
(Default: NULL, i.e., strictly positive random initialization).
Zeros in user-supplied initial factors stay locked at zero under
multiplicative updates.
}
\item{initG}{
model="tensor" only. Initial rank_individual x C x rank_isoform
non-negative array (Default: NULL).
}
\item{initW}{
model="tensor" only. Initial I x rank_isoform non-negative matrix
(Default: NULL).
}
\item{pseudocount}{
The pseudo count added to the denominators of the multiplicative
updates to avoid division by zero (Default: Machine Epsilon).
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

\section{Warnings}{
For model="tensor", two theoretically flagged regimes are allowed but
produce warnings (they are not input errors; see notes/tensor-model.md
in the GitHub repository):
rank_individual >= rank_isoform, where the individual-specific latent
tensor B is generically non-identifiable even when X is reconstructed
well (Theorem A/A'); and a number of individuals below the
necessary-condition threshold N* of Theorem C.
}

\value{
For model="matrix":
  B : A matrix which has C-rows and I-columns, the population-level
  cell-type-resolved isoform abundance.
  model, algorithm : "matrix", "Frobenius".
  RecError : The history of the reconstruction error ||X - AB||_F
  (0-th iteration is named "offset").
  RelChange : The history of the relative change of the reconstruction
  error.

For model="tensor":
  U : N x rank_individual factor matrix (individual patterns).
  G : rank_individual x C x rank_isoform core tensor.
  W : I x rank_isoform factor matrix (isoform programs).
  B : The reconstructed latent tensor (N x C x I),
  B_nci = sum_pq U_np G_pcq W_iq.
  Xhat : The reconstructed bulk matrix (N x I),
  Xhat_ni = sum_c A_nc B_nci.
  model, algorithm : "tensor", "Frobenius".
  rank_individual, rank_isoform : The specified ranks.
  RecError, RelChange : Histories as above (for ||X - Xhat||_F).
  NumIter : The number of performed iterations.
  Converged : TRUE if the relative change reached thr within num.iter
  iterations.

Columns of U and W are normalized to unit Euclidean norm with the
scales absorbed into G (numerical stabilization; B and Xhat are
unaffected).
}

\examples{
# Phase 1: matrix model
toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
set.seed(1234)
fit <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000, thr=1e-12)
max(abs(fit$B - toy$B))

# Phase 2: tensor model (identifiable candidate regime R_N < R_I)
tt <- isoToyModel(model="tensor", N=20, C=2, I=20,
    rank_individual=1, rank_isoform=2, A_variation="high", seed=10)
set.seed(10)
fit2 <- isoTensor(tt$X, tt$A, model="tensor",
    rank_individual=1, rank_isoform=2, num.iter=100)
sqrt(sum((fit2$B - tt$B)^2)) / sqrt(sum(tt$B^2))
}

\author{Koki Tsuyuzaki}

\seealso{\code{\link{isoToyModel}}}
\keyword{methods}
