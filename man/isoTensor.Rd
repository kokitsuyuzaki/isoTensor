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
not compressed) is estimated so that X_ni is approximated by
sum_c A_nc B_nci. Three optimizers are available for the same
Frobenius objective 1/2 ||X - Xhat||_F^2:

"lbfgs" (default): the non-negativity constraints are removed by the
log reparameterization U = exp(Theta_U) (and likewise for G, W) and
the unconstrained problem is solved by L-BFGS-B with the analytical
gradients (see notes/tensor-LBFGS.md in the GitHub repository). This
is the recommended optimizer: on synthetic benchmarks it reaches the
same or a better objective than the multiplicative updates in a
fraction of the time and recovers the latent tensor substantially
more accurately.

"mu": block multiplicative updates. Each block update is a Lee-Seung
non-negative least squares step, so the objective is provably
non-increasing (see notes/tensor-MU.md); kept as the theoretically
guaranteed monotone baseline.

"hybrid": the multiplicative updates first (controlled by num.iter and
thr), then L-BFGS refinement from the MU solution.
}
\usage{
isoTensor(X, A, model=c("matrix", "tensor"),
  rank_individual=NULL, rank_isoform=NULL,
  optimizer=c("lbfgs", "hybrid", "mu"),
  initB=NULL, initU=NULL, initG=NULL, initW=NULL,
  pseudocount=.Machine$double.eps,
  thr=1e-10, num.iter=100,
  optimizer.control=list(), verbose=FALSE)
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
\item{optimizer}{
model="tensor" only (ignored with a warning if explicitly supplied
with model="matrix", which is always fitted by multiplicative
updates). One of "lbfgs" (log-reparameterized L-BFGS, the default),
"hybrid" (multiplicative updates followed by L-BFGS refinement), or
"mu" (multiplicative updates only). The optimizer selects the
optimization method; the objective (algorithm = "Frobenius") is the
same for all three (Default: "lbfgs").
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
optimizer="mu" (multiplicative updates); under optimizer="lbfgs" and
in the L-BFGS stage of optimizer="hybrid" they are raised to 1e-10
before the log transform and can move away from zero.
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
updates to avoid division by zero. Used by the multiplicative updates
only (model="matrix", optimizer="mu", and the MU stage of
optimizer="hybrid") (Default: Machine Epsilon).
}
\item{thr}{
When the relative change of the reconstruction error is lower than
thr, the multiplicative-update iteration is terminated. Applies to
model="matrix", optimizer="mu", and the MU stage of
optimizer="hybrid"; thr=0 is allowed and simply runs all num.iter
iterations. Not used by the pure "lbfgs" path (Default: 1E-10).
}
\item{num.iter}{
The number of multiplicative-update iterations (same scope as thr)
(Default: 100).
}
\item{optimizer.control}{
model="tensor" with optimizer="lbfgs" or "hybrid" only (ignored with
a warning for optimizer="mu"). A named list of control values for the
L-BFGS stage: maxit (maximum number of L-BFGS iterations, Default:
5000) and factr (convergence tolerance of optim's L-BFGS-B in units
of machine epsilon, Default: 10). The defaults are the settings
validated on the synthetic benchmark grids (Default: list()).
}
\item{verbose}{
If verbose == TRUE, the relative change of the reconstruction error is
shown in each multiplicative-update iteration, and the L-BFGS stage
reports its convergence status (Default: FALSE).
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

For model="tensor" (all optimizers):
  U : N x rank_individual factor matrix (individual patterns).
  G : rank_individual x C x rank_isoform core tensor.
  W : I x rank_isoform factor matrix (isoform programs).
  B : The reconstructed latent tensor (N x C x I),
  B_nci = sum_pq U_np G_pcq W_iq.
  Xhat : The reconstructed bulk matrix (N x I),
  Xhat_ni = sum_c A_nc B_nci.
  model, algorithm, optimizer : "tensor", "Frobenius", and the used
  optimizer.
  rank_individual, rank_isoform : The specified ranks.
  objective : The final objective value 1/2 ||X - Xhat||_F^2.
  RecError : History of ||X - Xhat||_F. For optimizer="mu" one value
  per iteration ("offset" first); for "lbfgs" only the values at the
  start and end of the optimization (named "offset", "final"); for
  "hybrid" the MU history plus one "refined" value after L-BFGS.

Additionally for optimizer="mu" (and, describing the MU stage only,
for optimizer="hybrid"):
  RelChange : History of the relative change of the reconstruction
  error.
  NumIter : The number of performed MU iterations.
  Converged : TRUE if the relative change reached thr within num.iter
  iterations.
No per-iteration RelChange/NumIter is fabricated for the L-BFGS stage.

Additionally for optimizer="lbfgs" and "hybrid":
  lbfgs : A list with the optim() metadata of the L-BFGS stage:
  convergence (0 = converged, 1 = maxit reached), message, counts,
  maxit, factr, and failed (TRUE if optim() raised an error; for
  "hybrid" the MU-stage solution is then returned with a warning,
  while pure "lbfgs" stops with an error).

Columns of U and W are normalized to unit Euclidean norm with the
scales absorbed into G (numerical stabilization; B and Xhat are
unaffected). For the L-BFGS optimizers this normalization is applied
only once after the optimization finishes.
}

\examples{
# Phase 1: matrix model
toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
set.seed(1234)
fit <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000, thr=1e-12)
max(abs(fit$B - toy$B))

# Phase 2: tensor model (identifiable candidate regime R_N < R_I),
# fitted with the default log-reparameterized L-BFGS optimizer
tt <- isoToyModel(model="tensor", N=20, C=2, I=20,
    rank_individual=1, rank_isoform=2, A_variation="high", seed=10)
set.seed(10)
fit2 <- isoTensor(tt$X, tt$A, model="tensor",
    rank_individual=1, rank_isoform=2, optimizer="lbfgs")
sqrt(sum((fit2$B - tt$B)^2)) / sqrt(sum(tt$B^2))

# Multiplicative updates (monotone baseline) and MU -> L-BFGS hybrid
set.seed(10)
fit3 <- isoTensor(tt$X, tt$A, model="tensor",
    rank_individual=1, rank_isoform=2, optimizer="mu", num.iter=100)
set.seed(10)
fit4 <- isoTensor(tt$X, tt$A, model="tensor",
    rank_individual=1, rank_isoform=2, optimizer="hybrid", num.iter=100)
c(fit2$objective, fit3$objective, fit4$objective)
}

\author{Koki Tsuyuzaki}

\seealso{\code{\link{isoToyModel}}}
\keyword{methods}
