toy <- isoToyModel(model="tensor", N=20, C=2, I=20,
    rank_individual=1, rank_isoform=2,
    individual_effect=TRUE, A_variation="high", seed=10)

set.seed(10)
fit <- isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
    rank_individual=1, rank_isoform=2, num.iter=1500, thr=1e-12)

# Output structure
expect_true(is.list(fit))
expect_equivalent(names(fit),
    c("U", "G", "W", "B", "Xhat", "model", "algorithm", "optimizer",
      "rank_individual", "rank_isoform", "objective",
      "RecError", "RelChange", "NumIter", "Converged"))
expect_equivalent(fit$model, "tensor")
expect_equivalent(fit$optimizer, "mu")
expect_equal(fit$objective, 0.5 * tail(fit$RecError, 1)^2)
expect_equivalent(fit$rank_individual, 1)
expect_equivalent(fit$rank_isoform, 2)
expect_true(is.numeric(fit$NumIter))
expect_true(is.logical(fit$Converged))

# Dimensions
expect_equivalent(dim(fit$U), c(20, 1))
expect_equivalent(dim(fit$G), c(1, 2, 2))
expect_equivalent(dim(fit$W), c(20, 2))
expect_equivalent(dim(fit$B), c(20, 2, 20))
expect_equivalent(dim(fit$Xhat), c(20, 20))

# Non-negativity of all factors and reconstructions
expect_true(min(fit$U) >= 0)
expect_true(min(fit$G) >= 0)
expect_true(min(fit$W) >= 0)
expect_true(min(fit$B) >= 0)
expect_true(min(fit$Xhat) >= 0)

# Fixed seed reproducibility
set.seed(10)
fit2 <- isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
    rank_individual=1, rank_isoform=2, num.iter=1500, thr=1e-12)
expect_identical(fit$B, fit2$B)

# Full objective history never increases
expect_true(all(diff(fit$RecError[-1]) <= 1e-10))

# Positive control (identifiable candidate regime R_N < R_I):
# both X and the latent tensor B are recovered
relerr <- function(a, b) sqrt(sum((a - b)^2)) / sqrt(sum(b^2))
expect_true(relerr(fit$Xhat, toy$X) < 1e-3)
expect_true(relerr(fit$B, toy$B) < 0.05)

# A is unchanged after fitting
A_copy <- toy$A + 0
set.seed(10)
invisible(isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
    rank_individual=1, rank_isoform=2, num.iter=5))
expect_identical(toy$A, A_copy)

# thr=0 is a valid edge case: the sentinel stays positive, so exactly
# num.iter MU iterations are performed
set.seed(10)
fit_thr0 <- isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
    rank_individual=1, rank_isoform=2, num.iter=5, thr=0)
expect_equivalent(fit_thr0$NumIter, 5)
set.seed(1234)
toy_m0 <- isoToyModel(N=5, C=2, I=6, seed=1234)
fit_m0 <- isoTensor(toy_m0$X, toy_m0$A, model="matrix", num.iter=5, thr=0)
expect_equivalent(length(fit_m0$RecError), 6)

# Analytical gradient vs central finite differences (small config)
set.seed(123)
Ng <- 3; Cg <- 2; Ig <- 4; Rn <- 1; Ri <- 2
Ug <- matrix(runif(Ng*Rn, 0.5, 1.5), Ng, Rn)
Gg <- array(runif(Rn*Cg*Ri, 0.5, 1.5), dim=c(Rn, Cg, Ri))
Wg <- matrix(runif(Ig*Ri, 0.5, 1.5), Ig, Ri)
Ag <- matrix(runif(Ng*Cg, 0.2, 1), Ng, Cg); Ag <- Ag / rowSums(Ag)
Xg <- matrix(runif(Ng*Ig, 0, 2), Ng, Ig)
gr <- isoTensor:::.gradIsoTensor(Xg, Ag, Ug, Gg, Wg)
h <- 1e-6
fd_of <- function(theta, setfun){
    fd <- array(0, dim=dim(theta))
    for(j in seq_along(theta)){
        tp <- theta; tp[j] <- tp[j] + h
        tm <- theta; tm[j] <- tm[j] - h
        fd[j] <- (setfun(tp) - setfun(tm)) / (2 * h)
    }
    fd
}
fdU <- fd_of(Ug, function(t) isoTensor:::.objIsoTensor(Xg, Ag, t, Gg, Wg))
fdG <- fd_of(Gg, function(t) isoTensor:::.objIsoTensor(Xg, Ag, Ug, t, Wg))
fdW <- fd_of(Wg, function(t) isoTensor:::.objIsoTensor(Xg, Ag, Ug, Gg, t))
maxrel <- function(a, f) max(abs(a - f) / pmax(abs(a) + abs(f), 1e-8))
expect_true(maxrel(gr$U, fdU) < 1e-6)
expect_true(maxrel(gr$G, fdG) < 1e-6)
expect_true(maxrel(gr$W, fdW) < 1e-6)

# One-step objective non-increase for each block update separately
eps <- .Machine$double.eps
L0 <- isoTensor:::.objIsoTensor(Xg, Ag, Ug, Gg, Wg)
LU <- isoTensor:::.objIsoTensor(Xg, Ag,
    isoTensor:::.updateU_tensor(Xg, Ag, Ug, Gg, Wg, eps), Gg, Wg)
LG <- isoTensor:::.objIsoTensor(Xg, Ag, Ug,
    isoTensor:::.updateG_tensor(Xg, Ag, Ug, Gg, Wg, eps), Wg)
LW <- isoTensor:::.objIsoTensor(Xg, Ag, Ug, Gg,
    isoTensor:::.updateW_tensor(Xg, Ag, Ug, Gg, Wg, eps))
expect_true(LU <= L0 + 1e-12)
expect_true(LG <= L0 + 1e-12)
expect_true(LW <= L0 + 1e-12)

# Normalization leaves B and Xhat unchanged (numerical stabilization only)
nrm <- isoTensor:::.normalizeFactorsTensor(Ug, Gg, Wg)
expect_equal(isoTensor:::.recTensorPartialTucker(nrm$U, nrm$G, nrm$W),
    isoTensor:::.recTensorPartialTucker(Ug, Gg, Wg), tolerance=1e-12)
expect_equal(isoTensor:::.XhatTensor(Ag, nrm$U, nrm$G, nrm$W),
    isoTensor:::.XhatTensor(Ag, Ug, Gg, Wg), tolerance=1e-12)
expect_equal(sqrt(colSums(nrm$U^2)), rep(1, Rn))
expect_equal(sqrt(colSums(nrm$W^2)), rep(1, Ri))

# Rank validation
expect_error(isoTensor(toy$X, toy$A, model="tensor"))
expect_error(isoTensor(toy$X, toy$A, model="tensor", rank_individual=1))
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=0, rank_isoform=2))
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=1.5))
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=21, rank_isoform=2))
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=21))

# R_N >= R_I is NOT an error: it is a generically non-identifiable
# regime for the latent B, flagged by a warning
expect_warning(
    isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
        rank_individual=2, rank_isoform=2, num.iter=2),
    "non-identifiable")
# N below the Theorem C threshold N* is also warned
toy_small <- isoToyModel(model="tensor", N=5, C=3, I=12,
    rank_individual=2, rank_isoform=3, A_variation="high", seed=3)
expect_warning(
    isoTensor(toy_small$X, toy_small$A, model="tensor", optimizer="mu",
        rank_individual=2, rank_isoform=3, num.iter=2),
    "threshold")

# init factors are respected and validated
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2, initU=matrix(1, 3, 1)))
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2, initU=matrix(-1, 20, 1)))

# Matrix mode regression: unchanged behavior, rank args ignored
toy_m <- isoToyModel(N=5, C=2, I=6, seed=1234)
set.seed(1234)
fit_m1 <- isoTensor(toy_m$X, toy_m$A, model="matrix", num.iter=100)
set.seed(1234)
fit_m2 <- isoTensor(toy_m$X, toy_m$A, model="matrix",
    rank_individual=3, rank_isoform=5, num.iter=100)
expect_identical(fit_m1$B, fit_m2$B)
expect_equivalent(names(fit_m1),
    c("B", "model", "algorithm", "RecError", "RelChange"))

# Explicitly supplied optimizer is ignored (with a warning) in matrix mode
expect_warning(
    isoTensor(toy_m$X, toy_m$A, model="matrix", optimizer="lbfgs",
        num.iter=10),
    "ignored")
