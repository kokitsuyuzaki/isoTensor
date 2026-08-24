# optimizer="lbfgs" / "hybrid" for model="tensor"
# (log-reparameterized L-BFGS; see notes/tensor-LBFGS.md)

toy <- isoToyModel(model="tensor", N=20, C=2, I=20,
    rank_individual=1, rank_isoform=2,
    individual_effect=TRUE, A_variation="high", seed=10)
relerr <- function(a, b) sqrt(sum((a - b)^2)) / sqrt(sum(b^2))

## ---- optimizer="lbfgs" (the default) ------------------------------------
set.seed(10)
fit <- isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2)

# The default optimizer is lbfgs
expect_equivalent(fit$optimizer, "lbfgs")

# Output structure
expect_equivalent(names(fit),
    c("U", "G", "W", "B", "Xhat", "model", "algorithm", "optimizer",
      "rank_individual", "rank_isoform", "objective", "RecError", "lbfgs"))
expect_equivalent(names(fit$RecError), c("offset", "final"))

# Dimensions
expect_equivalent(dim(fit$U), c(20, 1))
expect_equivalent(dim(fit$G), c(1, 2, 2))
expect_equivalent(dim(fit$W), c(20, 2))
expect_equivalent(dim(fit$B), c(20, 2, 20))
expect_equivalent(dim(fit$Xhat), c(20, 20))

# Non-negativity (exp map guarantees it by construction)
expect_true(min(fit$U) >= 0)
expect_true(min(fit$G) >= 0)
expect_true(min(fit$W) >= 0)
expect_true(min(fit$B) >= 0)
expect_true(min(fit$Xhat) >= 0)

# Fixed seed reproducibility
set.seed(10)
fit2 <- isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2)
expect_identical(fit$B, fit2$B)

# Objective is reduced from the starting point and matches RecError
expect_true(fit$RecError["final"] < fit$RecError["offset"])
expect_equal(fit$objective, 0.5 * fit$RecError[["final"]]^2)

# Convergence metadata
expect_true(fit$lbfgs$convergence %in% c(0L, 1L))
expect_false(fit$lbfgs$failed)
expect_true(is.numeric(fit$lbfgs$counts))
expect_true(is.character(fit$lbfgs$message))

# Positive control (identifiable candidate regime R_N < R_I): both X and
# the latent tensor B are recovered
expect_true(relerr(fit$Xhat, toy$X) < 1e-3)
expect_true(relerr(fit$B, toy$B) < 0.05)

## ---- pack/unpack roundtrip ----------------------------------------------
set.seed(123)
Up <- matrix(runif(4 * 2, 0.5, 1.5), 4, 2)
Gp <- array(runif(2 * 3 * 2, 0.5, 1.5), dim=c(2, 3, 2))
Wp <- matrix(runif(5 * 2, 0.5, 1.5), 5, 2)
th <- isoTensor:::.packThetaTensor(Up, Gp, Wp)
expect_equivalent(length(th), 4*2 + 2*3*2 + 5*2)
p <- isoTensor:::.unpackThetaTensor(th, N=4, I=5, C=3, R_N=2, R_I=2)
expect_equivalent(dim(p$U), dim(Up))
expect_equivalent(dim(p$G), dim(Gp))
expect_equivalent(dim(p$W), dim(Wp))
expect_equal(p$U, Up, tolerance=1e-12)
expect_equal(p$G, Gp, tolerance=1e-12)
expect_equal(p$W, Wp, tolerance=1e-12)
# Reconstruction from packed/unpacked factors is consistent
expect_equal(isoTensor:::.recTensorPartialTucker(p$U, p$G, p$W),
    isoTensor:::.recTensorPartialTucker(Up, Gp, Wp), tolerance=1e-12)

# Analytical gradient in theta vs central finite differences
Ap <- matrix(runif(4 * 3, 0.2, 1), 4, 3); Ap <- Ap / rowSums(Ap)
Xp <- matrix(runif(4 * 5, 0, 2), 4, 5)
fn <- function(t){
    q <- isoTensor:::.unpackThetaTensor(t, N=4, I=5, C=3, R_N=2, R_I=2)
    isoTensor:::.objIsoTensor(Xp, Ap, q$U, q$G, q$W)
}
g_ana <- {
    q <- isoTensor:::.unpackThetaTensor(th, N=4, I=5, C=3, R_N=2, R_I=2)
    g <- isoTensor:::.gradIsoTensor(Xp, Ap, q$U, q$G, q$W)
    c(as.vector(g$U * q$U), as.vector(g$G * q$G), as.vector(g$W * q$W))
}
h <- 1e-6
g_fd <- sapply(seq_along(th), function(j){
    tp <- th; tp[j] <- tp[j] + h
    tm <- th; tm[j] <- tm[j] - h
    (fn(tp) - fn(tm)) / (2 * h)
})
expect_true(max(abs(g_ana - g_fd) / pmax(abs(g_ana) + abs(g_fd), 1e-8)) < 1e-6)

## ---- user-supplied initialization ---------------------------------------
initU <- matrix(1, 20, 1)
initG <- array(1, dim=c(1, 2, 2))
initW <- matrix(1, 20, 2)
fit_init <- isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2,
    initU=initU, initG=initG, initW=initW)
expect_true(relerr(fit_init$Xhat, toy$X) < 1e-3)
# deterministic (no RNG draw with a full user init)
fit_init2 <- isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2,
    initU=initU, initG=initG, initW=initW)
expect_identical(fit_init$B, fit_init2$B)

# Zeros in the init are raised to 1e-10 for lbfgs (NOT locked, unlike MU):
# entries starting at 0 can become positive
initW0 <- initW
initW0[1, ] <- 0
fit_z <- isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2,
    initU=initU, initG=initG, initW=initW0)
expect_true(all(is.finite(fit_z$objective)))
expect_true(max(fit_z$W[1, ]) > 1e-10)

## ---- optimizer.control ---------------------------------------------------
set.seed(10)
fit_ctl <- isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2,
    optimizer.control=list(maxit=10))
expect_equivalent(fit_ctl$lbfgs$maxit, 10)
expect_error(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2,
    optimizer.control=list(bogus=1)))
expect_warning(isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
    rank_individual=1, rank_isoform=2, num.iter=2,
    optimizer.control=list(maxit=10)), "ignored")

## ---- optimizer="hybrid" ---------------------------------------------------
set.seed(10)
fit_h <- isoTensor(toy$X, toy$A, model="tensor", optimizer="hybrid",
    rank_individual=1, rank_isoform=2, num.iter=200)
expect_equivalent(fit_h$optimizer, "hybrid")
expect_equivalent(names(fit_h),
    c("U", "G", "W", "B", "Xhat", "model", "algorithm", "optimizer",
      "rank_individual", "rank_isoform", "objective", "RecError", "lbfgs",
      "RelChange", "NumIter", "Converged"))
expect_equivalent(dim(fit_h$B), c(20, 2, 20))
expect_true(min(fit_h$U) >= 0 && min(fit_h$G) >= 0 && min(fit_h$W) >= 0)
expect_true(min(fit_h$B) >= 0)
# RecError = MU history + one refined value; refinement does not make the
# objective worse than the MU stage
expect_equivalent(tail(names(fit_h$RecError), 1), "refined")
mu_final <- fit_h$RecError[[length(fit_h$RecError) - 1]]
expect_true(fit_h$RecError[["refined"]] <= mu_final * (1 + 1e-8))

## ---- cross-optimizer comparison on the same easy data ---------------------
set.seed(10)
fit_mu <- isoTensor(toy$X, toy$A, model="tensor", optimizer="mu",
    rank_individual=1, rank_isoform=2, num.iter=200)
# lbfgs and hybrid reach at least as good an objective as MU
expect_true(fit$objective <= fit_mu$objective * (1 + 1e-8))
expect_true(fit_h$objective <= fit_mu$objective * (1 + 1e-8))
# and recover the ground-truth latent tensor on this identifiable dataset
expect_true(relerr(fit$B, toy$B) < 0.05)
expect_true(relerr(fit_h$B, toy$B) < 0.05)

# A is unchanged after lbfgs fitting
A_copy <- toy$A + 0
set.seed(10)
invisible(isoTensor(toy$X, toy$A, model="tensor",
    rank_individual=1, rank_isoform=2,
    optimizer.control=list(maxit=20)))
expect_identical(toy$A, A_copy)

# Warnings for theoretically flagged regimes also fire on the lbfgs path
expect_warning(
    isoTensor(toy$X, toy$A, model="tensor",
        rank_individual=2, rank_isoform=2,
        optimizer.control=list(maxit=5)),
    "non-identifiable")
