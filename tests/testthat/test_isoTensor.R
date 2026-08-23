toy <- isoToyModel(N=5, C=2, I=6, seed=1234)

set.seed(1234)
fit <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000, thr=1e-12)

# Output structure
expect_true(is.list(fit))
expect_equivalent(names(fit),
    c("B", "model", "algorithm", "RecError", "RelChange"))
expect_equivalent(fit$model, "matrix")

# Dimensions (B: cell type x isoform)
expect_equivalent(dim(fit$B), c(2, 6))

# Non-negativity
expect_true(min(fit$B) >= 0)

# Reproducibility with fixed seed
set.seed(1234)
fit_again <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000, thr=1e-12)
expect_identical(fit$B, fit_again$B)

# Objective (reconstruction error) never increases
expect_true(all(diff(fit$RecError[-1]) <= 1e-10))

# Noiseless synthetic data: B_true is recovered almost exactly
expect_true(max(abs(fit$B - toy$B)) < 1e-8)
expect_true(tail(fit$RecError, 1) < 1e-8)

# Noisy synthetic data: B is close to B_true, residual stays near noise level
toy_noisy <- isoToyModel(N=5, C=2, I=6, noise=0.05, seed=42)
set.seed(42)
fit_noisy <- isoTensor(toy_noisy$X, toy_noisy$A, num.iter=1000, thr=1e-12)
expect_true(norm(fit_noisy$B - toy_noisy$B, "F") / norm(toy_noisy$B, "F") < 0.2)
expect_true(all(diff(fit_noisy$RecError[-1]) <= 1e-10))

# Invalid dimensions
expect_error(isoTensor(toy$X, toy$A[seq_len(3), ]))
expect_error(isoTensor(toy$X, toy$A, initB=matrix(1, nrow=3, ncol=6)))
expect_error(isoTensor(toy$X, toy$A, initB=matrix(1, nrow=2, ncol=7)))

# Negative input rejection
X_neg <- toy$X
X_neg[1, 1] <- -1
expect_error(isoTensor(X_neg, toy$A))
A_neg <- toy$A
A_neg[1, 1] <- -1
expect_error(isoTensor(toy$X, A_neg))
expect_error(isoTensor(toy$X, toy$A, initB=matrix(-1, nrow=2, ncol=6)))

# Zero entries in X are allowed
X_zero <- toy$X
X_zero[1, 1] <- 0
X_zero[3, 4] <- 0
set.seed(1234)
fit_zero <- isoTensor(X_zero, toy$A, num.iter=100)
expect_true(min(fit_zero$B) >= 0)
expect_true(all(diff(fit_zero$RecError[-1]) <= 1e-10))

# Single cell type: closed-form least squares solution is recovered
toy1 <- isoToyModel(N=5, C=1, I=6, seed=1234)
set.seed(1234)
fit1 <- isoTensor(toy1$X, toy1$A, num.iter=1000, thr=1e-12)
B_ls <- crossprod(toy1$A, toy1$X) / sum(toy1$A * toy1$A)
expect_true(max(abs(fit1$B - B_ls)) < 1e-8)

# model="tensor" is not implemented yet
expect_error(isoTensor(toy$X, toy$A, model="tensor"))

# initB is respected
set.seed(1234)
fit_init <- isoTensor(toy$X, toy$A, initB=matrix(1, nrow=2, ncol=6),
    num.iter=1000, thr=1e-12)
expect_true(max(abs(fit_init$B - toy$B)) < 1e-8)
