# Dimensions
toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
expect_equivalent(dim(toy$X), c(5, 6))
expect_equivalent(dim(toy$A), c(5, 2))
expect_equivalent(dim(toy$B), c(2, 6))

# Non-negativity
expect_true(min(toy$X) >= 0)
expect_true(min(toy$A) >= 0)
expect_true(min(toy$B) >= 0)

# Cell type fractions sum to 1
expect_equal(rowSums(toy$A), rep(1, 5))

# Noiseless X is exactly A %*% B
expect_equal(toy$X, toy$A %*% toy$B)

# Reproducibility with fixed seed
toy2 <- isoToyModel(N=5, C=2, I=6, seed=1234)
expect_identical(toy$X, toy2$X)
expect_identical(toy$A, toy2$A)
expect_identical(toy$B, toy2$B)

# Noisy data is still non-negative and differs from noiseless
toy3 <- isoToyModel(N=5, C=2, I=6, noise=0.05, seed=1234)
expect_true(min(toy3$X) >= 0)
expect_false(identical(toy3$X, toy3$A %*% toy3$B))

# Single cell type
toy4 <- isoToyModel(N=5, C=1, I=6, seed=1234)
expect_equivalent(dim(toy4$A), c(5, 1))
expect_equal(as.vector(toy4$A), rep(1, 5))
