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

## model="tensor"

tt <- isoToyModel(model="tensor", N=20, C=3, I=20,
    rank_individual=2, rank_isoform=2, A_variation="moderate", seed=1)

# Dimensions
expect_equivalent(dim(tt$X), c(20, 20))
expect_equivalent(dim(tt$A), c(20, 3))
expect_equivalent(dim(tt$B), c(20, 3, 20))
expect_equivalent(dim(tt$U), c(20, 2))
expect_equivalent(dim(tt$G), c(2, 3, 2))
expect_equivalent(dim(tt$W), c(20, 2))

# Non-negativity
expect_true(min(tt$X) >= 0)
expect_true(min(tt$A) >= 0)
expect_true(min(tt$B) >= 0)

# Cell type fractions sum to 1
expect_equal(rowSums(tt$A), rep(1, 20))

# Noise-free X is exactly the contraction sum_c A_nc B_nci
X_rec <- matrix(0, nrow=20, ncol=20)
for(c in seq_len(3)){
    X_rec <- X_rec + tt$A[, c] * tt$B[, c, ]
}
expect_equal(tt$X, X_rec)

# B is exactly the partial Tucker reconstruction of (U, G, W)
B_rec <- array(0, dim=c(20, 3, 20))
for(c in seq_len(3)){
    B_rec[, c, ] <- tt$U %*% matrix(tt$G[, c, ], nrow=2, ncol=2) %*% t(tt$W)
}
expect_equal(tt$B, B_rec)

# Reproducibility with fixed seed
tt2 <- isoToyModel(model="tensor", N=20, C=3, I=20,
    rank_individual=2, rank_isoform=2, A_variation="moderate", seed=1)
expect_identical(tt$X, tt2$X)
expect_identical(tt$B, tt2$B)

# individual_effect=FALSE: B_nci does not depend on n
tt3 <- isoToyModel(model="tensor", N=20, C=3, I=20,
    individual_effect=FALSE, A_variation="moderate", seed=1)
expect_true(max(abs(sweep(tt3$B, c(2, 3), tt3$B[1, , ]))) == 0)

# effect_cell_types: individual variation restricted to specified cell types
tt4 <- isoToyModel(model="tensor", N=20, C=3, I=20,
    individual_effect=TRUE, effect_cell_types=1,
    A_variation="high", seed=2)
dev_by_celltype <- vapply(seq_len(3), function(c){
    max(abs(sweep(tt4$B[, c, , drop=FALSE], c(2, 3), tt4$B[1, c, ])))
}, numeric(1))
expect_true(dev_by_celltype[1] > 0)
expect_equal(dev_by_celltype[2], 0)
expect_equal(dev_by_celltype[3], 0)
expect_error(isoToyModel(model="tensor", N=20, C=3, I=20,
    rank_individual=1, effect_cell_types=1))

# A_variation levels: "none" gives identical rows / rank 1,
# variation increases monotonically with the level
tt5 <- isoToyModel(model="tensor", N=50, C=3, I=5, A_variation="none", seed=3)
expect_equal(tt5$A_summary$rank, 1)
expect_equal(tt5$A_summary$total_variation, 0)
tv <- vapply(c("low", "moderate", "high"), function(av){
    isoToyModel(model="tensor", N=50, C=3, I=5,
        A_variation=av, seed=3)$A_summary$total_variation
}, numeric(1))
expect_true(tv[1] < tv[2])
expect_true(tv[2] < tv[3])
# Numeric A_variation = total Dirichlet concentration
tt6 <- isoToyModel(model="tensor", N=50, C=3, I=5, A_variation=20, seed=3)
expect_true(tt6$A_summary$rank == 3)
expect_error(isoToyModel(model="tensor", N=10, C=3, I=5, A_variation="foo"))

# A_mean controls rare cell types
tt7 <- isoToyModel(model="tensor", N=50, C=3, I=5,
    A_variation="moderate", A_mean=c(0.9, 0.05, 0.05), seed=3)
expect_true(which.min(tt7$A_summary$mean_fraction) %in% c(2, 3))

# noise_A: A_obs is perturbed but stays a valid fraction matrix;
# X is generated from the true A
tt8 <- isoToyModel(model="tensor", N=10, C=3, I=5, noise_A=0.05, seed=4)
expect_false(identical(tt8$A, tt8$A_obs))
expect_true(min(tt8$A_obs) >= 0)
expect_equal(rowSums(tt8$A_obs), rep(1, 10))
tt9 <- isoToyModel(model="tensor", N=10, C=3, I=5, noise_A=0, seed=4)
expect_identical(tt9$A, tt9$A_obs)
