# =============================================================================
# Port verification: isoTensor(optimizer="lbfgs"/"hybrid") vs the validated
# implementation in ../isoTensor-experiments (src/lbfgs_utils.R et al.)
#
# 1. Numerical regression: refit the optimizer-calibration datasets
#    (easy_high / mod_a20 / hard_low, init seeds 1-5) with
#    isoTensor(optimizer="lbfgs") and compare objective / relX / relB
#    against the archived experiments results
#    (results/optimizer/alt_summary.tsv, method="lbfgs").
#    Note: the archived method="hybrid" rows are NOT a regression target —
#    in that calibration script the MU warm stage was called with thr=0,
#    which under the pre-fix sentinel ran 0 MU iterations, so those rows
#    are effectively "L-BFGS from a rescaled random init" (their ~0.9 s
#    runtimes confirm it). The full-grid / noise-stress hybrids used
#    thr=1e-8 and are genuine.
# 2. Runtime / accuracy comparison of mu / lbfgs / hybrid (item 17).
# 3. thr=0 edge case fix check.
#
# Run from the package root:
#   PATH=/home/koki/anaconda3/envs/r_4.3/bin:$PATH Rscript dev/lbfgs_port_check.R
# =============================================================================

for(f in list.files("R", full.names = TRUE)) source(f)
EXP <- "../isoTensor-experiments"
relerr <- function(a, b) sqrt(sum((a - b)^2)) / sqrt(sum(b^2))

## ---- 1a. Exact identity vs the experiments IMPLEMENTATION --------------
## The experiments L-BFGS code (src/lbfgs_utils.R) is executed locally on
## the same data / seeds and must match the package bitwise (same RNG
## draws, same arithmetic, same optim settings). This is the actual port
## check. The archived TSVs (1b) were produced inside the experiments
## container (different BLAS/R build), where L-BFGS trajectories diverge
## at floating-point level, so only statistical agreement is expected there.
source(file.path(EXP, "src/lbfgs_utils.R"))
cat("== 1a. bitwise identity vs locally-run experiments implementation ==\n")
worst_a <- 0
for(key in c("easy_high", "mod_a20", "hard_low")){
    dat <- readRDS(file.path(EXP, "results/optimizer/data",
        paste0(key, ".rds")))
    m <- dat$meta
    N <- nrow(dat$X); I <- ncol(dat$X); C <- ncol(dat$A_obs)
    R_N <- m$rank_individual; R_I <- m$rank_isoform
    prob <- make_lbfgs_problem(dat$X, dat$A_obs, R_N, R_I)
    for(s in 1:5){
        # experiments implementation, run here
        set.seed(s)
        U0 <- matrix(runif(N * R_N, 0.5, 1.5), N, R_N)
        G0 <- array(runif(R_N * C * R_I, 0.5, 1.5), dim = c(R_N, C, R_I))
        W0 <- matrix(runif(I * R_I, 0.5, 1.5), I, R_I)
        res <- run_lbfgs(prob, prob$pack(U0, G0, W0), 5000, 10)
        p <- prob$unpack(res$theta)
        B_exp <- .recTensorPartialTucker(p$U, p$G, p$W)
        # package port, same seed
        set.seed(s)
        fit <- suppressWarnings(isoTensor(dat$X, dat$A_obs, model = "tensor",
            rank_individual = R_N, rank_isoform = R_I,
            optimizer = "lbfgs",
            optimizer.control = list(maxit = 5000, factr = 10)))
        # package objective is recomputed after the final (invariant)
        # normalization, so compare on the scale of the data, not of the
        # (possibly ~1e-13) objective itself
        d_obj <- abs(fit$objective - res$objective) / (0.5 * sum(dat$X^2))
        d_B <- max(abs(fit$B - B_exp)) / max(abs(B_exp))
        worst_a <- max(worst_a, d_obj, d_B)
        cat(sprintf("%-9s seed=%d  obj %.6e (exp %.6e)  maxrel|dB| %.1e\n",
            key, s, fit$objective, res$objective, d_B))
    }
}
if(worst_a < 1e-12){
    cat("PASS: package lbfgs is numerically identical to the experiments implementation\n\n")
}else{
    stop("port identity check FAILED (max deviation ", worst_a, ")")
}

## ---- 1b. Statistical agreement with the archived container results -----
ref <- read.delim(file.path(EXP, "results/optimizer/alt_summary.tsv"))
ref <- ref[ref$method == "lbfgs", ]
cat("== 1b. vs archived container TSV (statistical agreement only) ==\n")
for(key in c("easy_high", "mod_a20", "hard_low")){
    dat <- readRDS(file.path(EXP, "results/optimizer/data",
        paste0(key, ".rds")))
    m <- dat$meta
    relBs <- sapply(1:5, function(s){
        set.seed(s)
        fit <- suppressWarnings(isoTensor(dat$X, dat$A_obs, model = "tensor",
            rank_individual = m$rank_individual,
            rank_isoform = m$rank_isoform, optimizer = "lbfgs"))
        relerr(fit$B, dat$B)
    })
    r <- ref[ref$dataset_key == key, ]
    cat(sprintf("%-9s median relB: package %.3e vs archived %.3e\n",
        key, median(relBs), median(r$relB)))
}
cat("\n")

## ---- 2. Runtime / accuracy: mu vs lbfgs vs hybrid (easy_high) ----------
cat("== 2. runtime comparison (easy_high, seed 1) ==\n")
dat <- readRDS(file.path(EXP, "results/optimizer/data/easy_high.rds"))
m <- dat$meta
run1 <- function(expr){
    t0 <- proc.time()
    fit <- suppressWarnings(expr)
    list(fit = fit, sec = (proc.time() - t0)[["elapsed"]])
}
set.seed(1)
r_mu <- run1(isoTensor(dat$X, dat$A_obs, model = "tensor",
    rank_individual = m$rank_individual, rank_isoform = m$rank_isoform,
    optimizer = "mu", num.iter = 20000, thr = 1e-8))
set.seed(1)
r_lb <- run1(isoTensor(dat$X, dat$A_obs, model = "tensor",
    rank_individual = m$rank_individual, rank_isoform = m$rank_isoform,
    optimizer = "lbfgs"))
set.seed(1)
r_hy <- run1(isoTensor(dat$X, dat$A_obs, model = "tensor",
    rank_individual = m$rank_individual, rank_isoform = m$rank_isoform,
    optimizer = "hybrid", num.iter = 2000, thr = 1e-8))
for(nm in c("mu", "lbfgs", "hybrid")){
    r <- switch(nm, mu = r_mu, lbfgs = r_lb, hybrid = r_hy)
    cat(sprintf("%-6s obj=%.3e relX=%.3e relB=%.3e  %.2fs\n",
        nm, r$fit$objective, relerr(r$fit$Xhat, dat$X),
        relerr(r$fit$B, dat$B), r$sec))
}
stopifnot(r_lb$fit$objective <= r_mu$fit$objective * (1 + 1e-6) ||
          r_lb$fit$objective < 1e-6)
stopifnot(r_hy$fit$objective <= r_mu$fit$objective * (1 + 1e-6) ||
          r_hy$fit$objective < 1e-6)
cat("\n")

## ---- 3. thr=0 edge case -------------------------------------------------
cat("== 3. thr=0 edge case ==\n")
toy <- isoToyModel(N = 5, C = 2, I = 6, seed = 1)
set.seed(1)
fm <- isoTensor(toy$X, toy$A, model = "matrix", thr = 0, num.iter = 5)
stopifnot(length(fm$RecError) == 6)  # offset + 5 iterations
tt <- isoToyModel(model = "tensor", N = 10, C = 2, I = 10,
    rank_individual = 1, rank_isoform = 2, A_variation = "high", seed = 2)
set.seed(2)
ft <- suppressWarnings(isoTensor(tt$X, tt$A, model = "tensor",
    rank_individual = 1, rank_isoform = 2, optimizer = "mu",
    thr = 0, num.iter = 5))
stopifnot(ft$NumIter == 5)
cat("PASS: thr=0 now runs the requested iterations\n")
