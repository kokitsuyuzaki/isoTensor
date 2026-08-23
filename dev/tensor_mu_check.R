# Step-by-step verification of the Phase 2 tensor MU (debugging policy:
# check small pieces in order; never debug by running full fitting only).
# Development script, not part of the package API.
#
# Steps (notes/tensor-MU.md section 7 / task instruction section 11):
#   1. tensor reconstruction        2. contraction with fixed A
#   3. objective                    4. analytical gradient
#   5. finite-difference gradient   6. MU 1 iteration
#   7. objective before/after       8. 10 iterations
#   9. convergence                 10. ground-truth B recovery
# plus positive/negative controls and multi-start.
#
# Run: Rscript dev/tensor_mu_check.R

for(f in list.files("R", full.names = TRUE)) source(f)

relerr <- function(a, b) norm(as.matrix(as.vector(a - b)), "F") /
    norm(as.matrix(as.vector(b)), "F")

## ---- Steps 1-3: reconstruction, contraction, objective ---------------
set.seed(1)
N <- 3; C <- 2; I <- 4; R_N <- 1; R_I <- 2
U <- matrix(runif(N * R_N, 0.5, 1.5), N, R_N)
G <- array(runif(R_N * C * R_I, 0.5, 1.5), dim = c(R_N, C, R_I))
W <- matrix(runif(I * R_I, 0.5, 1.5), I, R_I)
A <- matrix(runif(N * C, 0.2, 1), N, C); A <- A / rowSums(A)
X <- matrix(runif(N * I, 0, 2), N, I)

# 1) reconstruction: .recTensorPartialTucker vs naive triple loop
B1 <- .recTensorPartialTucker(U, G, W)
B0 <- array(0, dim = c(N, C, I))
for(n in seq_len(N)) for(c in seq_len(C)) for(i in seq_len(I)){
    for(p in seq_len(R_N)) for(q in seq_len(R_I)){
        B0[n, c, i] <- B0[n, c, i] + U[n, p] * G[p, c, q] * W[i, q]
    }
}
cat("step1 reconstruction max|diff| =", max(abs(B1 - B0)), "\n")

# 2) contraction with fixed A: .XhatTensor vs naive
Xh1 <- .XhatTensor(A, U, G, W)
Xh0 <- matrix(0, N, I)
for(n in seq_len(N)) for(i in seq_len(I)) for(c in seq_len(C)){
    Xh0[n, i] <- Xh0[n, i] + A[n, c] * B0[n, c, i]
}
cat("step2 contraction   max|diff| =", max(abs(Xh1 - Xh0)), "\n")

# 3) objective
L1 <- .objIsoTensor(X, A, U, G, W)
cat("step3 objective =", L1, "(naive:", 0.5 * sum((X - Xh0)^2), ")\n")

## ---- Steps 4-5: analytical vs finite-difference gradient -------------
fd_check <- function(X, A, U, G, W, h = 1e-6){
    g <- .gradIsoTensor(X, A, U, G, W)
    fd_one <- function(get, set){
        th <- get()
        fd <- array(0, dim = dim(th))
        for(j in seq_along(th)){
            tp <- th; tp[j] <- tp[j] + h
            tm <- th; tm[j] <- tm[j] - h
            fd[j] <- (set(tp) - set(tm)) / (2 * h)
        }
        fd
    }
    fdU <- fd_one(function() U, function(t) .objIsoTensor(X, A, t, G, W))
    fdG <- fd_one(function() G, function(t) .objIsoTensor(X, A, U, t, W))
    fdW <- fd_one(function() W, function(t) .objIsoTensor(X, A, U, G, t))
    rel <- function(a, f) max(abs(a - f) / pmax(abs(a) + abs(f), 1e-8))
    c(U = rel(g$U, fdU), G = rel(g$G, fdG), W = rel(g$W, fdW))
}
r1 <- fd_check(X, A, U, G, W)
cat("step4-5 gradient vs FD (N=3,C=2,I=4,R_N=1,R_I=2), max rel err:\n")
print(r1)
# second config (R_N=2, R_I=3, C=3)
set.seed(2)
N2 <- 4; C2 <- 3; I2 <- 5; Rn2 <- 2; Ri2 <- 3
U2 <- matrix(runif(N2 * Rn2, 0.5, 1.5), N2, Rn2)
G2 <- array(runif(Rn2 * C2 * Ri2, 0.5, 1.5), dim = c(Rn2, C2, Ri2))
W2 <- matrix(runif(I2 * Ri2, 0.5, 1.5), I2, Ri2)
A2 <- matrix(runif(N2 * C2, 0.2, 1), N2, C2); A2 <- A2 / rowSums(A2)
X2 <- matrix(runif(N2 * I2, 0, 2), N2, I2)
r2 <- fd_check(X2, A2, U2, G2, W2)
cat("second config (N=4,C=3,I=5,R_N=2,R_I=3), max rel err:\n")
print(r2)
stopifnot(all(c(r1, r2) < 1e-6))
cat("GRADIENT CHECK PASSED (tolerance 1e-6)\n\n")

## ---- Steps 6-9: MU single steps, 10 iterations, convergence ----------
if(!exists(".updateU_tensor")){
    cat("MU functions not implemented yet -- stopping after gradient check\n")
    quit(save = "no")
}
eps <- .Machine$double.eps
L0 <- .objIsoTensor(X, A, U, G, W)
Uu <- .updateU_tensor(X, A, U, G, W, eps)
LU <- .objIsoTensor(X, A, Uu, G, W)
Gg <- .updateG_tensor(X, A, U, G, W, eps)
LG <- .objIsoTensor(X, A, U, Gg, W)
Ww <- .updateW_tensor(X, A, U, G, W, eps)
LW <- .objIsoTensor(X, A, U, G, Ww)
cat(sprintf("step6-7 one-step objective: L0=%.6f  afterU=%.6f  afterG=%.6f  afterW=%.6f\n",
    L0, LU, LG, LW))
stopifnot(LU <= L0 + 1e-12, LG <= L0 + 1e-12, LW <= L0 + 1e-12)

# 8) 10 iterations monotone
Ui <- U; Gi <- G; Wi <- W
Ls <- .objIsoTensor(X, A, Ui, Gi, Wi)
for(it in 1:10){
    Ui <- .updateU_tensor(X, A, Ui, Gi, Wi, eps)
    Gi <- .updateG_tensor(X, A, Ui, Gi, Wi, eps)
    Wi <- .updateW_tensor(X, A, Ui, Gi, Wi, eps)
    Ls <- c(Ls, .objIsoTensor(X, A, Ui, Gi, Wi))
}
cat("step8 objective over 10 iters:", paste(sprintf("%.4f", Ls), collapse = " "), "\n")
stopifnot(all(diff(Ls) <= 1e-12))

# 9) convergence through the API
toy <- isoToyModel(model = "tensor", N = 20, C = 2, I = 20,
    rank_individual = 1, rank_isoform = 2,
    individual_effect = TRUE, A_variation = "high", seed = 10)
set.seed(10)
fit <- isoTensor(toy$X, toy$A, model = "tensor",
    rank_individual = 1, rank_isoform = 2, num.iter = 5000, thr = 1e-12)
cat(sprintf("step9 convergence: NumIter=%d  Converged=%s  final RecError=%.3e\n",
    fit$NumIter, fit$Converged, tail(fit$RecError, 1)))
stopifnot(all(diff(fit$RecError[-1]) <= 1e-10))

# 10) ground-truth B recovery (positive control, R_N < R_I)
cat(sprintf("step10 positive control (R_N=1 < R_I=2): relX=%.3e  relB=%.3e\n",
    relerr(fit$Xhat, toy$X), relerr(fit$B, toy$B)))

## ---- Negative control: proven non-identifiable regime ----------------
toyN <- isoToyModel(model = "tensor", N = 20, C = 3, I = 20,
    rank_individual = 2, rank_isoform = 2,
    individual_effect = TRUE, A_variation = "high", seed = 11)
set.seed(11)
fitN <- suppressWarnings(isoTensor(toyN$X, toyN$A, model = "tensor",
    rank_individual = 2, rank_isoform = 2, num.iter = 5000, thr = 1e-12))
cat(sprintf("negative control (R_N=R_I=2): relX=%.3e  relB=%.3e\n",
    relerr(fitN$Xhat, toyN$X), relerr(fitN$B, toyN$B)))

## ---- Multi-start ------------------------------------------------------
multistart <- function(toy, R_N, R_I, seeds, num.iter = 5000){
    t(vapply(seeds, function(s){
        set.seed(s)
        f <- suppressWarnings(isoTensor(toy$X, toy$A, model = "tensor",
            rank_individual = R_N, rank_isoform = R_I,
            num.iter = num.iter, thr = 1e-12))
        c(seed = s, obj = 0.5 * tail(f$RecError, 1)^2,
          relX = relerr(f$Xhat, toy$X), relB = relerr(f$B, toy$B))
    }, numeric(4)))
}
cat("\nmulti-start, identifiable candidate regime (R_N=1 < R_I=2):\n")
ms1 <- multistart(toy, 1, 2, seeds = 1:10)
print(round(ms1, 6))
cat("\nmulti-start, proved non-identifiable regime (R_N=R_I=2):\n")
ms2 <- multistart(toyN, 2, 2, seeds = 1:10)
print(round(ms2, 6))
cat(sprintf("\nsummary relB: identifiable regime  min=%.2e max=%.2e\n",
    min(ms1[, "relB"]), max(ms1[, "relB"])))
cat(sprintf("summary relB: non-identifiable    min=%.2e max=%.2e (relX all < %.1e)\n",
    min(ms2[, "relB"]), max(ms2[, "relB"]), max(ms2[, "relX"])))
