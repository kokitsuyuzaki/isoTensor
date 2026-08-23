# Numerical local-identifiability check for the partial Tucker model
#   B_nci = sum_pq U_np G_pcq W_iq,  X_ni = sum_c A_nc B_nci  (A known)
#
# Development script only (NOT part of the package API).
# Numerical Jacobian ranks are LOCAL / GENERIC EVIDENCE, not proofs.
#
# Method: exact Jacobians of theta = (U, G, W) -> vec(B) and -> vec(X).
#   - ker J_B contains the gauge directions (GL(R_N) x GL(R_I)),
#     so generically rank(J_B) = P - R_N^2 - R_I^2.
#   - B is locally identifiable from X at theta iff ker J_X = ker J_B,
#     i.e., deficiency := rank(J_B) - rank(J_X) == 0.
#   - deficiency > 0 certifies a tangent direction that changes B but
#     not X (local NON-identifiability of the latent tensor).
#
# Run: Rscript dev/identifiability_check.R

set.seed(1)

## ---- model construction ----------------------------------------------

make_theta <- function(N, C, I, R_N, R_I){
    list(U = matrix(runif(N * R_N, 0.5, 1.5), N, R_N),
         G = array(runif(R_N * C * R_I, 0.5, 1.5), dim = c(R_N, C, R_I)),
         W = matrix(runif(I * R_I, 0.5, 1.5), I, R_I))
}

make_A <- function(N, C, type = c("dirichlet", "identical", "rankdef"),
    alpha0 = 20){
    type <- match.arg(type)
    if(type == "identical"){
        return(matrix(1 / C, N, C))
    }
    if(type == "rankdef"){
        # rows in a 2-dimensional affine subspace of the simplex
        stopifnot(C >= 3)
        t1 <- runif(N)
        v1 <- c(0.7, 0.2, 0.1, rep(0, C - 3))
        v2 <- c(0.1, 0.3, 0.6, rep(0, C - 3))
        A <- outer(t1, v1) + outer(1 - t1, v2)
        return(A / rowSums(A))
    }
    A <- t(vapply(seq_len(N), function(n){
        g <- rgamma(C, shape = alpha0 / C)
        g / sum(g)
    }, numeric(C)))
    A
}

reconstruct_B <- function(U, G, W){
    N <- nrow(U); C <- dim(G)[2]; I <- nrow(W)
    R_N <- ncol(U); R_I <- ncol(W)
    B <- array(0, dim = c(N, C, I))
    for(c in seq_len(C)){
        Gc <- matrix(G[, c, ], R_N, R_I)
        B[, c, ] <- U %*% Gc %*% t(W)
    }
    B
}

contract_X <- function(A, B){
    N <- dim(B)[1]; C <- dim(B)[2]; I <- dim(B)[3]
    X <- matrix(0, N, I)
    for(c in seq_len(C)){
        X <- X + A[, c] * B[, c, ]
    }
    X
}

## ---- exact Jacobians --------------------------------------------------
## theta layout: (vec(U), vec(G), vec(W)) column-major

jacobians <- function(U, G, W, A){
    N <- nrow(U); R_N <- ncol(U)
    C <- dim(G)[2]; R_I <- dim(G)[3]
    I <- nrow(W)
    nU <- N * R_N; nG <- R_N * C * R_I; nW <- I * R_I
    P <- nU + nG + nW
    idxU <- function(n, p) n + (p - 1) * N
    idxG <- function(p, c, q) nU + p + (c - 1) * R_N + (q - 1) * R_N * C
    idxW <- function(i, q) nU + nG + i + (q - 1) * I
    # per-individual contractions K_n = G x_2 a_n, and H_nq = u_n' K_n
    Klist <- lapply(seq_len(N), function(n){
        K <- matrix(0, R_N, R_I)
        for(c in seq_len(C)){
            K <- K + A[n, c] * matrix(G[, c, ], R_N, R_I)
        }
        K
    })
    H <- t(vapply(seq_len(N), function(n) as.vector(U[n, ] %*% Klist[[n]]),
        numeric(R_I)))
    H <- matrix(H, N, R_I)

    ## J_X (NI x P)
    JX <- matrix(0, N * I, P)
    rowX <- function(n, i) n + (i - 1) * N
    for(n in seq_len(N)){
        KW <- Klist[[n]] %*% t(W)   # R_N x I
        for(i in seq_len(I)){
            r <- rowX(n, i)
            for(p in seq_len(R_N)){
                JX[r, idxU(n, p)] <- KW[p, i]
            }
            for(p in seq_len(R_N)) for(c in seq_len(C)) for(q in seq_len(R_I)){
                JX[r, idxG(p, c, q)] <- U[n, p] * A[n, c] * W[i, q]
            }
            for(q in seq_len(R_I)){
                JX[r, idxW(i, q)] <- H[n, q]
            }
        }
    }

    ## J_B (NCI x P)
    JB <- matrix(0, N * C * I, P)
    rowB <- function(n, c, i) n + (c - 1) * N + (i - 1) * N * C
    UG <- lapply(seq_len(C), function(c) U %*% matrix(G[, c, ], R_N, R_I)) # N x R_I
    GW <- lapply(seq_len(C), function(c) matrix(G[, c, ], R_N, R_I) %*% t(W)) # R_N x I
    for(n in seq_len(N)) for(c in seq_len(C)) for(i in seq_len(I)){
        r <- rowB(n, c, i)
        for(p in seq_len(R_N)){
            JB[r, idxU(n, p)] <- GW[[c]][p, i]
        }
        for(p in seq_len(R_N)) for(q in seq_len(R_I)){
            JB[r, idxG(p, c, q)] <- U[n, p] * W[i, q]
        }
        for(q in seq_len(R_I)){
            JB[r, idxW(i, q)] <- UG[[c]][n, q]
        }
    }
    list(JX = JX, JB = JB, P = P)
}

numrank <- function(J){
    d <- svd(J, nu = 0, nv = 0)$d
    if(max(d) == 0) return(0)
    sum(d > max(dim(J)) * .Machine$double.eps * max(d))
}

check_config <- function(N, C, I, R_N, R_I, Atype = "dirichlet",
    label = "", reps = 2){
    res <- vapply(seq_len(reps), function(r){
        th <- make_theta(N, C, I, R_N, R_I)
        A <- make_A(N, C, Atype)
        J <- jacobians(th$U, th$G, th$W, A)
        c(numrank(J$JB), numrank(J$JX), J$P)
    }, numeric(3))
    rankB <- max(res[1, ]); rankX <- max(res[2, ]); P <- res[3, 1]
    gauge <- R_N^2 + R_I^2
    # two-stage counting threshold on N (valid heuristic only for R_N < R_I)
    Nstar <- if(R_I > R_N) ceiling(R_N * (C * R_I - R_N) / (R_I - R_N)) else Inf
    cat(sprintf(
        "%-28s N=%2d C=%d I=%2d R_N=%d R_I=%d | P=%3d P-gauge=%3d NI=%3d | rk(J_B)=%3d rk(J_X)=%3d defic=%2d | N*=%s\n",
        label, N, C, I, R_N, R_I, P, P - gauge, N * I,
        rankB, rankX, rankB - rankX,
        ifelse(is.finite(Nstar), sprintf("%d", Nstar), "-")))
    invisible(rankB - rankX)
}

cat("== Local identifiability of the latent tensor B (deficiency = rank J_B - rank J_X;\n")
cat("   deficiency > 0 certifies local NON-identifiability; 0 = generic local evidence FOR) ==\n\n")

cat("-- R_N >= R_I (theorem predicts deficiency > 0 regardless of A) --\n")
check_config( 4, 2,  4, 1, 1, label = "R_N=R_I=1")
check_config( 6, 2,  6, 2, 1, label = "R_N>R_I (2,1)")
check_config( 8, 2,  8, 2, 2, label = "R_N=R_I=2, C=2")
check_config(20, 3, 20, 2, 2, label = "R_N=R_I=2, C=3 (old grid)")
check_config(10, 2, 10, 3, 2, label = "R_N>R_I (3,2)")

cat("\n-- R_N < R_I: counting threshold N* = R_N(C R_I - R_N)/(R_I - R_N) --\n")
check_config( 2, 2,  6, 1, 2, label = "(1,2) C=2, N<N*")
check_config( 4, 2,  6, 1, 2, label = "(1,2) C=2, N>N*")
check_config( 6, 2, 10, 2, 3, label = "(2,3) C=2, N<N*")
check_config(10, 2, 10, 2, 3, label = "(2,3) C=2, N>N*")
check_config(12, 3, 12, 2, 3, label = "(2,3) C=3, N<N*")
check_config(16, 3, 12, 2, 3, label = "(2,3) C=3, N>N*")

cat("\n-- Degenerate A (theorems predict deficiency > 0 even for R_N < R_I) --\n")
check_config(16, 3, 12, 2, 3, Atype = "identical", label = "(2,3) C=3, A identical")
check_config(16, 3, 12, 2, 3, Atype = "rankdef",  label = "(2,3) C=3, rank(A)=2")

## ---- explicit non-negative alternative solution (R_N = R_I regime) ----
cat("\n== Explicit alternative solution under non-negativity (R_N=R_I=2, C=3, N=20) ==\n")
N <- 20; C <- 3; I <- 20; R_N <- 2; R_I <- 2
th <- make_theta(N, C, I, R_N, R_I)
A <- make_A(N, C, "dirichlet")
B <- reconstruct_B(th$U, th$G, th$W)
X <- contract_X(A, B)
E <- array(runif(R_N * C * R_I, -1, 1), dim = c(R_N, C, R_I))
# u'_n -> u_n as eps -> 0 (each K_n invertible), so a small enough eps
# keeps U' > 0; how small depends on cond(K_n) (shrink adaptively)
eps <- 0.05
repeat{
    Gp <- th$G + eps * E
    Up <- matrix(0, N, R_N)
    for(n in seq_len(N)){
        K  <- matrix(0, R_N, R_I); Kp <- matrix(0, R_N, R_I)
        for(c in seq_len(C)){
            K  <- K  + A[n, c] * matrix(th$G[, c, ], R_N, R_I)
            Kp <- Kp + A[n, c] * matrix(Gp[, c, ],  R_N, R_I)
        }
        h <- as.vector(th$U[n, ] %*% K)
        # solve u' K' = h (K' square here)
        Up[n, ] <- solve(t(Kp), h)
    }
    if(min(Up) > 0 && min(Gp) > 0) break
    eps <- eps / 2
    if(eps < 1e-8) stop("no strictly positive perturbation found")
}
Bp <- reconstruct_B(Up, Gp, th$W)
Xp <- contract_X(A, Bp)
cat(sprintf("eps=%.2e  min(U')=%.4f  min(G')=%.4f  (both must be > 0)\n",
    eps, min(Up), min(Gp)))
cat(sprintf("||X'-X||_F / ||X||_F = %.3e   (must be ~ 0)\n",
    norm(Xp - X, "F") / norm(X, "F")))
cat(sprintf("||B'-B||_F / ||B||_F = %.3e   (non-negligible => B changed)\n",
    norm(matrix(Bp - B), "F") / norm(matrix(B), "F")))
