isoTensor <- function(X, A, model=c("matrix", "tensor"),
    rank_individual=NULL, rank_isoform=NULL,
    initB=NULL, initU=NULL, initG=NULL, initW=NULL,
    pseudocount=.Machine$double.eps,
    thr=1e-10, num.iter=100, verbose=FALSE){
    # Argument check
    model <- match.arg(model)
    .checkIsoTensor(X, A, initB, pseudocount, thr, num.iter, verbose)
    if(model == "tensor"){
        .checkIsoTensorTensor(X, A, rank_individual, rank_isoform,
            initU, initG, initW)
        return(.eachIsoTensorTensor(X, A, rank_individual, rank_isoform,
            initU, initG, initW, pseudocount, thr, num.iter, verbose))
    }
    # model="matrix": initialization of B
    int <- .initIsoTensor(X, A, initB)
    B <- int$B
    AtX <- int$AtX
    AtA <- int$AtA
    RecError <- c()
    RelChange <- c()
    RecError[1] <- .recError(X, A %*% B)
    RelChange[1] <- thr * 10
    # Iterative step
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        pre_Error <- RecError[iter]
        # Update B (A is always fixed)
        B <- .updateB(B, AtX, AtA, pseudocount)
        iter <- iter + 1
        RecError[iter] <- .recError(X, A %*% B)
        if(RecError[iter] == 0){
            RelChange[iter] <- 0
        }else{
            RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        }
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter,
                " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    names(RecError) <- c("offset", seq_len(iter-1))
    names(RelChange) <- c("offset", seq_len(iter-1))
    list(B = B, model = model, algorithm = "Frobenius",
        RecError = RecError, RelChange = RelChange)
}

.checkIsoTensor <- function(X, A, initB, pseudocount, thr, num.iter, verbose){
    stopifnot(is.matrix(X))
    stopifnot(is.matrix(A))
    stopifnot(is.numeric(X))
    stopifnot(is.numeric(A))
    if(nrow(X) != nrow(A)){
        stop("Please specify nrow(X) and nrow(A) are same (= number of individuals)")
    }
    if(any(is.na(X)) || any(is.na(A))){
        stop("NA is not supported in X and A")
    }
    if(min(X) < 0){
        stop("Please specify the non-negative X")
    }
    if(min(A) < 0){
        stop("Please specify the non-negative A")
    }
    if(!is.null(initB)){
        stopifnot(is.matrix(initB))
        if(nrow(initB) != ncol(A)){
            stop("Please specify nrow(initB) and ncol(A) are same (= number of cell types)")
        }
        if(ncol(initB) != ncol(X)){
            stop("Please specify ncol(initB) and ncol(X) are same (= number of isoforms)")
        }
        if(min(initB) < 0){
            stop("Please specify the non-negative initB")
        }
    }
    stopifnot(is.numeric(pseudocount))
    stopifnot(pseudocount >= 0)
    stopifnot(is.numeric(thr))
    stopifnot(is.numeric(num.iter))
    stopifnot(num.iter >= 1)
    stopifnot(is.logical(verbose))
}

.initIsoTensor <- function(X, A, initB){
    C <- ncol(A)
    I <- ncol(X)
    if(is.null(initB)){
        B <- matrix(runif(C * I), nrow = C, ncol = I)
        # Optimal global rescaling argmin_a ||X - a AB||_F
        AB <- A %*% B
        denom <- sum(AB * AB)
        if(denom > 0){
            B <- B * sum(X * AB) / denom
        }
    }else{
        B <- initB
    }
    # A is fixed, so these are computed only once
    AtX <- crossprod(A, X)
    AtA <- crossprod(A)
    list(B = B, AtX = AtX, AtA = AtA)
}

# MU rule derived from grad_B (1/2 ||X - AB||_F^2) = t(A) %*% A %*% B - t(A) %*% X
.updateB <- function(B, AtX, AtA, pseudocount){
    B * AtX / (AtA %*% B + pseudocount)
}

.recError <- function(X, X_bar){
    v <- as.vector(X_bar - X)
    sqrt(sum(v * v))
}

## ---- model="tensor": plain partial non-negative Tucker ----------------
## B_nci = sum_pq U_np G_pcq W_iq,  Xhat_ni = sum_c A_nc B_nci
## Derivation of gradients and MU: notes/tensor-MU.md

# H_nq = sum_{p,c} A_nc U_np G_pcq  (N x R_I); Xhat = H %*% t(W)
.HTensor <- function(A, U, G){
    R_N <- ncol(U)
    R_I <- dim(G)[3]
    C <- ncol(A)
    H <- matrix(0, nrow = nrow(U), ncol = R_I)
    for(c in seq_len(C)){
        Gc <- matrix(G[, c, ], nrow = R_N, ncol = R_I)
        H <- H + A[, c] * (U %*% Gc)
    }
    H
}

.XhatTensor <- function(A, U, G, W){
    .HTensor(A, U, G) %*% t(W)
}

.objIsoTensor <- function(X, A, U, G, W){
    0.5 * sum((X - .XhatTensor(A, U, G, W))^2)
}

# Analytical gradients of L = 1/2 ||X - Xhat||_F^2 (notes/tensor-MU.md)
#   grad_U = sum_c diag(A[,c]) (R W) Gc^T          (N x R_N)
#   grad_Gc = t(U) diag(A[,c]) (R W)               (R_N x R_I per c)
#   grad_W = t(R) H                                (I x R_I)
.gradIsoTensor <- function(X, A, U, G, W){
    R_N <- ncol(U)
    R_I <- ncol(W)
    C <- ncol(A)
    H <- .HTensor(A, U, G)
    R <- H %*% t(W) - X
    RW <- R %*% W
    gU <- matrix(0, nrow = nrow(U), ncol = R_N)
    gG <- array(0, dim = dim(G))
    for(c in seq_len(C)){
        Gc <- matrix(G[, c, ], nrow = R_N, ncol = R_I)
        gU <- gU + A[, c] * (RW %*% t(Gc))
        gG[, c, ] <- crossprod(U, A[, c] * RW)
    }
    gW <- crossprod(R, H)
    list(U = gU, G = gG, W = gW)
}

# MU for U: U * [sum_c diag(A[,c]) (XW) Gc^T] / [sum_c diag(A[,c]) (XhatW) Gc^T]
# (row-wise Lee-Seung NNLS updates with designs W Kn^T; notes/tensor-MU.md)
.updateU_tensor <- function(X, A, U, G, W, pseudocount){
    R_N <- ncol(U)
    R_I <- ncol(W)
    C <- ncol(A)
    XW <- X %*% W
    XhatW <- .XhatTensor(A, U, G, W) %*% W
    numer <- matrix(0, nrow = nrow(U), ncol = R_N)
    denom <- matrix(0, nrow = nrow(U), ncol = R_N)
    for(c in seq_len(C)){
        Gct <- t(matrix(G[, c, ], nrow = R_N, ncol = R_I))
        numer <- numer + A[, c] * (XW %*% Gct)
        denom <- denom + A[, c] * (XhatW %*% Gct)
    }
    U * numer / (denom + pseudocount)
}

# MU for G (slice-wise): Gc * [t(U) diag(A[,c]) XW] / [t(U) diag(A[,c]) XhatW]
# (single NNLS in vec(G) with design W kron (U facesplit A))
.updateG_tensor <- function(X, A, U, G, W, pseudocount){
    R_N <- ncol(U)
    R_I <- ncol(W)
    C <- ncol(A)
    XW <- X %*% W
    XhatW <- .XhatTensor(A, U, G, W) %*% W
    Gnew <- G
    for(c in seq_len(C)){
        Gc <- matrix(G[, c, ], nrow = R_N, ncol = R_I)
        numer <- crossprod(U, A[, c] * XW)
        denom <- crossprod(U, A[, c] * XhatW)
        Gnew[, c, ] <- Gc * numer / (denom + pseudocount)
    }
    Gnew
}

# MU for W: W * (t(X) H) / (W (t(H) H))  (classical NMF one-sided update
# on X^T ~ W H^T with H fixed)
.updateW_tensor <- function(X, A, U, G, W, pseudocount){
    H <- .HTensor(A, U, G)
    numer <- crossprod(X, H)
    denom <- W %*% crossprod(H)
    W * numer / (denom + pseudocount)
}

# Absorb column scales of U and W into G (leaves B and Xhat unchanged;
# numerical stabilization only, see notes/tensor-MU.md section 5)
.normalizeFactorsTensor <- function(U, G, W){
    R_N <- ncol(U)
    R_I <- ncol(W)
    su <- sqrt(colSums(U * U))
    su[su == 0] <- 1
    U <- sweep(U, 2, su, "/")
    for(p in seq_len(R_N)){
        G[p, , ] <- G[p, , ] * su[p]
    }
    sw <- sqrt(colSums(W * W))
    sw[sw == 0] <- 1
    W <- sweep(W, 2, sw, "/")
    for(q in seq_len(R_I)){
        G[, , q] <- G[, , q] * sw[q]
    }
    list(U = U, G = G, W = W)
}

.checkIsoTensorTensor <- function(X, A, rank_individual, rank_isoform,
    initU, initG, initW){
    N <- nrow(X)
    I <- ncol(X)
    C <- ncol(A)
    if(is.null(rank_individual) || is.null(rank_isoform)){
        stop("Please specify both rank_individual and rank_isoform for model='tensor'")
    }
    stopifnot(is.numeric(rank_individual), length(rank_individual) == 1)
    stopifnot(is.numeric(rank_isoform), length(rank_isoform) == 1)
    if(rank_individual != round(rank_individual) || rank_individual < 1){
        stop("rank_individual must be a positive integer")
    }
    if(rank_isoform != round(rank_isoform) || rank_isoform < 1){
        stop("rank_isoform must be a positive integer")
    }
    if(rank_individual > N){
        stop("rank_individual must not exceed the number of individuals nrow(X)")
    }
    if(rank_isoform > I){
        stop("rank_isoform must not exceed the number of isoforms ncol(X)")
    }
    R_N <- rank_individual
    R_I <- rank_isoform
    # Not errors: mathematically allowed but theoretically flagged regimes
    if(R_N >= R_I && C >= 2){
        warning(paste("rank_individual >= rank_isoform:",
            "the individual-specific latent tensor B is generically",
            "non-identifiable under this model even when X is",
            "reconstructed well (Theorem A/A' in notes/tensor-model.md).",
            "Interpret the estimated B with caution."))
    }else if(R_I > R_N && N * (R_I - R_N) < R_N * (C * R_I - R_N)){
        warning(paste0("The number of individuals nrow(X) = ", N,
            " is below the necessary-condition threshold N* = ",
            ceiling(R_N * (C * R_I - R_N) / (R_I - R_N)),
            " for these ranks (Theorem C in notes/tensor-model.md);",
            " the latent tensor B is locally non-identifiable."))
    }
    if(!is.null(initU)){
        stopifnot(is.matrix(initU))
        if(nrow(initU) != N || ncol(initU) != R_N){
            stop("initU must be a nrow(X) x rank_individual matrix")
        }
        if(min(initU) < 0) stop("Please specify the non-negative initU")
    }
    if(!is.null(initG)){
        stopifnot(is.array(initG))
        if(!identical(dim(initG), as.integer(c(R_N, C, R_I)))){
            stop("initG must be a rank_individual x ncol(A) x rank_isoform array")
        }
        if(min(initG) < 0) stop("Please specify the non-negative initG")
    }
    if(!is.null(initW)){
        stopifnot(is.matrix(initW))
        if(nrow(initW) != I || ncol(initW) != R_I){
            stop("initW must be a ncol(X) x rank_isoform matrix")
        }
        if(min(initW) < 0) stop("Please specify the non-negative initW")
    }
}

.eachIsoTensorTensor <- function(X, A, R_N, R_I, initU, initG, initW,
    pseudocount, thr, num.iter, verbose){
    N <- nrow(X)
    I <- ncol(X)
    C <- ncol(A)
    # Strictly positive initialization (avoids zero-locking of MU);
    # user-supplied zeros in init* stay locked at 0 (documented)
    if(is.null(initU)){
        U <- matrix(runif(N * R_N, 0.5, 1.5), nrow = N, ncol = R_N)
    }else{
        U <- initU
    }
    if(is.null(initG)){
        G <- array(runif(R_N * C * R_I, 0.5, 1.5), dim = c(R_N, C, R_I))
    }else{
        G <- initG
    }
    if(is.null(initW)){
        W <- matrix(runif(I * R_I, 0.5, 1.5), nrow = I, ncol = R_I)
    }else{
        W <- initW
    }
    # Optimal global rescaling absorbed into the core
    Xhat <- .XhatTensor(A, U, G, W)
    d <- sum(Xhat * Xhat)
    if(d > 0){
        G <- G * sum(X * Xhat) / d
        Xhat <- Xhat * sum(X * Xhat) / d
    }
    RecError <- c()
    RelChange <- c()
    RecError[1] <- .recError(X, .XhatTensor(A, U, G, W))
    RelChange[1] <- thr * 10
    iter <- 1
    while ((RelChange[iter] > thr) && (iter <= num.iter)) {
        pre_Error <- RecError[iter]
        # A is fixed; each block update is a Lee-Seung NNLS MU step
        U <- .updateU_tensor(X, A, U, G, W, pseudocount)
        G <- .updateG_tensor(X, A, U, G, W, pseudocount)
        W <- .updateW_tensor(X, A, U, G, W, pseudocount)
        # scale absorption into the core (B and Xhat invariant)
        nrm <- .normalizeFactorsTensor(U, G, W)
        U <- nrm$U
        G <- nrm$G
        W <- nrm$W
        iter <- iter + 1
        RecError[iter] <- .recError(X, .XhatTensor(A, U, G, W))
        if(RecError[iter] == 0){
            RelChange[iter] <- 0
        }else{
            RelChange[iter] <- abs(pre_Error - RecError[iter]) / RecError[iter]
        }
        if (verbose) {
            cat(paste0(iter-1, " / ", num.iter,
                " |Previous Error - Error| / Error = ",
                RelChange[iter], "\n"))
        }
        if (is.nan(RelChange[iter])) {
            stop("NaN is generated. Please run again or change the parameters.\n")
        }
    }
    names(RecError) <- c("offset", seq_len(iter-1))
    names(RelChange) <- c("offset", seq_len(iter-1))
    list(U = U, G = G, W = W,
        B = .recTensorPartialTucker(U, G, W),
        Xhat = .XhatTensor(A, U, G, W),
        model = "tensor", algorithm = "Frobenius",
        rank_individual = R_N, rank_isoform = R_I,
        RecError = RecError, RelChange = RelChange,
        NumIter = iter - 1,
        Converged = (RelChange[iter] <= thr))
}
