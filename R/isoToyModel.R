isoToyModel <- function(model=c("matrix", "tensor"), N=5, C=2, I=6,
    noise=0, seed=NULL,
    rank_individual=2, rank_isoform=2,
    individual_effect=TRUE, effect_cell_types=NULL,
    A_variation="moderate", A_mean=NULL, noise_A=0){
    model <- match.arg(model)
    stopifnot(is.numeric(N), N >= 1)
    stopifnot(is.numeric(C), C >= 1)
    stopifnot(is.numeric(I), I >= 1)
    stopifnot(is.numeric(noise), noise >= 0)
    if(!is.null(seed)){
        set.seed(seed)
    }
    if(model == "matrix"){
        .isoToyModelMatrix(N, C, I, noise)
    }else{
        .isoToyModelTensor(N, C, I, rank_individual, rank_isoform,
            individual_effect, effect_cell_types,
            A_variation, A_mean, noise, noise_A)
    }
}

.isoToyModelMatrix <- function(N, C, I, noise){
    # A: individual x cell type fractions (rows sum to 1)
    A <- matrix(runif(N * C, min = 0.1, max = 1), nrow = N, ncol = C)
    A <- A / rowSums(A)
    # B: cell type x isoform abundance with cell-type-specific isoforms
    B <- matrix(runif(C * I, min = 0.5, max = 2), nrow = C, ncol = I)
    for(c in seq_len(C)){
        specific <- which(seq_len(I) %% C == (c - 1))
        B[c, specific] <- B[c, specific] * 5
    }
    X <- A %*% B
    if(noise > 0){
        X <- X + matrix(rnorm(N * I, mean = 0, sd = noise),
            nrow = N, ncol = I)
        X[X < 0] <- 0
    }
    list(X = X, A = A, B = B, N = N, C = C, I = I, noise = noise)
}

.isoToyModelTensor <- function(N, C, I, R_N, R_I,
    individual_effect, effect_cell_types,
    A_variation, A_mean, noise, noise_A){
    stopifnot(is.numeric(R_N), R_N >= 1)
    stopifnot(is.numeric(R_I), R_I >= 1)
    stopifnot(is.logical(individual_effect))
    stopifnot(is.numeric(noise_A), noise_A >= 0)
    if(!is.null(effect_cell_types)){
        stopifnot(all(effect_cell_types %in% seq_len(C)))
        if(R_N < 2){
            stop("effect_cell_types requires rank_individual >= 2 (one shared component + individual-varying components)")
        }
    }
    # A: individual x cell type fractions (Dirichlet rows)
    A <- .generateA(N, C, A_variation, A_mean)
    # Ground truth factors of the partial Tucker model
    # B_nci = sum_pq U_np G_pcq W_iq (cell type mode is not compressed)
    U <- matrix(runif(N * R_N, min = 0.5, max = 1.5), nrow = N, ncol = R_N)
    G <- array(runif(R_N * C * R_I, min = 0.5, max = 1.5),
        dim = c(R_N, C, R_I))
    W <- matrix(runif(I * R_I, min = 0, max = 1), nrow = I, ncol = R_I)
    if(!individual_effect){
        # Identical rows of U so that B_nci does not depend on n
        U <- matrix(rep(U[1, ], each = N), nrow = N, ncol = R_N)
    }else if(!is.null(effect_cell_types)){
        # Component 1 is shared across individuals; components >= 2
        # carry individual variation restricted to effect_cell_types
        U[, 1] <- 1
        others <- setdiff(seq_len(C), effect_cell_types)
        G[seq_len(R_N)[-1], others, ] <- 0
    }
    B <- .recTensorPartialTucker(U, G, W)
    # X_ni = sum_c A_nc B_nci
    X <- matrix(0, nrow = N, ncol = I)
    for(c in seq_len(C)){
        X <- X + A[, c] * B[, c, ]
    }
    if(noise > 0){
        X <- X + matrix(rnorm(N * I, mean = 0, sd = noise),
            nrow = N, ncol = I)
        X[X < 0] <- 0
    }
    # Observed (misspecified) cell fractions handed to the fitting step
    A_obs <- A
    if(noise_A > 0){
        A_obs <- A + matrix(rnorm(N * C, mean = 0, sd = noise_A),
            nrow = N, ncol = C)
        A_obs[A_obs < 0] <- 0
        A_obs <- A_obs / rowSums(A_obs)
    }
    list(X = X, A = A, A_obs = A_obs, B = B, U = U, G = G, W = W,
        N = N, C = C, I = I,
        rank_individual = R_N, rank_isoform = R_I,
        individual_effect = individual_effect,
        effect_cell_types = effect_cell_types,
        A_variation = A_variation, noise = noise, noise_A = noise_A,
        A_summary = .summarizeA(A))
}

# Total Dirichlet concentration alpha0 controlling the between-individual
# variation of cell fractions: Var(A_nc) = m_c (1 - m_c) / (alpha0 + 1)
.A_concentration <- function(A_variation){
    if(is.numeric(A_variation)){
        stopifnot(A_variation > 0)
        return(A_variation)
    }
    switch(A_variation,
        "none" = Inf,
        "low" = 200,
        "moderate" = 20,
        "high" = 2,
        stop("A_variation must be 'none', 'low', 'moderate', 'high', or a positive number (total Dirichlet concentration)"))
}

.generateA <- function(N, C, A_variation, A_mean){
    if(is.null(A_mean)){
        A_mean <- rep(1 / C, C)
    }
    stopifnot(length(A_mean) == C)
    stopifnot(all(A_mean > 0))
    A_mean <- A_mean / sum(A_mean)
    alpha0 <- .A_concentration(A_variation)
    if(is.infinite(alpha0)){
        A <- matrix(A_mean, nrow = N, ncol = C, byrow = TRUE)
    }else{
        A <- t(vapply(seq_len(N), function(n){
            g <- rgamma(C, shape = alpha0 * A_mean)
            if(sum(g) == 0){
                A_mean
            }else{
                g / sum(g)
            }
        }, numeric(C)))
    }
    A
}

# B_nci = sum_pq U_np G_pcq W_iq
.recTensorPartialTucker <- function(U, G, W){
    N <- nrow(U)
    R_N <- ncol(U)
    C <- dim(G)[2]
    R_I <- dim(G)[3]
    I <- nrow(W)
    B <- array(0, dim = c(N, C, I))
    for(c in seq_len(C)){
        Gc <- matrix(G[, c, ], nrow = R_N, ncol = R_I)
        B[, c, ] <- U %*% Gc %*% t(W)
    }
    B
}

# Diagnostics of A related to recoverability (recorded per simulation run)
.summarizeA <- function(A){
    d <- svd(A, nu = 0, nv = 0)$d
    tol <- max(dim(A)) * .Machine$double.eps * max(d)
    rankA <- sum(d > tol)
    condA <- if(min(d) > 0) max(d) / min(d) else Inf
    corA <- suppressWarnings(cor(A))
    mean_fraction <- colMeans(A)
    total_variation <- sum(apply(A, 2, var))
    list(rank = rankA, condition = condA, cor = corA,
        mean_fraction = mean_fraction,
        min_mean_fraction = min(mean_fraction),
        total_variation = total_variation)
}
