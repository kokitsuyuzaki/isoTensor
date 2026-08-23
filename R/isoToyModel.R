isoToyModel <- function(model=c("matrix"), N=5, C=2, I=6,
    noise=0, seed=NULL){
    model <- match.arg(model)
    stopifnot(is.numeric(N), N >= 1)
    stopifnot(is.numeric(C), C >= 1)
    stopifnot(is.numeric(I), I >= 1)
    stopifnot(is.numeric(noise), noise >= 0)
    if(!is.null(seed)){
        set.seed(seed)
    }
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
