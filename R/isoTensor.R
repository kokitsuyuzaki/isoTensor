isoTensor <- function(X, A, model=c("matrix", "tensor"),
    rank_individual=NULL, rank_isoform=NULL,
    initB=NULL, pseudocount=.Machine$double.eps,
    thr=1e-10, num.iter=100, verbose=FALSE){
    # Argument check
    model <- match.arg(model)
    if(model == "tensor"){
        stop(paste("model='tensor' is not implemented yet.",
            "Phase 2 (individual-specific cell-type-resolved isoform",
            "abundance) is under design. See notes/tensor-model.md",
            "in the GitHub repository."))
    }
    .checkIsoTensor(X, A, initB, pseudocount, thr, num.iter, verbose)
    # Initialization of B
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
