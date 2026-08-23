\name{isoToyModel}
\docType{methods}
\alias{isoToyModel}

\title{
Toy Model Data for isoTensor
}
\description{
isoToyModel generates a small synthetic dataset for isoTensor. A ground
truth cell type fraction matrix A (individual x cell type, rows sum to
1) and a ground truth cell-type-by-isoform matrix B are generated, and
the bulk isoform abundance matrix X = A B is returned, optionally with
additive Gaussian noise truncated at 0.
}
\usage{
isoToyModel(model=c("matrix"), N=5, C=2, I=6, noise=0, seed=NULL)
}

\arguments{
\item{model}{
Only "matrix" is available for now (Default: "matrix").
}
\item{N}{
The number of individuals (Default: 5).
}
\item{C}{
The number of cell types (Default: 2).
}
\item{I}{
The number of isoforms (Default: 6).
}
\item{noise}{
The standard deviation of the additive Gaussian noise added to X.
Negative values after adding noise are truncated to 0 (Default: 0).
}
\item{seed}{
If not NULL, set.seed(seed) is called before the data generation for
reproducibility (Default: NULL).
}
}

\value{
A list containing X (N x I), A (N x C), B (C x I), and the
specified N, C, I, and noise.
}

\examples{
toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
dim(toy$X)
}

\author{Koki Tsuyuzaki}

\seealso{\code{\link{isoTensor}}}
\keyword{methods}
