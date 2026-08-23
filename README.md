# isoTensor
Individual-aware cell type-resolved isoform deconvolution using non-negative tensor decomposition.

`isoTensor` takes a bulk isoform abundance matrix `X` (individual x isoform)
and a cell type fraction matrix `A` (individual x cell type, estimated by
external reference-based deconvolution tools such as Bisque, and kept fixed),
and estimates cell-type-resolved isoform abundance.

The current release implements the population-level matrix model
`X ~ A B` (`model="matrix"`), solved by multiplicative updates.
An individual-specific tensor model (`model="tensor"`) is under design;
see `notes/tensor-model.md` and `PLAN.md`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kokitsuyuzaki/isoTensor")
```

## Usage

```r
library("isoTensor")

toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
fit <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000)
fit$B          # cell type x isoform abundance
fit$RecError   # objective history
```

## License

MIT (see `LICENSE.md`).

## Authors

- Koki Tsuyuzaki
