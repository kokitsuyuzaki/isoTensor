# isoTensor
Individual-aware cell type-resolved isoform deconvolution using non-negative tensor decomposition.

`isoTensor` takes a bulk isoform abundance matrix `X` (individual x isoform)
and a cell type fraction matrix `A` (individual x cell type, estimated by
external reference-based deconvolution tools such as Bisque, and kept fixed),
and estimates cell-type-resolved isoform abundance.

Two models are implemented:

- `model="matrix"`: the population-level model `X ~ A B`
  (`B`: cell type x isoform), solved by multiplicative updates.
- `model="tensor"`: an individual-specific latent tensor
  `B` (individual x cell type x isoform) with a partial non-negative
  Tucker structure `B_nci = sum_pq U_np G_pcq W_iq` (the cell type
  mode is not compressed); see `notes/tensor-model.md` and `PLAN.md`.

For the tensor model, three optimizers minimize the same Frobenius
objective (`optimizer=`): `"lbfgs"` (default) is a
log-reparameterized L-BFGS with analytical gradients — in synthetic
benchmarks it is several times faster than the multiplicative updates
and recovers the latent tensor far more accurately; `"mu"` is the
multiplicative-update baseline whose objective is provably
non-increasing (Lee-Seung); `"hybrid"` runs MU first and then refines
with L-BFGS. See `notes/tensor-MU.md` and `notes/tensor-LBFGS.md`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("kokitsuyuzaki/isoTensor")
```

## Usage

```r
library("isoTensor")

# Population-level matrix model
toy <- isoToyModel(N=5, C=2, I=6, seed=1234)
fit <- isoTensor(toy$X, toy$A, model="matrix", num.iter=1000)
fit$B          # cell type x isoform abundance
fit$RecError   # objective history

# Individual-specific tensor model
tt <- isoToyModel(model="tensor", N=20, C=2, I=20,
    rank_individual=1, rank_isoform=2, A_variation="high", seed=10)
fit <- isoTensor(
    tt$X,
    tt$A,
    model = "tensor",
    rank_individual = 1,
    rank_isoform = 2,
    optimizer = "lbfgs")
fit$B          # individual x cell type x isoform latent tensor
fit$objective  # final 1/2 ||X - Xhat||_F^2
```

## Limitations

- The tensor mode rests on strong structural assumptions (low-rank partial
  Tucker with the cell type mode uncompressed, fixed `A`).
- For `rank_individual >= rank_isoform` the latent tensor is generically
  non-identifiable under the current dense model (a warning is raised);
  `rank_individual < rank_isoform` avoids that proved condition but does
  not by itself guarantee identifiability.
- Good reconstruction of `X` does not guarantee recovery of the latent
  tensor `B` — evaluate both.
- Individual-level tensor estimates can be highly sensitive to isoform
  measurement noise and to a poorly conditioned or low-variation `A`
  (rare cell types are especially unstable); see the recoverability and
  noise studies in the sibling `isoTensor-experiments` repository
  (`reports/FINAL_SUMMARY.md`).
- The matrix mode is the safer population-level option and is equivalent
  to per-isoform NNLS.

## License

MIT (see `LICENSE.md`).

## Authors

- Koki Tsuyuzaki
