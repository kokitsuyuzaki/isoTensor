# isoTensor development plan

## Scope

`isoTensor` estimates cell-type-resolved isoform abundance from bulk RNA-seq.

Inputs:

- `X` (N x I, non-negative): bulk isoform abundance (individual x isoform)
- `A` (N x C, non-negative): cell type fractions (individual x cell type),
  estimated by external reference-based deconvolution tools (e.g., Bisque)

Design principles:

- Estimating cell type fractions is **out of scope**. `A` is always fixed
  and never updated.
- The cell type mode is a biologically explicit axis given by the
  single-cell reference; it is in principle **not compressed** in any model.
- This repository contains the R/Bioconductor package only. Real-data
  analyses (GTEx / ROSMAP), benchmarking, and paper figure reproduction
  live in the separate `isoTensor-experiments` repository.

## Phase 1: Matrix model (implemented)

Population-level isoform deconvolution:

    X ≈ A B,   B ∈ R_+^{C x I}
    min_{B >= 0} (1/2) ||X - A B||_F^2

solved by the multiplicative update

    B <- B ⊙ (A^T X) / (A^T A B + ε)

See `notes/matrix-model.md` for the derivation and convergence properties.

Status: implemented in `R/isoTensor.R` (`model="matrix"`), tested in
`tests/testthat/test_isoTensor.R`, documented in `vignettes/isoTensor.Rmd`.

Possible follow-ups (not yet scheduled):

- KL / beta-divergence objectives (RNA-seq counts are heteroscedastic)
- Row masks / weights (e.g., per-individual library size weights)
- Isoform-group (per-gene) normalization or simplex constraints on B
- Rank-free diagnostics: residual analysis per isoform / individual

## Phase 2: Tensor model (design stage — DO NOT implement before review)

Individual-specific cell-type-resolved isoform abundance:

    X_{ni} ≈ Σ_c A_{nc} B_{nci},   B ∈ R_+^{N x C x I}

This is underdetermined (N·C·I parameters vs N·I observations), so a
low-rank structure on B is required. Candidate structures (CP, full
Tucker, partial Tucker, shared + deviation) are compared in
`notes/tensor-model.md`. The current recommendation is a **partial
non-negative Tucker decomposition that leaves the cell type mode
uncompressed**, with the matrix model as its exact special case.

Gate: the model choice must be reviewed and approved before any Phase 2
code is written.

## Phase 3 and beyond (tentative)

- Bioconductor submission (BiocCheck, S4/SummarizedExperiment input
  support if requested by reviewers)
- Statistical extensions: uncertainty quantification, cross-validated
  rank selection (Gabriel-style holdout as in nnTensor)

## Engineering conventions

- Follow nnTensor coding style: plain-function API, 4-space indent,
  dot-prefixed internal helpers, hand-written `man/*.Rd`, testthat with
  explicit `test_file()` listing, Rmd vignettes.
- Return value style follows nnTensor: a plain list with factor
  matrices plus `RecError` / `RelChange` histories (element 1 is the
  "offset" = value at initialization).
- Local test environment: conda env `r_4.3`
  (`/home/koki/anaconda3/envs/r_4.3/bin/R`, R 4.3.3). The `Depends:
  R (>= 4.3.0)` line must be raised to the Bioc-devel requirement at
  submission time.
