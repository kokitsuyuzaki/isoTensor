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

## Phase 2: Tensor model (design stage — DO NOT implement fitting before review)

Individual-specific cell-type-resolved isoform abundance:

    X_{ni} ≈ Σ_c A_{nc} B_{nci},   B ∈ R_+^{N x C x I}

This is underdetermined (N·C·I parameters vs N·I observations), so a
low-rank structure on B is required. The current candidate is the
**plain partial non-negative Tucker decomposition** with the cell type
mode uncompressed:

    B_nci = Σ_pq U_np G_pcq W_iq

with no anchored/fixed components (the earlier "population component +
non-negative deviation" anchoring was rejected — upward-only deviations;
see "Rejected designs" in `notes/tensor-model.md`). `model="matrix"`
and `model="tensor"` are two different models; they are compared
outside the models (simple and A-weighted averages of B over
individuals vs the matrix-model B).

Before any fitting/MU code: study **recoverability of the latent
tensor by simulation**. The ground-truth generator is implemented as
`isoToyModel(model="tensor", ...)` (partial Tucker factors, Dirichlet
cell fractions with variation controlled by a concentration dial,
optional restriction of individual effects to specific cell types,
noise in X, error in A). Design, grid, metrics, and the
provable-vs-simulation split of identifiability statements:
`notes/tensor-simulation.md` and the "Parameter count and
identifiability" section of `notes/tensor-model.md`.

Status (2026-08-23): identifiability verified
(`notes/tensor-model.md`, `dev/identifiability_check.R`), MU derived
(`notes/tensor-MU.md`) and minimally implemented
(`isoTensor(model="tensor")`): gradient-checked, per-block monotone
(Lee-Seung NNLS reduction), positive control (R_N < R_I) recovers the
latent tensor on toy data, negative control (R_N >= R_I) shows good X
fit with poor/unstable B — as the theory predicts. Non-error warnings
flag the non-identifiable regimes.

Status (2026-08-24): optimizer calibration, the full recoverability
grid and the realistic-noise stress test were run in
`isoTensor-experiments`; the log-reparameterized L-BFGS optimizer
validated there was ported into the package
(`optimizer = c("lbfgs", "hybrid", "mu")`, default "lbfgs";
`notes/tensor-LBFGS.md`, port verified numerically in
`dev/lbfgs_port_check.R`). MU is kept as the provably monotone
baseline. The thr=0 sentinel edge case found by the experiments was
fixed. `algorithm` (objective; currently Frobenius only) and
`optimizer` (method) are kept conceptually separate for future
KL / beta-divergence objectives.

Next gates (not started): rank selection, KL divergence, A
uncertainty / spillover handling, multi-start package API,
Bioconductor submission prep. Real-data analyses stay in
`isoTensor-experiments`.

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
