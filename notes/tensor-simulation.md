# Phase 2: Simulation framework for tensor recoverability

Purpose: before implementing any fitting algorithm (no MU yet), decide
by simulation **under which conditions the individual × cell type ×
isoform latent tensor $\mathcal{B}$ is recoverable** from
$X_{ni} = \sum_c A_{nc} B_{nci}$ under the plain partial Tucker
constraint (see `notes/tensor-model.md` for the model and for what is
provable without simulation).

Division of labor:

- **This repository**: the ground-truth generator
  (`isoToyModel(model="tensor", ...)`) and its unit tests, plus this
  design note.
- **`isoTensor-experiments`**: actually running simulation grids,
  benchmarking, and figures.

## 1. Ground-truth generator

Implemented as `isoToyModel(model="tensor", ...)` — one exported
function shared with Phase 1 rather than a new API
(`model=` selects the generator, mirroring `isoTensor(model=)`).

Generation order:

1. Factors $U_{\mathrm{true}} \in \mathbb{R}_+^{N \times R_N}$
   (entries $\mathrm{Unif}(0.5, 1.5)$),
   $\mathcal{G}_{\mathrm{true}} \in \mathbb{R}_+^{R_N \times C \times R_I}$
   ($\mathrm{Unif}(0.5, 1.5)$),
   $W_{\mathrm{true}} \in \mathbb{R}_+^{I \times R_I}$
   ($\mathrm{Unif}(0, 1)$).
2. $\mathcal{B}_{\mathrm{true}} = \mathcal{G}_{\mathrm{true}} \times_1
   U_{\mathrm{true}} \times_3 W_{\mathrm{true}}$
   (cell type mode uncompressed).
3. $A_{\mathrm{true}}$: rows drawn i.i.d. from
   $\mathrm{Dirichlet}(\alpha_0\, m)$ (see §2).
4. $X_{ni} = \sum_c A_{nc} B_{nci}$ (noise-free), then optionally
   $X \leftarrow \max(X + \mathcal{N}(0, \texttt{noise}^2), 0)$.
5. Optionally a *perturbed* fraction matrix `A_obs`
   ($A + \mathcal{N}(0, \texttt{noise\_A}^2)$, clipped at 0,
   rows re-normalized) to mimic deconvolution error: $X$ is always
   generated from the true $A$, while the fitting step is fed `A_obs`.

Switches for the baseline comparisons (§5):

- `individual_effect=FALSE`: all rows of $U$ are made identical, so
  $B_{nci} = B_{ci}$ (Simulation 1).
- `effect_cell_types=c(...)`: component 1 of $U$ is a shared component
  and the core slices of components $p \ge 2$ are zeroed outside the
  listed cell types, so individual variation is restricted to specific
  cell types (Simulation 2). Restricting to specific isoform programs
  is achieved the same way on the $q$ axis if needed (not yet exposed
  as an argument).

Everything (factors, $\mathcal{B}_{\mathrm{true}}$, $A$, `A_obs`,
diagnostics) is returned so that any future fitting result can be
scored against the ground truth.

## 2. Variation of A: generation and quantification

Hypothesis to test: if cell fractions barely vary across individuals,
the information to separate the cell type direction is missing and
recoverability degrades; sufficient variation plus shared low-rank
structure should improve it. (The extreme case is provable: identical
rows ⇒ the cell type axis is completely unidentified; see
`notes/tensor-model.md` statement 1.)

Generation: rows of $A$ are $\mathrm{Dirichlet}(\alpha_0\, m)$ where
$m$ (`A_mean`) is the mean composition and $\alpha_0$ the total
concentration; then $\mathrm{Var}(A_{nc}) = m_c(1-m_c)/(\alpha_0+1)$,
so $\alpha_0$ is a single dial for between-individual variation:

| level (`A_variation`) | $\alpha_0$ | intent |
|---|---|---|
| `"none"` | $\infty$ (identical rows) | provably unrecoverable control |
| `"low"` | 200 | near-degenerate regime |
| `"moderate"` | 20 | realistic bulk cohorts |
| `"high"` | 2 | strongly varying composition |

A numeric `A_variation` is passed through as $\alpha_0$ directly, so
grids can sweep it continuously. Rare cell types are created via a
skewed `A_mean` (e.g., `c(0.9, 0.05, 0.05)`).

Quantification (recorded per run in `A_summary`):

- $\mathrm{rank}(A)$ (numerical, SVD-based)
- condition number $\sigma_{\max}/\sigma_{\min}$
- pairwise correlation matrix of the fraction columns
- mean fractions and minimum mean fraction (rare cell type indicator)
- total variance $\sum_c \mathrm{Var}_n(A_{nc})$ (overall variation
  dial, monotone in $1/\alpha_0$)

## 3. Simulation grid

Axes to vary (full grid lives in `isoTensor-experiments`):

- $N$ (individuals), $C$ (cell types), $I$ (isoforms)
- $R_N$, $R_I$ (true ranks; later also misspecified fitting ranks)
- variation of $A$ ($\alpha_0$), rare cell types (`A_mean`)
- noise in $X$ (`noise`), error in $A$ (`noise_A`)
- similarity between cell-type-specific profiles (via the spread of
  $\mathcal{G}$ slices across $c$; not yet an argument — to be added
  when the grid runs start)

**Start small.** First stage (eyeball + numerical checks, runs in
seconds):

- $N = 10\text{–}30$, $C = 2\text{–}3$, $I = 10\text{–}30$,
  $R_N = 2$, $R_I = 2$, noise-free.

Two flags from the identifiability analysis
(`notes/tensor-model.md` statements 2–3):

- $R_N = R_I$ is a regime where sign-unconstrained recovery provably
  fails; keep it in the grid *as the probe of what non-negativity
  buys*, but add $R_N < R_I$ settings (e.g., $R_N=2$, $R_I=3, 4$) as
  the primary recoverable regime.
- The counting heuristic $N(R_I - R_N) \gtrsim R_N(C R_I - R_N)$
  predicts a transition around $N \approx 14$ for
  $(R_N, R_I, C) = (2, 3, 3)$ — the $N = 10\text{–}30$ range straddles
  it deliberately; check whether the empirical transition tracks the
  prediction.

## 4. Evaluation metrics (for the future fitting step)

Primary target: the reconstructed latent tensor $\mathcal{B}$, **not**
the factors ($U, \mathcal{G}, W$ carry Tucker indeterminacies —
invertible/monomial transforms leave $\mathcal{B}$ unchanged — so
element-wise factor agreement must not be a primary metric).

1. **Observed-data reconstruction**: $\|X - \hat{X}\|_F$ (and relative
   version). Note: small $X$-residual with large $\mathcal{B}$-error is
   exactly the non-identifiability signature we are looking for —
   always report both.
2. **Latent tensor recovery**:
   $\|\mathcal{B}_{\mathrm{true}} - \hat{\mathcal{B}}\|_F$ and the
   relative error
   $\|\mathcal{B}_{\mathrm{true}} - \hat{\mathcal{B}}\|_F /
   \|\mathcal{B}_{\mathrm{true}}\|_F$.
3. **Cell-type-specific recovery**: per cell type $c$, relative
   Frobenius error and correlation between
   $B_{\mathrm{true},\cdot c \cdot}$ and $\hat{B}_{\cdot c \cdot}$
   (expected to degrade first for rare / collinear cell types).
4. **Individual-specific recovery**: per individual $n$,
   correlation / error of the latent profile $B_{n \cdot \cdot}$;
   also the recovery of the *deviation from the population mean*
   $B_{nci} - \bar{B}_{ci}$, which is the quantity that distinguishes
   the tensor model from the matrix model.
5. **Population-level agreement**: $\bar{B}^{\mathrm{tensor}}$ (simple
   mean over $n$) and the fraction-weighted
   $\bar{B}^{\mathrm{tensor,w}}_{ci} = \sum_n A_{nc} B_{nci} / \sum_n
   A_{nc}$ versus $\hat{B}^{\mathrm{matrix}}$ on the same data (the
   weighted version matches the matrix model's least squares estimand;
   see `notes/tensor-model.md`).
6. **Multi-start stability** (empirical identifiability diagnostic):
   dispersion of $\hat{\mathcal{B}}$ across random restarts; high
   dispersion at near-zero $X$-residual indicates a flat,
   non-identified regime.

## 5. Baseline simulations against the matrix model

The matrix model (Phase 1, already implemented) runs on every synthetic
dataset as the baseline.

**Simulation 1 — no individual effect** ($B_{nci} = B_{ci}$;
`individual_effect=FALSE`). Expectation: the matrix model is
sufficient; the tensor model is unnecessary (and must not do worse).
Already verified for the matrix side in the unit tests: with
`A_variation="moderate"` (so $\mathrm{rank}(A) = C$) the matrix model
recovers the shared slice to $<10^{-6}$ on noise-free data.

**Simulation 2 — individual effect present**, restrictable to specific
cell types (`effect_cell_types`) and, later, specific isoform programs.
Expectation: the matrix model collapses everything to a population
average (its estimand is the $A$-weighted mean); the tensor model may
recover the individual variation — this is precisely the recoverability
question, scored with metrics 2–4.

**Simulation 3 — phenotype-associated effect (future).** Two groups of
individuals differing only in a specific isoform program of a specific
cell type; evaluate whether the group difference is recovered from
$\hat{\mathcal{B}}$ (or the individual representation $\hat{U}$). To be
designed in `isoTensor-experiments` once fitting exists; the generator
will need a `group` argument (not yet implemented).

## 6. Per-run record

Each simulation run should record: all generator arguments (incl.
seed), `A_summary` (§2), the metrics of §4, iteration counts, and the
number of restarts — so that recoverability maps (e.g., error vs
$\alpha_0$ × $N$) can be drawn without re-running.

## 7. Gate

MU derivation and fitting implementation start only after this
simulation design (and the revised model note) is reviewed and the
small-grid recoverability results justify the model choice.
