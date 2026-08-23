# Phase 2: Tensor model (design note — no implementation yet)

Revised 2026-08-23 after review: the anchored
"population component + non-negative deviation" variant was **rejected**
(see "Rejected designs" at the end). The current candidate is the
**plain partial non-negative Tucker decomposition** below. Before any
fitting/MU code is written, the recoverability of the latent tensor
must be studied by simulation (see `notes/tensor-simulation.md`).

## Problem

Individual-specific cell-type-resolved isoform abundance
$\mathcal{B} \in \mathbb{R}_+^{N \times C \times I}$ with the
observation model

$$
X_{ni} \approx \sum_c A_{nc} B_{nci},
$$

where $A \in \mathbb{R}_+^{N \times C}$ is given by an external
deconvolution tool and fixed. The data has $N I$ entries while a free
$\mathcal{B}$ has $N C I$ parameters: for each (individual, isoform)
pair there is 1 observation and $C$ latent values, so without low-rank
structure the problem is plainly underdetermined. All identification
must come from structure shared across individuals interacting with the
variation of $A$ across individuals.

Hard design constraint: the cell type mode is a biologically explicit
axis given by the single-cell reference and is **not compressed**.

## Current candidate: plain partial non-negative Tucker

$$
B_{nci} = \sum_{p=1}^{R_N} \sum_{q=1}^{R_I} U_{np}\, G_{pcq}\, W_{iq},
$$

$$
U \in \mathbb{R}_+^{N \times R_N},\qquad
\mathcal{G} \in \mathbb{R}_+^{R_N \times C \times R_I},\qquad
W \in \mathbb{R}_+^{I \times R_I}.
$$

Equivalently $\mathcal{B} = \mathcal{G} \times_1 U \times_3 W$ (mode-2,
the cell type mode, is left uncompressed). No factor is anchored or
fixed; there is no built-in "population component". `model="matrix"`
and `model="tensor"` are treated as two different models with different
assumptions — algebraic nesting of the matrix model inside the tensor
model is **not** a design requirement.

Interpretation: $W$ gives $R_I$ non-negative isoform programs, $U$
gives $R_N$ individual patterns, and for each cell type $c$ the slice
$G_{\cdot c \cdot}$ states how individual patterns use isoform programs
in that cell type.

## Relationship to the matrix model

The matrix model ($X \approx A B$, $B \in \mathbb{R}_+^{C \times I}$)
estimates a population-level profile; the tensor model estimates
$\mathcal{B} \in \mathbb{R}_+^{N \times C \times I}$. They are compared
**outside** the models, e.g. by

$$
\bar{B}^{\mathrm{tensor}}_{ci} = \frac{1}{N} \sum_n B^{\mathrm{tensor}}_{nci}
\qquad \text{(simple average)}
$$

versus $\hat{B}^{\mathrm{matrix}}_{ci}$. A weighted summary should also
be considered, because the matrix model's least squares estimand is an
$A$-weighted object (the normal equations
$A^\top A\, \hat{B} = A^\top X$ weight individual $n$ by its fractions
$a_n$): a natural counterpart is the fraction-weighted average

$$
\bar{B}^{\mathrm{tensor,w}}_{ci}
= \frac{\sum_n A_{nc}\, B^{\mathrm{tensor}}_{nci}}{\sum_n A_{nc}} .
$$

When individual effects are absent ($B_{nci} = B_{ci}$), both summaries
coincide with $B_{ci}$ and the matrix model is sufficient (this is
Simulation 1 in `notes/tensor-simulation.md`; verified numerically in
the unit tests: the matrix model recovers the shared slice exactly on
noise-free tensor-generated data with `individual_effect=FALSE`).

## Parameter count and identifiability

### Counting

Raw factor parameters:

$$
P = N R_N + R_N C R_I + I R_I
\qquad \text{vs} \qquad
N I \ \text{observations}.
$$

$P < NI$ is easy to satisfy (e.g., $N=500$, $C=10$, $I=10^5$, $R_N=3$,
$R_I=50$: $P \approx 5.0\times10^6$ vs $NI = 5\times10^7$), but **this
is only a necessary-type condition, not sufficient** for
identifiability of $\mathcal{B}$.

Two indeterminacies must be distinguished:

- **Factor indeterminacy** (harmless for our target): for any
  invertible $S \in \mathrm{GL}(R_N)$, $T \in \mathrm{GL}(R_I)$,
  $(U S,\ \mathcal{G} \times_1 S^{-1} \times_3 T^{-1},\ W T)$ gives the
  same $\mathcal{B}$. This is why factor-wise element comparison is not
  a primary evaluation metric; the evaluation target is the
  reconstructed $\mathcal{B}$ itself.
- **$\mathcal{B}$-indeterminacy** (fatal): different $\mathcal{B}$
  within the model class producing the same $X$. This is what the
  simulation study must characterize.

### Structure of the observation map

Contracting over $c$ first,

$$
X_{ni} = \sum_q H_{nq} W_{iq},
\qquad
H_{nq} = \sum_p U_{np} \big[\mathcal{G} \times_2 a_n\big]_{pq},
$$

where $a_n$ is the $n$-th row of $A$ and
$\mathcal{G} \times_2 a_n \in \mathbb{R}^{R_N \times R_I}$. So the data
factorizes as $X = H W^\top$ (an NMF of rank $R_I$), and the
cell type information enters only through the per-individual
contractions $\mathcal{G} \times_2 a_n$.

### What can be stated without simulation

These follow from the linear-algebraic structure (genericity assumed
where noted); they are claims about what is *impossible*, so they hold
regardless of the fitting algorithm:

1. **No variation in $A$ implies non-recoverability (provable).** If
   all rows are equal ($a_n = a$, $\mathrm{rank}(A)=1$) and $C > 1$,
   only the contraction $\sum_c a_c B_{nci}$ is observable. Explicitly:
   for any $\mathcal{E}$ with $\sum_c a_c E_{pcq} = 0$ (a
   $(C-1) R_N R_I$-dimensional space), replacing $\mathcal{G}$ by
   $\mathcal{G} + \mathcal{E}$ changes $\mathcal{B}$ but not $X$; for
   strictly positive $\mathcal{G}$ a small such $\mathcal{E}$ keeps
   non-negativity. The cell type axis is completely unidentified.
2. **$R_N \ge R_I$ breaks sign-unconstrained identifiability
   (generic).** Given the factorization $X = H W^\top$, recovery of
   $(U, \mathcal{G})$ from $H$ requires solving
   $H_{n\cdot} = u_n^\top (\mathcal{G} \times_2 a_n)$ per individual.
   If $R_N \ge R_I$, then for a *generic alternative* core
   $\mathcal{G}'$ the matrix $\mathcal{G}' \times_2 a_n$ has full row
   space $\mathbb{R}^{R_I}$, so some $u'_n$ reproduces $H$ exactly for
   every $n$ — the latent tensor is not determined by $X$ if signs are
   unconstrained. Whether non-negativity (inequality constraints)
   rescues this regime is a **simulation question**. Consequence for
   the study design: include $R_N < R_I$ settings as the primary
   regime, and treat $R_N = R_I$ as a probe of the effect of
   non-negativity.
3. **Dimension-count threshold on $N$ (heuristic necessary
   condition).** In the same second stage, unknowns modulo
   $\mathrm{GL}(R_N)$ number $N R_N + R_N C R_I - R_N^2$ against
   $N R_I$ equations, giving

   $$
   N (R_I - R_N) \ \gtrsim\ R_N (C R_I - R_N).
   $$

   E.g., $R_N=2$, $R_I=3$, $C=3$: $N \gtrsim 14$. This motivates the
   small-grid choice $N = 10\text{–}30$: the transition should be
   visible there. This is a counting heuristic (ignores
   non-negativity, assumes the $W$-stage resolved), not a theorem.
4. **The $W$-stage is an NMF.** Uniqueness of $X = H W^\top$ up to
   monomial transforms requires standard NMF-type conditions
   (separability / sufficiently-scattered factors). These conditions
   are not verifiable a priori on real data; in simulation they
   translate into "sparser / more distinct isoform programs are easier"
   and should be varied empirically.

### What only simulation can tell us

- Whether and when **non-negativity** restores recoverability in
  regimes that are unidentifiable without sign constraints (including
  $R_N = R_I$).
- The **quantitative** dependence of recovery error on the variation of
  $A$ (Dirichlet concentration), $\mathrm{rank}(A)$ /
  condition number / column correlations, rare cell type frequency,
  noise in $X$, and error in $A$.
- Robustness to **rank misspecification** and to model misspecification
  (true $\mathcal{B}$ not exactly partial Tucker).

No claims beyond the four numbered statements above should be made
without either a proof or supporting simulation.

## Comparison of candidate structures (kept for the record)

| | (a) CP | (b) full Tucker | (c) plain partial Tucker | (d) shared + deviation |
|---|---|---|---|---|
| Identifiability | best (Kruskal) | worst | intermediate; see above | needs centering ⇒ conflicts with $\ge 0$ |
| Cell type axis kept explicit | no (mixed in $V$) | no ($V$ compresses) | **yes** | yes |
| Parameters | $R(N{+}C{+}I)$ | $NR_N{+}CR_C{+}IR_I{+}R_NR_CR_I$ | $NR_N{+}R_NCR_I{+}IR_I$ | $CI$ + dev. |
| MU derivable | yes | yes | **yes** (to be derived after the simulation study) | no (signed / non-polynomial) |

- **(a) CP**: strongest uniqueness theory, but a single component mixes
  cell types through its loading vector, effectively entangling the
  cell type mode — against the design constraint.
- **(b) full Tucker**: compresses the cell type mode; rejected outright.
- **(c) plain partial Tucker**: the current candidate (above).
- **(d) shared + deviation** ($B_{nci} = B^0_{ci} + D_{nci}$):
  separating $B^0$ from $D$ requires a centering constraint
  incompatible with $D \ge 0$; a non-negative $D$ only models upward
  deviations; a signed $D$ leaves the non-negative MU framework.
  Rejected as a model; the population-level comparison is instead done
  outside the model (see "Relationship to the matrix model").

## Rejected designs

**Anchored partial Tucker (rejected 2026-08-23).** An earlier draft of
this note recommended fixing $U_{\cdot 1} = \mathbf{1}_N$, warm-starting
the core slice $G_{1\cdot\cdot}$ from the Phase 1 $\hat{B}$, and
interpreting components $p \ge 2$ as non-negative individual
deviations, so that the matrix model was exactly nested in the tensor
model. This is **not adopted**:

- with the remaining components constrained non-negative, individual
  effects are limited to *upward* deviations from the population
  component, which is an unacceptable modeling bias;
- exact algebraic nesting of the matrix model is not a design
  requirement; the two models are compared externally instead.

(The *generator* in `isoToyModel(model="tensor", effect_cell_types=...)`
does use a fixed shared component internally to construct ground truth
with cell-type-restricted individual effects — that is a simulation
construction, not a constraint of the fitting model.)

## API note (vs nnTensor)

`nnTensor::NTD` expresses ranks as `rank=c(R1,R2,R3)` plus `modes=` to
select which modes are compressed. For `isoTensor` the cell type mode
is *never* compressed by design, so exposing `modes` would be
misleading; the explicit named arguments `rank_individual` /
`rank_isoform` (already reserved in the `isoTensor()` signature and
used by `isoToyModel(model="tensor")`) are the recommended API.
Everything else follows nnTensor conventions (`initB`, `pseudocount`,
`thr`, `num.iter`, `verbose`, `RecError`/`RelChange` histories).

## Open questions for review

1. Simulation design in `notes/tensor-simulation.md`: is the grid and
   the set of recorded diagnostics sufficient before MU derivation?
2. Objective for Phase 2: stay with Frobenius, or move to KL/beta
   divergence for count-like abundances (affects the MU derivation but
   not the model choice)?
3. Should $W$ (isoform programs) be constrained per gene (isoforms of
   one gene forming a simplex) already in Phase 2, or later?
