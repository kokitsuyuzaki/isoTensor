# Phase 2: Tensor model candidates (design note — no implementation yet)

## Problem

Individual-specific cell-type-resolved isoform abundance
$\mathcal{B} \in \mathbb{R}_+^{N \times C \times I}$ with the
observation model

$$
X_{ni} \approx \sum_c A_{nc} B_{nci},
$$

where $A$ ($N \times C$) is fixed. The data has $N I$ entries while a
free $\mathcal{B}$ has $N C I$ parameters, so the problem is
underdetermined by a factor of $C$ for **every** individual: without
pooling across individuals, each slice $B_{n \cdot \cdot}$ is a single
linear equation per isoform ($C$ unknowns, 1 observation). All
identification must therefore come from low-rank structure shared
across individuals, interacting with the variation of $A$ across
individuals. A hard design constraint: the cell type mode is a
biologically explicit axis given by the single-cell reference and is in
principle **not compressed**.

## Candidates

### (a) CP decomposition (rank $R$)

$$
B_{nci} = \sum_{r=1}^{R} U_{nr} V_{cr} W_{ir},
\qquad U \in \mathbb{R}_+^{N \times R},\;
V \in \mathbb{R}_+^{C \times R},\;
W \in \mathbb{R}_+^{I \times R}.
$$

1. **Identifiability**: the strongest of all candidates — CP is
   essentially unique (up to permutation and scaling) under Kruskal's
   condition, and non-negativity further relaxes the requirements.
2. **Interpretability**: each component is a triplet (individual
   loading, cell type loading, isoform loading). However, the cell type
   axis is expressed only through the loadings $V$: a single component
   mixes cell types, and individual-mode and isoform-mode variation are
   forced to covary within one component. This effectively *does*
   compress/entangle the cell type mode, against the design constraint.
3. **Parameters**: $R(N + C + I)$ — the most parsimonious.
4. **MU**: standard non-negative CP MU; with $A$ fixed the model is
   $X_{ni} \approx \sum_r U_{nr} (\sum_c A_{nc} V_{cr}) W_{ir}$ — but
   note the $n$-index appears in both $U$ and $A$, so the effective
   design couples them row-wise (Khatri–Rao-like). MU is derivable
   (each factor enters linearly, so block-wise NNLS + MM applies).
5. **Relation to matrix model**: poor nesting. Recovering Phase 1
   requires $R = C$, $V = I_C$, $U = \mathbf{1}\,$-like structure —
   contrived, and rank is then tied to $C$.

### (b) Full Tucker ($R_N, R_C, R_I$)

$$
B_{nci} = \sum_{p,s,q} U_{np}\, V_{cs}\, G_{psq}\, W_{iq}.
$$

1. **Identifiability**: the weakest — every mode has a
   rotation/monomial indeterminacy absorbable into $\mathcal{G}$;
   non-negativity alone does not restore uniqueness.
2. **Interpretability**: compresses the cell type mode via $V$
   ($C \times R_C$) — directly violates the design constraint. If
   $R_C = C$, $V$ is a square non-negative (near-permutation) matrix
   with no benefit over not having it.
3. **Parameters**: $N R_N + C R_C + I R_I + R_N R_C R_I$.
4. **MU**: available (cf. `nnTensor::NTD`).
5. **Relation to matrix model**: recoverable only through the same
   contortions as (c), with extra indeterminacy on top.

Rejected primarily by criterion 2: the whole point of `isoTensor` is
that the cell type axis is externally anchored by $A$.

### (c) Partial Tucker — cell type mode uncompressed (current favorite)

$$
B_{nci} = \sum_{p=1}^{R_N} \sum_{q=1}^{R_I} U_{np}\, G_{pcq}\, W_{iq},
\qquad
U \in \mathbb{R}_+^{N \times R_N},\;
\mathcal{G} \in \mathbb{R}_+^{R_N \times C \times R_I},\;
W \in \mathbb{R}_+^{I \times R_I}.
$$

1. **Identifiability**: intermediate. Tucker-type indeterminacy remains
   in the two compressed modes (invertible non-negative transforms of
   $U$ and $W$ absorbable into $\mathcal{G}$; effectively monomial —
   permutation + scaling — when factors are sparse). Practical
   mitigations: (i) normalize columns of $U$ and $W$ to unit norm and
   push scale into $\mathcal{G}$; (ii) **fix the first column of $U$ to
   $\mathbf{1}_N$** so that the slice $G_{1 \cdot \cdot}$ is pinned as
   the population component (see criterion 5). Uniqueness should be
   assessed empirically (multi-start stability), as in standard NTD
   practice.
2. **Interpretability**: good. The cell type axis of $\mathcal{B}$ is
   kept intact. $W$ gives $R_I$ non-negative "isoform programs"; $U$
   gives $R_N$ "individual patterns" (e.g., disease groups, genotype
   strata); for each cell type $c$, the matrix $G_{\cdot c \cdot}$
   states how individual patterns use isoform programs *in that cell
   type*. Cell-type-specific individual variation is exactly the object
   of interest and is directly readable from $\mathcal{G}$.
3. **Parameters**: $N R_N + R_N C R_I + I R_I$ — dominated by $I R_I$
   for realistic sizes (e.g., $N=500$, $C=10$, $I=10^5$, $R_N=3$,
   $R_I=50$: $\approx 5.0 \times 10^6$ vs $N I = 5 \times 10^7$
   observations — comfortably overdetermined).
4. **MU**: natural. Each factor block ($U$, $\mathcal{G}$, $W$) enters
   the model linearly, so with the other blocks fixed the subproblem is
   NNLS with a structured design and Lee–Seung MM gives multiplicative
   updates with the usual numerator/denominator split. The fixed $A$
   contracts with $\mathcal{G}$ over $c$ and (row-wise) with $U$ over
   $n$; updates need per-individual weighting by $A$ but stay in
   matrix-product form after mode-wise unfolding. To be derived
   explicitly (gradient split $\nabla^+ - \nabla^-$ per block) before
   implementation.
5. **Relation to matrix model**: exact nesting. With $R_N = 1$,
   $U = \mathbf{1}_N$: $B_{nci} = \sum_q G_{1cq} W_{iq}$, independent
   of $n$; with $R_I = I$, $W = I$ this is literally Phase 1
   ($B = G_{1\cdot\cdot}$). Hence Phase 1's $\hat{B}$ provides a
   warm start: fix $U_{\cdot 1} = \mathbf{1}_N$, initialize
   $G_{1 \cdot \cdot} \hat{=}$ (program-space projection of) $\hat{B}$,
   and let the remaining $R_N - 1$ components capture individual
   deviations.

### (d) Shared population component + individual deviation

$$
B_{nci} = B^0_{ci} + D_{nci}, \qquad D \text{ low-rank},
$$

(or multiplicative variants $B_{nci} = B^0_{ci}\,\exp(d_{nci})$).

1. **Identifiability**: requires a centering constraint (e.g.,
   $\sum_n D_{nci} = 0$) to separate $B^0$ from $D$ — but centering is
   **incompatible with $D \ge 0$**. A non-negative additive $D$ can
   only model upward deviations, silently redefining $B^0$ as a lower
   envelope (noise-sensitive). A signed $D$ leaves the non-negative MU
   framework and needs projection to keep $B \ge 0$.
2. **Interpretability**: the best in principle — "population profile +
   individual deviation" is exactly the biological question, and $B^0$
   is the Phase 1 output.
3. **Parameters**: $CI$ + deviation parameters (comparable to (c)).
4. **MU**: not natural. Additive signed deviation breaks MU;
   multiplicative/log-linear deviation makes the objective
   non-polynomial in the factors.
5. **Relation to matrix model**: perfect ($B^0 = $ Phase 1 solution).

Not viable as-is within the non-negative MU framework, **but** its
interpretation is largely recovered inside candidate (c): with
$U_{\cdot 1} = \mathbf{1}_N$ fixed, $G_{1\cdot\cdot} \times_q W$ is the
shared component and the remaining components are (non-negative,
upward) individual deviations. The "upward-only deviation" limitation
should be stated honestly; if signed deviations become essential, that
is a model change (e.g., HALS with projection, or a hierarchical
shrinkage model à la TCA/bMIND at gene level) to be discussed then.

## Summary and recommendation

| | (a) CP | (b) full Tucker | (c) partial Tucker | (d) shared + deviation |
|---|---|---|---|---|
| Identifiability | **best** (Kruskal) | worst | intermediate, fixable by anchoring | needs centering ⇒ conflicts with $\ge 0$ |
| Cell type axis kept explicit | no (mixed in $V$) | no ($V$ compresses) | **yes** | **yes** |
| Parameters | $R(N{+}C{+}I)$ | $NR_N{+}CR_C{+}IR_I{+}R_NR_CR_I$ | $NR_N{+}R_NCR_I{+}IR_I$ | $CI$ + dev. |
| MU derivable | yes | yes | **yes** | no (signed / non-polynomial) |
| Nests matrix model | contrived | contrived | **exactly** ($R_N{=}1$) | exactly |

**Recommendation: (c) partial non-negative Tucker with the cell type
mode uncompressed**, augmented with two ideas borrowed from the other
candidates:

- from (d): fix $U_{\cdot 1} = \mathbf{1}_N$ and warm-start the
  corresponding core slice from the Phase 1 $\hat{B}$, so the model
  reads "population component + non-negative individual-specific
  components" and Phase 1/Phase 2 results are directly comparable;
- from (a): keep $R_N$, $R_I$ small and check multi-start stability as
  an empirical identifiability diagnostic.

This matches the proposed API
`isoTensor(X, A, model="tensor", rank_individual, rank_isoform)`.

## API note (vs nnTensor)

`nnTensor::NTD` expresses ranks as `rank=c(R1,R2,R3)` plus
`modes=1:3` to select which modes are compressed — partial Tucker is
selected by dropping a mode from `modes`. For `isoTensor` the
cell type mode is *never* compressed by design, so exposing
`modes` would be misleading; the explicit named arguments
`rank_individual` / `rank_isoform` (already reserved in the function
signature) are clearer and are the recommended API. Everything else
follows nnTensor conventions (`initB`, `pseudocount`, `thr`,
`num.iter`, `verbose`, `RecError`/`RelChange` histories).

## Open questions for review

1. Is the upward-only nature of the non-negative individual deviations
   (under the fixed-$U_{\cdot 1}$ anchoring) acceptable for the
   intended biology, or are signed deviations required from the start?
2. Should $W$ (isoform programs) be constrained per gene (isoforms of
   one gene forming a simplex) already in Phase 2, or later?
3. Objective: stay with Frobenius, or move to KL/beta divergence for
   count-like abundances (affects the MU derivation but not the model
   choice)?
