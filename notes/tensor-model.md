# Phase 2: Tensor model (design note — no implementation yet)

Revised 2026-08-23 (2nd revision): the identifiability section was
re-derived rigorously and verified numerically
(`dev/identifiability_check.R`); every claim is now classified as
PROVED / NECESSARY CONDITION / LOCAL-GENERIC EVIDENCE / EMPIRICAL
HYPOTHESIS / UNSUPPORTED. Headline change: the earlier
"$R_N \ge R_I$ breaks sign-unconstrained identifiability" claim is
**upgraded to PROVED by explicit construction, and strengthened**: it
also holds *under non-negativity* whenever the true factors are
strictly positive (interior). The simulation grid is therefore NOT
restricted to $R_N < R_I$; all three regimes are tested (see
`notes/tensor-simulation.md`).

## Problem

Individual-specific cell-type-resolved isoform abundance
$\mathcal{B} \in \mathbb{R}_+^{N \times C \times I}$ with observation

$$
X_{ni} \approx \sum_c A_{nc} B_{nci},
$$

$A \in \mathbb{R}_+^{N \times C}$ known and fixed (rows are cell
fractions, non-zero). Hard design constraint: the cell type mode is
not compressed.

## Current candidate: plain partial non-negative Tucker

$$
B_{nci} = \sum_{p=1}^{R_N} \sum_{q=1}^{R_I} U_{np}\, G_{pcq}\, W_{iq},
\qquad
U \in \mathbb{R}_+^{N \times R_N},\;
\mathcal{G} \in \mathbb{R}_+^{R_N \times C \times R_I},\;
W \in \mathbb{R}_+^{I \times R_I},
$$

i.e., $\mathcal{B} = \mathcal{G} \times_1 U \times_3 W$. No anchored or
fixed components. `model="matrix"` and `model="tensor"` are two
different models; algebraic nesting is not a requirement.

## Relationship to the matrix model

Compared **outside** the models: the simple average
$\bar{B}^{\mathrm{tensor}}_{ci} = \frac{1}{N}\sum_n B_{nci}$ and the
fraction-weighted
$\bar{B}^{\mathrm{tensor,w}}_{ci} = \sum_n A_{nc} B_{nci} / \sum_n A_{nc}$
(matching the matrix model's least squares estimand, whose normal
equations $A^\top A \hat{B} = A^\top X$ weight individuals by their
fractions) versus $\hat{B}^{\mathrm{matrix}}$.

---

# Identifiability (verified 2026-08-23)

## 1. Exact rewrites of the observation model

Cell-type slice form:

$$
B_{\cdot c \cdot} = U G_c W^\top,
\qquad G_c \in \mathbb{R}^{R_N \times R_I},
$$

$$
X = \sum_{c=1}^{C} \mathrm{diag}(A_{\cdot c})\, U G_c W^\top .
$$

Row (per-individual) form, with $a_n$ the $n$-th row of $A$ and
$K_n := \mathcal{G} \times_2 a_n = \sum_c A_{nc} G_c
\in \mathbb{R}^{R_N \times R_I}$:

$$
x_n^\top = u_n^\top K_n W^\top .
$$

Global matrix form via the **row-wise Khatri–Rao (face-splitting)
product** $U \odot_r A \in \mathbb{R}^{N \times R_N C}$ (rows
$u_n \otimes a_n$), with $\mathcal{G}$ unfolded to
$G_{(pc) \times q} \in \mathbb{R}^{R_N C \times R_I}$:

$$
X = (U \odot_r A)\; G_{(pc) \times q}\; W^\top,
\qquad
\mathrm{vec}(X) = \big(W \otimes (U \odot_r A)\big)\, \mathrm{vec}(G).
$$

How $A$ contracts the modes — the single most important structural
fact: **$A$ does not act as an ordinary left factor** (as in the matrix
model $X = AB$). It is contracted *jointly with the unknown $U$ over
the shared individual index $n$* (row-wise Kronecker). Useful rank
bound (rows $u_n \otimes a_n$ lie in
$\mathbb{R}^{R_N} \otimes \mathrm{rowspace}(A)$):

$$
\mathrm{rank}(U \odot_r A) \le \min\big(N,\; R_N\,\mathrm{rank}(A)\big).
$$

Also define $H := (U \odot_r A)\, G_{(pc)\times q}
\in \mathbb{R}^{N \times R_I}$ (rows $h_n^\top = u_n^\top K_n$), so
that $X = H W^\top$ — an exact rank-$R_I$ (NMF-type) factorization.

## 2. Three notions of identifiability

- **(A) Factor identifiability** — $(U, \mathcal{G}, W)$ unique. Never
  holds: the gauge group $\mathrm{GL}(R_N) \times \mathrm{GL}(R_I)$
  ($U \mapsto US$, $W \mapsto WT$,
  $\mathcal{G} \mapsto \mathcal{G} \times_1 S^{-1} \times_3 T^{-1}$)
  preserves $\mathcal{B}$. Harmless for isoTensor's purpose.
- **(B) Latent tensor identifiability** — $\mathcal{B}$ unique given
  $(X, A)$ within the model class. **This is the notion that matters**;
  all claims below are about (B) unless stated.
- **(C) Downstream quantities** — implied by (B) when it holds; some
  survive even when (B) fails (trivially, the observed contraction
  $\sum_c A_{nc} B_{nci} = X_{ni}$; also the matrix-model estimand).
  Group differences / individual programs generally need (B) or at
  least identifiability of the relevant projections — do not assume
  they survive (B)-failure without checking.

## 3. Conditional identifiability (Cases 1–4)

**Case 1 — $U, W$ known, $\mathcal{G}$ unknown.** Linear inverse
problem $\mathrm{vec}(X) = (W \otimes M)\,\mathrm{vec}(G)$ with
$M := U \odot_r A$. $\mathcal{G}$ (equivalently $\mathcal{B}$, since
$U, W$ full column rank make $\mathcal{B} \leftrightarrow \mathcal{G}$)
is identifiable **iff** $\mathrm{rank}(W) = R_I$ and
$\mathrm{rank}(M) = R_N C$. Necessary for the latter:
$\mathrm{rank}(A) = C$, $N \ge R_N C$, and joint genericity of the
pairs $(u_n, a_n)$. [PROVED — linear algebra]

**Case 2 — $W$ known, $(U, \mathcal{G})$ unknown.** $H = X W (W^\top
W)^{-1}$ is recovered exactly; remains the bilinear problem
$h_n^\top = u_n^\top (\mathcal{G} \times_2 a_n)$ with the
$\mathrm{GL}(R_N)$ gauge. Theorem A below (non-identifiability for
$R_N \ge R_I$) lives entirely in this case; for $R_N < R_I$ Theorem C
gives the counting threshold. Any non-identifiability here transfers
a fortiori to Case 4.

**Case 3 — $U$ known, $(\mathcal{G}, W)$ unknown.** If
$M = U \odot_r A$ has full column rank $R_N C$, then
$Y := G_{(pc)\times q} W^\top$ is recovered linearly from $X = M Y$,
and $\mathcal{B}$ is a *linear* function of $Y$:
$B_{nci} = \sum_p U_{np} Y_{(pc), i}$. Hence $\mathcal{B}$ is
identifiable — **regardless of $R_N$ vs $R_I$** — while $(\mathcal{G},
W)$ retain the harmless $\mathrm{GL}(R_I)$ ambiguity. [PROVED —
sufficient direction; the converse is not claimed]
Asymmetry worth noting: knowing $U$ resolves the problem, knowing $W$
does not — the difficulty is concentrated in the entanglement of $U$
with $A$ on the shared index $n$.

**Case 4 — all unknown (actual isoTensor).** Adds to Case 2 the
factorization ambiguity of $X = H W^\top$ (at best $\mathrm{GL}(R_I)$;
under non-negativity, uniqueness of this NMF stage needs
separability / sufficiently-scattered-type conditions). All Case-2
negative results apply. Positive results require simultaneously
resolving the $W$-stage and the Case-2 stage; no sufficient condition
is proved here.

## 4. Main results

### Theorem A ($R_N \ge R_I$ ⇒ $\mathcal{B}$ not identifiable) — PROVED

*Assume $C \ge 2$, rows of $A$ non-zero, true $(U, \mathcal{G}, W)$
with $W$ full column rank. If $R_N \ge R_I$, there is a continuum of
parameters reproducing $X$ exactly with different $\mathcal{B}$.*

Construction: choose **any** alternative core $\mathcal{G}'$ such that
every $K'_n = \mathcal{G}' \times_2 a_n$ has rank $R_I$ (a generic
condition, satisfiable because $R_N \ge R_I$ and $a_n \neq 0$). Then
$\mathrm{rowspace}(K'_n) = \mathbb{R}^{R_I}$, so for each $n$ there is
$u'_n$ with $u'^\top_n K'_n = u_n^\top K_n = h_n^\top$. The triple
$(U', \mathcal{G}', W)$ reproduces $H$, hence $X$, exactly; for generic
$\mathcal{G}'$ the slices $u'^\top_n G'_c W^\top$ differ from
$u^\top_n G_c W^\top$, i.e., $\mathcal{B}' \ne \mathcal{B}$.

Remarks. (i) **No assumption on the variation of $A$** — arbitrary,
even maximal, variation of cell fractions does *not* repair this
regime. (ii) This is an explicit construction, not parameter counting.
(iii) $C \ge 2$ is essential: for $C = 1$, matching the contraction is
matching $\mathcal{B}$ itself.

### Theorem A′ (non-negativity does not rescue interior truth) — PROVED (generic interior truth)

*If in addition $U > 0$ and $\mathcal{G} > 0$ entrywise (strictly) and
each $K_n$ has rank $R_I$, then the construction can be carried out
inside the non-negative cone:* take $\mathcal{G}' = \mathcal{G} +
\varepsilon \mathcal{E}$; then $u'_n \to u_n$ as $\varepsilon \to 0$
(each $K'_n$ invertible on its row space), so small $\varepsilon$
keeps $U' > 0$, $\mathcal{G}' > 0$.

Numerical instance ($R_N = R_I = 2$, $C = 3$, $N = I = 20$, Dirichlet
$A$, $\varepsilon = 1.25 \times 10^{-2}$): all factors strictly
positive, $\|X' - X\|_F / \|X\|_F = 1.4 \times 10^{-16}$,
$\|\mathcal{B}' - \mathcal{B}\|_F / \|\mathcal{B}\|_F = 1.9 \times
10^{-2}$ [`dev/identifiability_check.R`].

Consequences: for **dense (strictly positive) ground truth,
$R_N < R_I$ is necessary** for recoverability, with or without
non-negativity. Non-negativity can only help through *active*
constraints — zeros/sparsity in the true factors (EMPIRICAL
HYPOTHESIS, §8).

### Theorem B (degenerate $A$) — PROVED

*If $\mathrm{rank}(A) = r < C$, then any core perturbation
$\mathcal{E}$ with $\mathcal{E} \times_2 a_n = 0\ \forall n$ — a space
of dimension $\ge (C - r) R_N R_I > 0$ — changes $\mathcal{B}$ but not
$X$.* Interior truth version as in A′. Hence $\mathrm{rank}(A) = C$ is
**necessary**. The extreme $r = 1$ (identical rows) loses the cell
type axis entirely. $\mathrm{rank}(A) = C$ is **not sufficient**
(Theorem A holds for any $A$).

### Theorem C (counting threshold for $R_N < R_I$) — NECESSARY CONDITION

*Conditional on the $W$-stage being resolved (Case 2): unknowns modulo
the $\mathrm{GL}(R_N)$ gauge number $N R_N + R_N C R_I - R_N^2$
against $N R_I$ equations. If*

$$
N (R_I - R_N) \;<\; R_N (C R_I - R_N),
$$

*then at a generic smooth truth the solution fiber has dimension
exceeding the gauge orbit (fiber-dimension / constant-rank argument),
so $\mathcal{B}$ is locally non-identifiable, with predicted local
deficiency $R_N(C R_I - R_N) - N(R_I - R_N)$.* Since Case 2
non-identifiability implies Case 4 non-identifiability, the threshold
$N^* = \lceil R_N (C R_I - R_N) / (R_I - R_N) \rceil$ is a
**necessary**-type sample size condition for the full problem. It is
*not* claimed sufficient.

## 5. Small exact example (pencil and paper): $R_N = R_I = 1$

$B_{nci} = u_n g_c w_i$, so $X = \tilde{u} w^\top$ with
$\tilde{u}_n = u_n (a_n^\top g)$. Pick **any** $g' \in
\mathbb{R}^C_{>0}$ and set $u'_n = u_n (a_n^\top g)/(a_n^\top g')$:
$X$ is unchanged, while

$$
\frac{B'_{nci}}{B_{nci}} = \frac{g'_c\,(a_n^\top g)}{g_c\,(a_n^\top g')}
\neq 1 \quad \text{unless } g' \propto g .
$$

Even the simplest instance of the model is unidentifiable for
$C \ge 2$, with positivity for free — Theorem A at the smallest scale.
(All four $(R_N, R_I) \in \{1,2\}^2$ combinations on $N=3, C=2, I=3$:
$(1,1), (2,1), (2,2)$ fall under Theorem A; $(1,2)$ is the only
candidate-identifiable one, with Theorem C threshold $N^* = 3$.)

## 6. Numerical local identifiability (Jacobian ranks) — LOCAL/GENERIC EVIDENCE

Method (`dev/identifiability_check.R`): exact Jacobians of
$\theta \to \mathcal{B}$ and $\theta \to X$;
$\mathrm{defic} := \mathrm{rank}(J_B) - \mathrm{rank}(J_X) > 0$
certifies a tangent direction changing $\mathcal{B}$ but not $X$
(local non-identifiability); $\mathrm{defic} = 0$ is generic *local*
evidence for identifiability. In every run
$\mathrm{rank}(J_B) = P - R_N^2 - R_I^2$ (gauge dimensions as
expected; $P$ = raw parameter count).

| config | $N$ | $C$ | $I$ | $R_N$ | $R_I$ | $A$ | defic | predicted |
|---|---|---|---|---|---|---|---|---|
| Thm A | 4 | 2 | 4 | 1 | 1 | generic | 1 | 1 |
| Thm A | 6 | 2 | 6 | 2 | 1 | generic | 6 | 6 |
| Thm A | 8 | 2 | 8 | 2 | 2 | generic | 4 | 4 |
| Thm A (old grid!) | 20 | 3 | 20 | 2 | 2 | generic | 8 | 8 |
| Thm A | 10 | 2 | 10 | 3 | 2 | generic | 13 | 13 |
| Thm C, $N < N^*{=}3$ | 2 | 2 | 6 | 1 | 2 | generic | 1 | 1 |
| $N > N^*$ | 4 | 2 | 6 | 1 | 2 | generic | 0 | 0 |
| Thm C, $N < N^*{=}8$ | 6 | 2 | 10 | 2 | 3 | generic | 2 | 2 |
| $N > N^*$ | 10 | 2 | 10 | 2 | 3 | generic | 0 | 0 |
| Thm C, $N < N^*{=}14$ | 12 | 3 | 12 | 2 | 3 | generic | 2 | 2 |
| $N > N^*$ | 16 | 3 | 12 | 2 | 3 | generic | 0 | 0 |
| Thm B ($r=1$) | 16 | 3 | 12 | 2 | 3 | identical rows | 21 | $\ge 12$ |
| Thm B ($r=2$) | 16 | 3 | 12 | 2 | 3 | rank-deficient | 6 | 6 |

Predictions: Theorem A rows use the construction's dimension
$R_N C R_I + N(R_N - R_I) - R_N^2$ — the observed deficiency matches
**exactly in every case**, i.e., the constructed family accounts for
the *entire* local ambiguity. Theorem C rows use
$R_N(C R_I - R_N) - N(R_I - R_N)$ — again exact. Theorem B ($r=2$)
matches $(C-r) R_N R_I = 6$ exactly.

Caveats: local full rank $\ne$ global uniqueness (defic $= 0$ does not
prove global identifiability — that remains EMPIRICAL, probed by
multi-start simulations); Jacobian ranks at interior points ignore the
inequality constraints (fine here: at interior truth they are
inactive, which is exactly Theorem A′'s point).

## 7. Role of $A$

| scenario | status |
|---|---|
| identical rows ($a_n \equiv a$) | unrecoverable for $C \ge 2$ [PROVED, Thm B $r{=}1$; observed defic 21] |
| variation but $\mathrm{rank}(A) = r < C$ | $(C-r) R_N R_I$ unrecoverable directions [PROVED, Thm B; observed exactly] |
| full column rank | **necessary but not sufficient** (Thm A kills $R_N \ge R_I$ at any variation); with $R_N < R_I$, $N \ge N^*$, generic $\theta$: locally identifiable in all tested configs [LOCAL/GENERIC EVIDENCE] |
| strongly correlated columns | identifiability can hold exactly but the problem becomes ill-conditioned — noise amplification; identifiability $\ne$ well-posedness [EMPIRICAL quantification via $\kappa(A)$, column correlations] |
| rare cell type ($A_{\cdot c} \approx 0$) | exact 0: that cell type's slices drop out of $X$ entirely — unidentifiable [PROVED, limiting case]; near 0: error blow-up expected $\propto$ 1/frequency [EMPIRICAL] |

$A$'s columns interact with $U$ through $U \odot_r A$ (§1): what must
be well-conditioned is not $A$ alone but the set
$\{u_n \otimes a_n\}_n$ — a joint condition on cell-fraction variation
*and* individual-loading variation.

## 8. Non-negativity

Analyzed real-valued first (§4–6); what $\ge 0$ can add:

- **At interior truth: nothing, locally.** Constraints are inactive;
  Theorems A′/B interior versions show the ambiguity survives inside
  the cone. [PROVED]
- **Gauge reduction**: with sufficiently sparse factors the absorbable
  transforms shrink from $\mathrm{GL}$ to monomial
  (permutation × scaling) — an NMF-type phenomenon. This affects
  factor identifiability (A-level), not $\mathcal{B}$-level directly.
- **Potential $\mathcal{B}$-level gains only via active (zero)
  constraints**: sparse/boundary ground truth might rescue
  $R_N \ge R_I$ or sub-threshold $N$; plausibly requires
  separability / sufficiently-scattered-type conditions on the
  factors (and on the NMF stage $X = H W^\top$). No proof for this
  model. [EMPIRICAL HYPOTHESIS — sparse-truth settings are a grid
  axis; see `notes/tensor-simulation.md`]

## 9. Classification of all claims

| claim | classification |
|---|---|
| $R_N \ge R_I$ ($C \ge 2$, real factors) ⇒ $\mathcal{B}$ not identifiable | **PROVED** (Thm A, constructive) |
| ... still true under non-negativity for strictly positive truth | **PROVED** (Thm A′, generic interior truth; numerical instance) |
| $R_N \ge R_I$ with sparse/boundary truth may be recoverable | EMPIRICAL HYPOTHESIS |
| $\mathrm{rank}(A) = C$ necessary | **PROVED** (Thm B) |
| some variation of $A$ across individuals necessary | **PROVED** (Thm B, $r=1$ case) |
| full-rank / high-variation $A$ sufficient | **FALSE** (disproved by Thm A) |
| $N(R_I - R_N) \ge R_N(C R_I - R_N)$ (for $R_N < R_I$) | NECESSARY CONDITION (Case-2-conditional proof sketch; numerically exact) |
| $R_N < R_I$ + threshold + generic $A, \theta$ ⇒ locally identifiable | LOCAL/GENERIC EVIDENCE (all tested configs, defic $=0$) |
| ... ⇒ globally identifiable | EMPIRICAL HYPOTHESIS (multi-start simulation) |
| Case 1 / Case 3 conditions | **PROVED** |
| noise robustness, $A$-error propagation, rare-cell-type error scaling | EMPIRICAL |

**Verdict on the original claim** ("$R_N \ge R_I$ breaks
sign-unconstrained identifiability, generic"): correct, and now
**PROVED** — the previous derivation took an unnecessary two-stage
detour ("given the factorization $X = HW^\top$…") and was flagged as
generic evidence; the direct construction needs no $W$-stage
assumption. The claim is also *strengthened*: the earlier note left
"whether non-negativity rescues this regime" fully open — for interior
truth it provably does not; only sparse/boundary truth remains open.

---

## Comparison of candidate structures (kept for the record)

| | (a) CP | (b) full Tucker | (c) plain partial Tucker | (d) shared + deviation |
|---|---|---|---|---|
| Identifiability | best (Kruskal) | worst | see §4–9 above | needs centering ⇒ conflicts with $\ge 0$ |
| Cell type axis kept explicit | no (mixed in $V$) | no ($V$ compresses) | **yes** | yes |
| Parameters | $R(N{+}C{+}I)$ | $NR_N{+}CR_C{+}IR_I{+}R_NR_CR_I$ | $NR_N{+}R_NCR_I{+}IR_I$ | $CI$ + dev. |
| MU derivable | yes | yes | **yes** (to be derived after the simulation study) | no (signed / non-polynomial) |

- **(a) CP**: strongest uniqueness theory, but a single component mixes
  cell types through its loading vector — against the design
  constraint.
- **(b) full Tucker**: compresses the cell type mode; rejected.
- **(c) plain partial Tucker**: current candidate.
- **(d) shared + deviation**: centering incompatible with $D \ge 0$;
  rejected as a model (population comparison done externally instead).

## Rejected designs

**Anchored partial Tucker (rejected 2026-08-23).** Fixing
$U_{\cdot 1} = \mathbf{1}_N$, warm-starting $G_{1\cdot\cdot}$ from the
Phase 1 $\hat{B}$, and reading components $p \ge 2$ as non-negative
individual deviations. Not adopted: (i) upward-only deviations are an
unacceptable modeling bias; (ii) exact nesting of the matrix model is
not a design requirement. (The *generator*
`isoToyModel(model="tensor", effect_cell_types=...)` uses a fixed
shared component internally to construct ground truth — a simulation
construction, not a fitting-model constraint.)

## API note (vs nnTensor)

`nnTensor::NTD` uses `rank=c(R1,R2,R3)` plus `modes=` to choose
compressed modes. For isoTensor the cell type mode is never
compressed, so `rank_individual` / `rank_isoform` (already reserved in
`isoTensor()` and used by `isoToyModel(model="tensor")`) are the
recommended API. Everything else follows nnTensor conventions.

## Open questions for review

1. Do the theory-driven grid changes in `notes/tensor-simulation.md`
   (all three rank regimes; dense $R_N \ge R_I$ as negative control;
   sparsity axis) cover what should be checked before MU derivation?
2. Objective for Phase 2: Frobenius vs KL/beta divergence.
3. Per-gene simplex constraints on $W$: Phase 2 or later.
