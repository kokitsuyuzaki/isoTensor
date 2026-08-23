# Phase 2: Multiplicative updates for the plain partial non-negative Tucker model

Derivation note (written BEFORE the implementation, 2026-08-23).
Model and identifiability: `notes/tensor-model.md`. Simulation design:
`notes/tensor-simulation.md`.

## 1. Objective and reconstruction

Data $X \in \mathbb{R}_+^{N \times I}$, fixed
$A \in \mathbb{R}_+^{N \times C}$ (never updated). Factors

$$
U \in \mathbb{R}_+^{N \times R_N}, \qquad
\mathcal{G} \in \mathbb{R}_+^{R_N \times C \times R_I}, \qquad
W \in \mathbb{R}_+^{I \times R_I},
$$

latent tensor and reconstruction:

$$
B_{nci} = \sum_{p,q} U_{np} G_{pcq} W_{iq},
\qquad
\hat{X}_{ni} = \sum_{c} A_{nc} B_{nci}
= \sum_{c,p,q} A_{nc} U_{np} G_{pcq} W_{iq}.
$$

Objective:

$$
L(U, \mathcal{G}, W)
= \frac{1}{2} \| X - \hat{X} \|_F^2
= \frac{1}{2} \sum_{n,i} \Big( X_{ni} - \sum_{c,p,q} A_{nc} U_{np} G_{pcq} W_{iq} \Big)^2,
\qquad U, \mathcal{G}, W \ge 0 .
$$

Useful intermediates (all re-used by the updates):

| symbol | definition | size |
|---|---|---|
| $G_c$ | $c$-th cell type slice of $\mathcal{G}$ | $R_N \times R_I$ |
| $K_n$ | $\mathcal{G} \times_2 a_n = \sum_c A_{nc} G_c$ | $R_N \times R_I$ |
| $H$ | $H_{nq} = \sum_{p,c} A_{nc} U_{np} G_{pcq}$, i.e., rows $h_n^\top = u_n^\top K_n$; equivalently $H = \sum_c \mathrm{diag}(A_{\cdot c})\, U G_c$ | $N \times R_I$ |
| $\hat{X}$ | $H W^\top$ | $N \times I$ |
| $R$ | $\hat{X} - X$ (residual) | $N \times I$ |

## 2. Gradients (index derivation, then matrix form)

Write $R_{ni} = \hat{X}_{ni} - X_{ni}$, so
$\partial L / \partial \theta_j = \sum_{n,i} R_{ni}\,
\partial \hat{X}_{ni} / \partial \theta_j$.

### 2.1 Gradient w.r.t. $U$

$$
\frac{\partial \hat{X}_{ni}}{\partial U_{mp}}
= \delta_{nm} \sum_{c,q} A_{nc} G_{pcq} W_{iq}
= \delta_{nm}\, (K_n W^\top)_{p i},
$$

$$
\frac{\partial L}{\partial U_{mp}}
= \sum_{i} R_{mi} \sum_{c,q} A_{mc} G_{pcq} W_{iq}
= \sum_{c} A_{mc} \sum_{q} \Big(\sum_i R_{mi} W_{iq}\Big) G_{pcq}
= \sum_{c} A_{mc} \big[ (R W)\, G_c^\top \big]_{mp}.
$$

Matrix form ($N \times R_N$):

$$
\nabla_U L = \sum_{c} \mathrm{diag}(A_{\cdot c})\, (R W)\, G_c^\top .
$$

Since $R = \hat{X} - X$ and $A, W, G_c \ge 0$, the split into
non-negative parts is

$$
\nabla_U^+ = \sum_{c} \mathrm{diag}(A_{\cdot c})\, (\hat{X} W)\, G_c^\top,
\qquad
\nabla_U^- = \sum_{c} \mathrm{diag}(A_{\cdot c})\, (X W)\, G_c^\top .
$$

### 2.2 Gradient w.r.t. $\mathcal{G}$

$$
\frac{\partial \hat{X}_{ni}}{\partial G_{pcq}} = A_{nc} U_{np} W_{iq},
\qquad
\frac{\partial L}{\partial G_{pcq}}
= \sum_{n,i} R_{ni} A_{nc} U_{np} W_{iq}
= \sum_{n} U_{np} A_{nc} (R W)_{nq}
= \big[ U^\top \mathrm{diag}(A_{\cdot c})\, (R W) \big]_{pq}.
$$

Matrix form, per cell type slice ($R_N \times R_I$ each):

$$
\nabla_{G_c} L = U^\top \mathrm{diag}(A_{\cdot c})\, (R W),
\qquad
\nabla_{G_c}^{\pm} = U^\top \mathrm{diag}(A_{\cdot c})\, (\hat{X} W \ \text{or}\ X W).
$$

### 2.3 Gradient w.r.t. $W$

$$
\frac{\partial \hat{X}_{ni}}{\partial W_{jq}}
= \delta_{ij} \sum_{p,c} A_{nc} U_{np} G_{pcq}
= \delta_{ij} H_{nq},
\qquad
\frac{\partial L}{\partial W_{jq}}
= \sum_{n} R_{nj} H_{nq}
= (R^\top H)_{jq}.
$$

Matrix form ($I \times R_I$):

$$
\nabla_W L = R^\top H,
\qquad
\nabla_W^+ = \hat{X}^\top H = W (H^\top H),
\qquad
\nabla_W^- = X^\top H .
$$

(The identity $\hat{X}^\top H = W H^\top H$ follows from
$\hat{X} = H W^\top$ and gives the cheaper form.)

## 3. Multiplicative updates

With the generic MU scheme
$\theta \leftarrow \theta \odot \nabla^-_\theta L / (\nabla^+_\theta L + \varepsilon)$:

$$
U \leftarrow U \odot
\frac{\sum_c \mathrm{diag}(A_{\cdot c})\,(X W)\, G_c^\top}
     {\sum_c \mathrm{diag}(A_{\cdot c})\,(\hat{X} W)\, G_c^\top + \varepsilon}
$$

$$
G_c \leftarrow G_c \odot
\frac{U^\top \mathrm{diag}(A_{\cdot c})\,(X W)}
     {U^\top \mathrm{diag}(A_{\cdot c})\,(\hat{X} W) + \varepsilon}
\qquad (c = 1, \dots, C)
$$

$$
W \leftarrow W \odot
\frac{X^\top H}
     {W (H^\top H) + \varepsilon}
$$

$A$ never appears on the left-hand side: it is fixed.
$\hat{X}$ (and $H$) are recomputed from the current factors before each
block update, and the blocks are updated in the order
$U \to \mathcal{G} \to W$ (each block update is individually monotone,
§5, so the order affects the path but not the monotonicity).

### Dimensions of every term

| term | size |
|---|---|
| $XW$, $\hat{X}W$ | $N \times R_I$ |
| $(XW) G_c^\top$ | $N \times R_N$ |
| $\mathrm{diag}(A_{\cdot c}) \cdot$ (row scaling by $A_{\cdot c}$) | $N \times \cdot$ |
| numerator/denominator for $U$ | $N \times R_N$ |
| $U^\top \mathrm{diag}(A_{\cdot c}) (XW)$ | $R_N \times R_I$ |
| numerator/denominator for $G_c$ | $R_N \times R_I$ |
| $H$, $X^\top H$ | $N \times R_I$, $I \times R_I$ |
| $H^\top H$ | $R_I \times R_I$ |
| numerator/denominator for $W$ | $I \times R_I$ |

### Efficient computation (as implemented)

- `diag(A[,c]) %*% M` is computed as `A[,c] * M` (R column-major
  recycling = row scaling); no diagonal matrix is formed.
- $XW$ is computed once per $U$/$\mathcal{G}$ update; $H$ once per
  $W$ update / reconstruction; $\hat{X} = H W^\top$.
- Per iteration cost: $O(N I R_I)$ for the $X W$/$\hat X W$/$X^\top H$
  products, plus $O(C N R_N R_I)$ for the $c$-loops — linear in every
  dimension, no $N \times C \times I$ tensor is ever materialized
  during fitting (the returned `B` is formed once at the end).

## 4. Monotonicity (MM) — what is guaranteed and why

$\hat{X}$ is **linear in each block** with the other blocks fixed, with
non-negative design matrices:

- $\mathcal{G}$-block: $\mathrm{vec}(\hat{X}) = \big(W \otimes (U
  \odot_r A)\big) \mathrm{vec}(\mathcal{G})$ — a single non-negative
  least squares (NNLS) problem in $\mathcal{G}$; the update in §3 is
  exactly the Lee–Seung Frobenius MU
  $\theta \odot (D^\top y)/(D^\top D \theta)$ for this design
  ($D^\top \mathrm{vec}(X)$ and $D^\top \mathrm{vec}(\hat X)$ are
  precisely the numerator/denominator in §2.2).
- $W$-block: $X^\top \approx W H^\top$ — the classical one-sided NMF
  update with the other factor ($H$) fixed.
- $U$-block: row-wise, $x_n \approx (W K_n^\top)\, u_n$ — $N$
  independent NNLS problems with designs $W K_n^\top \ge 0$; §3's $U$
  update is the Lee–Seung MU applied to each row.

Therefore each block update is an exact instance of the Lee–Seung
multiplicative update for an NNLS subproblem, and the classical
auxiliary-function (majorize–minimize) argument applies verbatim per
block: **$L$ is non-increasing under each of the three updates, hence
under the full sweep — this is theoretically guaranteed (for
$\varepsilon = 0$), not merely empirical.** The structured contraction
with the fixed $A$ does not break the guarantee, because $A$ enters
only through the fixed non-negative design matrices
($W \otimes (U \odot_r A)$, $H$, $W K_n^\top$).

Caveat on $\varepsilon$: the implementation adds
$\varepsilon = $ `pseudocount` (machine epsilon by default) to the
denominators to avoid division by zero. The classical proof is for
$\varepsilon = 0$; the $\varepsilon$-guard shortens each step by a
relative $O(\varepsilon)$ amount and its effect is below numerical
noise. Monotonicity with the guard is verified by unit tests
(non-increase up to a $10^{-10}$ tolerance) — we do not claim a
separate theorem for $\varepsilon > 0$.

Fixed points: stationary points of the sweep satisfy the KKT
conditions of each NNLS subproblem. As with all MU, entries that reach
exactly 0 stay 0 (zero-locking):

- default initialization is **strictly positive**
  ($\mathrm{Unif}(0.5, 1.5)$ entries, then a single global rescale
  absorbed into $\mathcal{G}$: $\alpha = \langle X, \hat{X}\rangle /
  \|\hat{X}\|_F^2$), so no entry starts locked;
- exact zeros arising later from zero numerators are legitimate KKT
  boundary points, not bugs;
- user-supplied `initU`/`initG`/`initW` containing zeros will keep
  those zeros locked — documented, not "fixed".

## 5. Normalization (numerical stabilization only)

Partial Tucker has a scaling ambiguity: column $p$ of $U$ can be
multiplied by $s_p > 0$ if slice $G_{p \cdot \cdot}$ is divided by
$s_p$ (same for $W$ columns vs $G_{\cdot \cdot q}$), leaving both
$\mathcal{B}$ and $\hat{X}$ unchanged. To keep the blocks on
comparable scales, after every sweep the columns of $U$ and $W$ are
rescaled to unit $\ell_2$ norm and the scales are absorbed into
$\mathcal{G}$:

$$
U_{\cdot p} \leftarrow U_{\cdot p}/s_p,\quad
G_{p \cdot \cdot} \leftarrow s_p\, G_{p\cdot\cdot};\qquad
W_{\cdot q} \leftarrow W_{\cdot q}/t_q,\quad
G_{\cdot\cdot q} \leftarrow t_q\, G_{\cdot\cdot q}.
$$

A unit test verifies that this leaves $\mathcal{B}$ and $\hat{X}$
unchanged to machine precision. This is **numerical stabilization
only** — it does not (and cannot) improve identifiability of
$\mathcal{B}$; see `notes/tensor-model.md`.

## 6. Initialization

Plain strictly-positive random initialization (reproducible via the
caller's `set.seed`), as in §4. The anchored initialization from the
Phase 1 population-level $\hat{B}$ is **not used**; recorded here only
as a possible future `init` option (it was rejected as a *model
constraint*; as a mere starting point it would be legitimate but is
out of scope for the minimal implementation).

## 7. Verification protocol (done before/while implementing)

1. Analytical gradients of §2 checked element-wise against central
   finite differences on $N{=}3, C{=}2, I{=}4, R_N{=}1, R_I{=}2$
   (and a second config), tolerance: max relative error $< 10^{-6}$
   with step $h = 10^{-6}$. MU is implemented only after this passes.
2. One-step objective non-increase for each block update separately.
3. Full-history non-increase of `RecError`.
4. Positive control ($R_N < R_I$, identifiable regime) and negative
   control ($R_N \ge R_I$, proven non-identifiable) on
   `isoToyModel(model="tensor")`, plus multi-start —
   see `dev/tensor_mu_check.R` and the report.
