# Phase 1: Matrix model — derivation of the multiplicative update

## Model

- $X \in \mathbb{R}_+^{N \times I}$: bulk isoform abundance
  (individual $\times$ isoform), observed.
- $A \in \mathbb{R}_+^{N \times C}$: cell type fractions
  (individual $\times$ cell type), externally estimated and **fixed**.
- $B \in \mathbb{R}_+^{C \times I}$: population-level cell-type-resolved
  isoform abundance, to be estimated.

$$
X \approx A B,
\qquad
\min_{B \ge 0} \; f(B) = \frac{1}{2}\|X - AB\|_F^2 .
$$

Because $A$ is fixed, this is a non-negative least squares (NNLS) problem
in $B$; each column of $B$ (one isoform) is an independent NNLS problem
with the shared design matrix $A$. $f$ is convex in $B$, and strictly
convex iff $A$ has full column rank ($\mathrm{rank}(A) = C$, requiring
$N \ge C$).

## Gradient

$$
\nabla_B f
= A^\top (AB - X)
= \underbrace{A^\top A B}_{\text{non-negative}}
- \underbrace{A^\top X}_{\text{non-negative}} .
$$

Element-wise:

$$
\frac{\partial f}{\partial B_{ci}}
= [A^\top A B]_{ci} - [A^\top X]_{ci} .
$$

Both terms are non-negative whenever $A, B, X \ge 0$, so the gradient
splits naturally as $\nabla_B f = \nabla^+ - \nabla^-$ with
$\nabla^+ = A^\top A B$ and $\nabla^- = A^\top X$.

## Multiplicative update

The standard MU heuristic (equivalently, gradient descent with the
adaptive per-element step size $\eta_{ci} = B_{ci} / [A^\top A B]_{ci}$)
gives

$$
B_{ci} \leftarrow B_{ci} \cdot
\frac{[\nabla^-]_{ci}}{[\nabla^+]_{ci}}
= B_{ci} \cdot
\frac{[A^\top X]_{ci}}{[A^\top A B]_{ci} + \varepsilon},
$$

i.e., in matrix notation

$$
B \leftarrow B \odot \frac{A^\top X}{A^\top A B + \varepsilon},
$$

where $\varepsilon$ (machine epsilon, argument `pseudocount`) guards
against division by zero. This coincides with the update proposed in the
task description, and is exactly the $V$-update of Lee–Seung Frobenius
NMF with the other factor fixed (cf. `nnTensor::NMF` with
`algorithm="Frobenius", fixU=TRUE`, where the roles are
$U = A$, $V = B^\top$).

## Properties

- **Monotonicity.** By the Lee–Seung auxiliary (majorize–minimize)
  argument, each update does not increase $f$. Since each column of $B$
  is an independent NNLS problem, monotonicity holds column-wise and
  hence globally. Verified numerically in the unit tests
  (`RecError` history is non-increasing).
- **Fixed points.** Stationary points of the update satisfy the KKT
  conditions of the NNLS problem: $B_{ci} > 0 \Rightarrow
  [\nabla_B f]_{ci} = 0$; entries at $0$ stay at $0$ (absorbing), which
  is why the default initialization is strictly positive
  (`runif`, followed by an optimal global rescaling of $B$ to speed up
  early iterations).
- **Recovery.** If $X = A B_{\mathrm{true}}$ exactly and $A$ has full
  column rank, the unique global minimizer is $B_{\mathrm{true}}$ (the
  unconstrained least squares solution, which is feasible), and MU from
  a positive initialization converges to it. Verified on the toy data
  ($N=5$, $C=2$, $I=6$): $\max|B - B_{\mathrm{true}}| \approx 5 \times
  10^{-15}$ after ~400 iterations.
- **Non-identifiability caveat.** If $\mathrm{rank}(A) < C$ (e.g.,
  collinear cell type fractions across individuals), the minimizer is
  not unique and the recovered $B$ depends on the initialization. This
  is a property of the data, not the algorithm; a warning-level
  diagnostic may be added later.

## Implementation notes

- $A^\top X$ and $A^\top A$ are constant because $A$ is fixed; they are
  computed once before the loop (`.initIsoTensor`), so each iteration
  costs only $O(C^2 I)$ — independent of $N$.
- Convergence criterion: relative change of the reconstruction error
  $\mathrm{RecError} = \|X - AB\|_F$ falls below `thr`, or `num.iter`
  iterations are reached. The objective history is available as
  `fit$RecError` ($f = \mathrm{RecError}^2 / 2$).
- Zeros in $X$ are allowed as-is (the Frobenius objective does not
  require positivity of $X$, unlike KL/IS divergences).
