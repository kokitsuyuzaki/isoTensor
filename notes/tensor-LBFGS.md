# Log-reparameterized L-BFGS for the partial non-negative Tucker model

Written BEFORE the implementation (same policy as `tensor-MU.md`).
This is a **port of the implementation validated in the separate
`isoTensor-experiments` repository** (`src/lbfgs_utils.R`,
`src/optimizer_alt_fit.R`; reports `optimizer-calibration.md`,
`2026-08-24-full-recoverability.md`, `2026-08-24-realistic-noise-stress.md`),
not a new algorithm. Where this note and the experiments code could
diverge, the experiments code wins; §8 records the numerical
equivalence check.

## 1. Objective (unchanged)

Same model and objective as `tensor-MU.md`:

    B_nci = sum_{p,q} U_np G_pcq W_iq,
    Xhat_ni = sum_c A_nc B_nci,
    L(U, G, W) = 1/2 ||X - Xhat||_F^2,

with U (N x R_N), G (R_N x C x R_I), W (I x R_I) all non-negative and
A (N x C) fixed. With G_c = G[, c, ] and
H = sum_c diag(A[,c]) U G_c (N x R_I), Xhat = H W^T.

`algorithm` (the objective / divergence; currently only "Frobenius")
and `optimizer` (the method used to minimize it; "mu", "lbfgs",
"hybrid") are conceptually orthogonal. This note is about the
optimizer only; the objective is identical for all three.

## 2. Log reparameterization

The non-negativity constraints are removed by the elementwise
exponential map

    U = exp(Theta_U),   G = exp(Theta_G),   W = exp(Theta_W),

with Theta_U in R^{N x R_N}, Theta_G in R^{R_N x C x R_I},
Theta_W in R^{I x R_I} unconstrained. The image is the strictly
positive orthant, so every iterate is feasible by construction (the
closure U >= 0 is reached only in the limit Theta -> -inf; in practice
entries can become numerically 0 after exp underflow, which is
harmless for the objective).

The unconstrained problem min_Theta L(exp(Theta_U), exp(Theta_G),
exp(Theta_W)) is solved with `stats::optim(method = "L-BFGS-B")`
(no box constraints are supplied; L-BFGS-B without bounds is plain
L-BFGS). This is exactly what the experiments implementation does.

## 3. Gradient (chain rule)

The analytical gradients of L in the original parameters are already
derived and finite-difference-verified in `tensor-MU.md` §2 and
implemented in `.gradIsoTensor`:

    grad_U L  = sum_c diag(A[,c]) (R W) G_c^T        (N x R_N)
    grad_Gc L = U^T diag(A[,c]) (R W)                (R_N x R_I per c)
    grad_W L  = R^T H                                (I x R_I)

with R = Xhat - X. Since dU_np/dTheta_U,np = exp(Theta_U,np) = U_np
(elementwise, diagonal Jacobian), the chain rule gives

    dL/dTheta_U = U  (Hadamard) grad_U L
    dL/dTheta_G = G  (Hadamard) grad_G L
    dL/dTheta_W = W  (Hadamard) grad_W L.

No new gradient derivation is needed; the implementation multiplies
`.gradIsoTensor` output elementwise by the current factors. (The
finite-difference check of `.gradIsoTensor` in the unit tests, plus a
direct FD check of the packed gradient in theta, covers this.)

## 4. Parameter packing

`optim()` takes a single numeric vector. Packing convention (identical
to `isoTensor-experiments/src/lbfgs_utils.R`):

    theta = ( vec(Theta_U),  vec(Theta_G),  vec(Theta_W) )
    lengths:  N*R_N          R_N*C*R_I      I*R_I

where vec() is R's native column-major vectorization: `as.vector()`
of the matrix / array. For the 3-d array G the linear index runs p
fastest, then c, then q, i.e. `as.vector(G)` of `dim = c(R_N, C, R_I)`.
Unpacking reverses this with `matrix(exp(.), N, R_N)`,
`array(exp(.), dim = c(R_N, C, R_I))`, `matrix(exp(.), I, R_I)`, so
pack -> unpack is the identity on strictly positive factors (up to
exp(log(x)) rounding, ~1 ulp); a unit test asserts this roundtrip and
the dimension consistency of the reconstructed factors.

Packing applies `log(pmax(x, eps))` with eps = 1e-10 (the experiments
value): zeros in a starting point are raised to 1e-10 before the log.
Consequence for user-supplied initial factors:

- optimizer = "mu": zeros stay locked at zero (multiplicative updates).
- optimizer = "lbfgs" / "hybrid" (refinement stage): zeros are replaced
  by 1e-10 and CAN move away from zero (no zero-locking).

Both behaviors are documented in the man page.

## 5. Numerical stability of exp / log

Checked against the experiments implementation, which uses **no
clipping, no safe-exp, no rescaling inside the optimization**:

- Overflow: exp(Theta) overflows only for Theta > ~709. Starting from
  log(Unif(0.5, 1.5)) (|Theta| < 0.7) and minimizing a coercive-in-
  scale objective, iterates never approach this; overflow would make
  L = Inf and L-BFGS-B's line search rejects such steps.
- Underflow: exp(Theta) -> 0 for Theta < ~-745 flushes to 0 gracefully
  (the gradient dL/dTheta = x * grad_x also -> 0; the parameter just
  stops moving, which is the correct limit at the boundary).
- Robustness: `optim()` can still fail in rare degenerate line-search
  states (experiments: 4/450 = 0.9% of random-init runs on the full
  grid, 1/270 = 0.4% under noise). The port wraps `optim()` in
  `tryCatch`: optimizer="lbfgs" stops with an informative error
  (suggesting a different seed / "hybrid" / "mu"); optimizer="hybrid"
  falls back to the MU-stage solution with a warning (the MU solution
  is always feasible).

No processing that would change the objective is added.

## 6. Scaling ambiguity / normalization

The partial Tucker scaling gauge (U, W column scales absorbable into
G) is NOT touched during the optimization: renormalizing between
`optim()` iterations would change the parameterization mid-run and
break gradient consistency, and `optim()` offers no safe hook anyway.
The experiments implementation likewise never normalizes inside
L-BFGS and was stable on 620 datasets (0 refine failures). Policy:

- initialization: factors are used as-is (for "hybrid" the MU stage
  already returns normalized U, W by its own per-sweep normalization);
- after `optim()` returns: one final `.normalizeFactorsTensor` call
  absorbs the column scales of U and W into G (B and Xhat provably
  unchanged; unit-tested). This keeps the returned factors on the same
  gauge convention as the MU path. No identifiability improvement is
  claimed.

Note: the "lbfgs" random-init path deliberately does NOT apply the MU
path's initial global rescaling of G, to remain numerically identical
to the validated experiments implementation (§8).

## 7. Initialization and the three optimizers

- optimizer = "lbfgs": strictly positive random init from the same
  distribution and draw order as the MU path (U, G, W ~ Unif(0.5, 1.5),
  drawn U -> G -> W), or user-supplied initU/initG/initW (zeros raised
  to 1e-10, §4). Then a single L-BFGS run.
- optimizer = "hybrid": the existing MU loop (controlled by num.iter /
  thr / pseudocount, unchanged) runs first from the same init rules as
  optimizer="mu"; its solution (zeros raised to 1e-10) is the L-BFGS
  starting point. One init -> MU -> L-BFGS; multi-start stays in
  isoTensor-experiments.
- optimizer = "mu": the existing MU loop, unchanged (monotone
  Lee-Seung baseline, `tensor-MU.md`).

L-BFGS stage controls: `optimizer.control = list(maxit = 5000,
factr = 10)` (defaults = the experiments calibration values;
`factr` is optim's relative tolerance in units of machine epsilon,
so factr = 10 is a tight ~2e-15 relative tolerance). `num.iter` and
`thr` keep their MU meaning and are ignored by the pure "lbfgs" path.

## 8. Why this is a port, and the equivalence check

Evidence from isoTensor-experiments (all on datasets in the candidate
recoverable regime R_N < R_I):

- calibration (3 datasets x 5 inits): L-BFGS reaches relB 5-100x
  better than MU 20k in ~1/10 the time (easy_high 0.037 -> 1.5e-4,
  mod_a20 0.206 -> 0.013, hard_low 0.242 -> 0.044);
- full recoverability grid (620 datasets): hybrid two-stage with this
  L-BFGS refinement, 0 failures in 3,100 refine jobs; L-BFGS alone
  is equal-or-better than hybrid (15 datasets x 30 inits);
- realistic-noise stress: selected relB of L-BFGS alone and hybrid
  identical on all 9 datasets, L-BFGS alone slightly faster.

The port is verified numerically in `dev/lbfgs_port_check.R`: the
same calibration datasets (easy_high, mod_a20, hard_low; init seeds
1-5) are re-fit with `isoTensor(optimizer="lbfgs")` and the resulting
objective / relX / relB are compared against the archived experiments
results (`results/optimizer/alt_summary.tsv`); with identical RNG
draws, identical arithmetic and identical optim settings the
trajectories must coincide.
