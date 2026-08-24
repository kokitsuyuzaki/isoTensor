# isoTensor — Project Status

作成日: 2026-08-24 (project close-out 時点、commit bc5e451 以降)。
科学的結論の横断整理は `isoTensor-experiments/reports/FINAL_SUMMARY.md`。

## Implemented

- `model="matrix"`: population-level X ~ AB、multiplicative updates
  (per-isoform NNLS と同値な凸問題)
- `model="tensor"`: partial non-negative Tucker
  (B_nci = Σ U_np G_pcq W_iq、cell type mode 非圧縮、A 固定)
- optimizer = `"lbfgs"` (default; log 再パラメータ化 L-BFGS-B, 解析勾配,
  `optimizer.control=list(maxit, factr)`) / `"mu"` (単調 MU baseline) /
  `"hybrid"` (MU → L-BFGS refine)
- `isoToyModel()`: matrix / tensor の合成データ生成 (Dirichlet A,
  A_mean / A_variation dial, individual effect の制御, noise / noise_A)
- 診断・警告: R_N ≥ R_I の generic non-identifiability warning、
  Theorem C の sample-size warning、非負性・次元チェック
- tests: testthat (matrix / tensor / lbfgs / toy generator)、
  vignette、hand-written man、NEWS

## Validated (isoTensor-experiments による)

- matrix mode と Sciege-style per-isoform NNLS の同値性
  (synthetic で機械精度、GTEx-scale 実データで cor 0.999999)
- tensor mode の synthetic recovery (noise-free candidate regime で
  relB 1e-3〜1e-15; full grid 620 datasets)
- 理論的 non-identifiability regime (R_N ≥ R_I) の負対照挙動
  (relX 小 / relB 不定 / objective–relB 無相関)
- optimizer 比較 (L-BFGS ≫ MU; hybrid は L-BFGS 単独と同等;
  thr=1e-8 が MU の実用停止点; 検証済み default 値 maxit=5000, factr=10)

## Experimental / limited

- individual tensor の real-data 適用: GTEx Heart LV / Muscle の診断では
  **Red / Yellow 下位** (cond(A) 高、measurement noise、rare cell type、
  解の不安定性)。raw TPM スケールでは L-BFGS の line-search 失敗が増える
  (入力スケーリング未実装)。population-level (matrix) は実データで Green。

## Not implemented (将来課題)

- automatic rank selection (R_N, R_I の選択; Gabriel-style holdout 等)
- multi-start API (`n.start` など; 現在は呼び出し側の責務)
- KL / beta-divergence 目的関数
- spillover-aware individual-level inference
- A の不確実性の伝播 (A 固定の緩和)
- gene–isoform 階層 (isoform group / simplex 制約)
- 実データ向け入力スケーリング (L-BFGS の数値安定化)
