# Phase 2 設計改訂報告 (2026-08-23)

指示どおりanchoring案を廃止し、plain partial Tuckerを第一候補として
`notes/tensor-model.md` を改訂、`notes/tensor-simulation.md` を新規作成、
ground-truth generatorを `isoToyModel(model="tensor")` として実装した。
fitting algorithm / MUは実装していない。

変更ファイル: `notes/tensor-model.md`(改訂)、`notes/tensor-simulation.md`(新規)、
`PLAN.md`、`R/isoToyModel.R`、`man/isoToyModel.Rd`、`NAMESPACE`、
`tests/testthat/test_isoToyModel.R`、`tests/testthat/test_isoTensor.R`。
`R CMD check --no-manual`: 0 errors / 0 warnings / 0 notes。

## 1. plain partial Tuckerの最終的な数式

    B_nci = Σ_{p=1}^{R_N} Σ_{q=1}^{R_I} U_np G_pcq W_iq
    U ∈ R_+^{N×R_N},  G ∈ R_+^{R_N×C×R_I},  W ∈ R_+^{I×R_I}

すなわち B = G ×₁ U ×₃ W(cell type mode = mode 2 は非圧縮)。
観測モデルは X_ni ≈ Σ_c A_nc B_nci、A は固定。
anchoringなし(U の第1列固定・population component・Phase 1 warm start は
すべて削除し、tensor-model.md の "Rejected designs" に理由付きで記録)。

## 2. Unknown parameter count

生パラメータ数: P = N·R_N + R_N·C·R_I + I·R_I(観測数 N·I)。
ただし P < NI は必要条件的なcountに過ぎず十分条件ではない。
factor不定性(S∈GL(R_N), T∈GL(R_I) が G に吸収され B 不変)は
B-recoveryには無害。致命的なのは「同じ X を与える異なる B」で、
これがsimulationの調査対象(詳細は tensor-model.md の
"Parameter count and identifiability")。

## 3. Matrix modelとの役割の違い

- matrix: population-level B ∈ R_+^{C×I} を推定。推定量は正規方程式
  AᵀA B̂ = AᵀX による A-weighted な population profile。
- tensor: individual-specific B ∈ R_+^{N×C×I} を推定。
- nestingは要求しない。比較はモデルの外で行う:
  単純平均 B̄_ci = (1/N)Σ_n B_nci と、matrix modelの推定対象に対応する
  weighted版 B̄^w_ci = Σ_n A_nc B_nci / Σ_n A_nc の両方を用意。

## 4. Simulation generatorの設計

新規APIは作らず `isoToyModel(model="tensor", ...)` に統合
(main APIの `isoTensor(model=)` と対応)。生成手順:
U, G, W(runif)→ B = G ×₁ U ×₃ W → A(Dirichlet行)→
X_ni = Σ_c A_nc B_nci(noise-free)→ 任意で X にGaussian noise、
A に摂動(A_obs として返却、X は常に真の A から生成)。

切替: `individual_effect=FALSE` で B_nci = B_ci(Simulation 1)、
`effect_cell_types=` で個人効果を特定cell typeに限定(Simulation 2、
generator内部でのみ共有component構造を使用 — fitting modelの制約ではない)。
Simulation 3(phenotype 2群)は将来 `group` 引数として設計のみ記載。

検証済み(unit tests): 次元・非負性・行和1・X=収縮の厳密一致・
B=(U,G,W)再構成の厳密一致・seed再現性・individual_effect=FALSEで
Bがn不変・effect_cell_types=1でcell type 1のみ変動・
matrix modelがSimulation 1データで共有sliceを<1e-6で復元。

## 5. A variationの生成と定量化

行を Dirichlet(α₀·m) から生成。Var(A_nc) = m_c(1−m_c)/(α₀+1) なので
α₀ 一つでvariationを制御:
none(∞, 同一行)/ low(200)/ moderate(20)/ high(2)。
数値指定で連続sweepも可能。rare cell typeは `A_mean`(例 c(0.9,0.05,0.05))。
毎回 `A_summary` に記録: rank(A)(SVD)、condition number、
列間相関行列、平均fraction・最小平均fraction、total variance Σ_c Var_n(A_nc)。
動作確認: none→TV=0/rank=1、low→0.0025、moderate→0.028、high→0.22。

## 6. Recoverabilityについて理論的に言えること

X = H Wᵀ(rank R_I のNMF)、H_nq = u_nᵀ(G ×₂ a_n) という構造から:

1. **[証明可能] A無変動なら復元不能**: 全行 a_n = a (C>1) のとき
   Σ_c a_c E_pcq = 0 なる摂動 E((C−1)R_N R_I 次元)が X を変えずに
   B を変える。cell type軸は完全に不定(interior G なら非負性も保持)。
2. **[証明可能・generic] R_N ≥ R_I では符号制約なしに非識別**:
   任意のgenericな代替core G' に対し G'×₂a_n が行フルランクとなり、
   各個人で u'_n が存在して H を厳密再現。非負性が救うかはsimulation課題。
   → R_N < R_I を主要regimeにすべき、という設計上の重要な帰結。
3. **[counting heuristic] N の閾値**: N(R_I − R_N) ≳ R_N(C·R_I − R_N)。
   (R_N,R_I,C)=(2,3,3) で N ≈ 14 に遷移予測 → 小gridの N=10–30 で観測可能。
4. **[既知理論の適用] W段はNMF**: X = HWᵀ の一意性はseparability /
   sufficiently scattered型条件に依存。実データでは検証不能、
   simulationでは「programが疎/相異なほど易しい」として現れるはず。

過剰な主張はせず、上記以外(非負性の効果、noise・A誤差・rare cell type・
rank誤指定の定量的影響)はすべてsimulation-onlyと明記した。

## 7. Simulationで確認すべきこと

- 非負性が符号なし非識別regime(特に R_N = R_I)を救うか
- 復元誤差 vs α₀(A variation)× N のrecoverability map、
  遷移点がcounting予測(N≈14)と一致するか
- rank(A)/condition/列相関/rare cell typeの影響
- X noise・A誤差(noise_A)への頑健性、rank誤指定の影響
- multi-start dispersion(X残差小×B分散大 = 非識別regimeの署名)
- Simulation 1: tensor不要の確認(matrix十分・tensorが悪化しないこと)
- Simulation 2: matrixはpopulation平均に潰れ、tensorが個人変動を復元するか
- 常に X 残差と B 誤差を両方報告(小残差・大B誤差が非識別の証拠)

## 8. 次にMUを導出してよい状態か

**まだ**。コード面の準備(generator・診断・baseline)は整ったが、
gateとして (i) 本設計(tensor-model.md / tensor-simulation.md)のレビュー、
(ii) 小規模grid(N=10–30, C=2–3, I=10–30, R_N=2, R_I∈{2,3,4})での
recoverability確認(isoTensor-experiments側で実行)が先。
特に理論的帰結 6-2 により、R_N < R_I regimeの確認を経ずに
MUを実装するのは危険。simulation結果がpartial Tuckerを正当化した時点で
MU導出(block-wise NNLS + MM)に進むのが妥当。

## 状態

Phase 2 fitting/MU未実装のまま停止。レビュー待ち:
`notes/tensor-model.md`(改訂版)と `notes/tensor-simulation.md`。
