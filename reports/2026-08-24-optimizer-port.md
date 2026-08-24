# L-BFGS optimizer移植報告 (2026-08-24)

isoTensor-experimentsで検証済みのlog-reparameterized L-BFGSを
package本体へ移植した。変更: `notes/tensor-LBFGS.md`(新規、実装前に
数式を整理)、`R/isoTensor.R`、`NAMESPACE`(stats::optim)、
`tests/testthat/test_isoTensor_lbfgs.R`(新規)、
`test_isoTensor_tensor.R`(names更新・MU明示化・thr=0 test)、
`dev/lbfgs_port_check.R`(新規、移植検証)、`man/isoTensor.Rd`、
`README.md`、`vignettes/isoTensor.Rmd`(tensor節追加)、`inst/NEWS`、
`PLAN.md`、`notes/tensor-MU.md`(位置づけ更新)。

## 1. 採用したoptimizer API

```r
isoTensor(X, A, model="tensor",
    rank_individual=..., rank_isoform=...,
    optimizer=c("lbfgs", "hybrid", "mu"),
    optimizer.control=list(maxit=5000, factr=10))
```

単一enum引数(experiments報告Part Kの推奨形)。`algorithm`(目的関数、
現状Frobeniusのみ、返り値に記録)と`optimizer`(解法)は概念的に分離し、
将来のKL/beta-divergence追加時は`algorithm`引数を足す設計
(今回は引数を追加しない)。`num.iter`/`thr`/`pseudocount`はMU段
(mu・hybridのstage 1)のみに作用、L-BFGS段は`optimizer.control`
(maxit, factr — experimentsの較正値がdefault)。matrix modeで
optimizerを明示指定した場合はwarning付きで無視(matrixは常にMU)。

## 2. Default optimizer = "lbfgs"(理由)

- 0.99.x未リリースでdefault変更のコストが最小の時期
- full recoverability grid(15 datasets×30 inits)でL-BFGS単独は
  hybridと同等以上・やや高速、realistic noise stressでは全datasetで
  selected relBが一致しL-BFGS単独が高速
- calibrationでMU 20kの約1/10の時間でrelBを5–100倍改善
- 完全な旧挙動は`optimizer="mu"`で維持(コードは無変更・削除なし)

## 3. L-BFGS数式(notes/tensor-LBFGS.md)

U=exp(Θ_U), G=exp(Θ_G), W=exp(Θ_W)で非負制約を除去し、
L=½‖X−X̂‖²_Fを`optim(method="L-BFGS-B")`(bounds無し=plain L-BFGS、
experimentsと同一)で最小化。勾配は既存の解析gradient
(`.gradIsoTensor`、FD検証済み)にchain ruleで因子をHadamard積:
∂L/∂Θ_U = U⊙∇_U L(G, Wも同様)。packingは
θ=(vec Θ_U, vec Θ_G, vec Θ_W)(column-major、Gはp最速→c→q)、
pack時に`log(pmax(x, 1e-10))`。pack→unpack完全一致・θでの
gradient FD一致はunit test化。

## 4. Numerical stability対策

experiments実装に忠実(clipping・safe-exp・目的関数を変える処理なし):
overflowはline searchが拒否、underflowは境界への正しい極限として無害。
まれなoptim例外(experiments実測0.4–0.9%)はtryCatchで処理 —
lbfgs単独は情報付きerror(別seed/hybrid/muを提案)、hybridは
warning付きでMU段の解を返す。scaling ambiguityへの正規化は
optim中には行わず(gradient consistencyを壊すため)、終了後に1回のみ
`.normalizeFactorsTensor`(B・X̂不変、test済み)。lbfgsのrandom init
はMUと同一分布・同一draw順(experimentsと数値一致させるため、MU init
の global rescaleは適用しない)。init中のzeroはlbfgs系では1e-10へ
引き上げ(lockされない)、MUではlockのまま — optimizer別に文書化。

## 5. Experiments実装との一致(dev/lbfgs_port_check.R)

- **bitwise一致**: experimentsの`src/lbfgs_utils.R`を同一環境で実行し、
  3較正datasets(easy_high/mod_a20/hard_low)×5 seedsの全15 runで
  objectiveが全桁一致、Bのmax相対差 1.4e-16〜3.3e-16(機械精度)。
  「似たアルゴリズム」ではなく同一実装であることを数値確認。
- archived TSV(container内で生成、別BLAS/R build)とはL-BFGS軌跡の
  fp発散により per-run一致はしないが、median relBは同水準
  (easy 2.6e-5 vs 1.5e-4 / mod 1.19e-2 vs 1.33e-2 / hard 6.9e-2 vs 4.4e-2)。
- 付随して判明: calibration A6/A7のhybrid行はwarm段が`thr=0`
  (当時のsentinel仕様でMU 0回)だったため実質「rescaled init からの
  L-BFGS」(runtime 0.9sとも整合)。full grid・noise stressのhybrid
  (thr=1e-8)は正常で、結論(lbfgs≧hybrid)は不変。

## 6. mu / lbfgs / hybrid比較(easy_high, seed 1)

| optimizer | objective | relX | relB | runtime |
|---|---|---|---|---|
| mu (20k iter, thr=1e-8) | 1.3e-4 | 6.5e-4 | 7.7e-2 | 5.0 s |
| lbfgs (default) | 1.4e-8 | 6.6e-6 | **2.6e-5** | **0.54 s** |
| hybrid (MU 2k → L-BFGS) | 6.1e-13 | 4.4e-8 | 3.6e-6 | 0.58 s |

experimentsの「L-BFGSはMUより高速かつ高精度」がpackage本体でも再現。
unit testでもlbfgs/hybridのobjective ≤ muを同一easy dataで確認。

## 7. thr=0修正

初期sentinelを`RelChange[1] <- max(thr*10, .Machine$double.xmin)`に
変更(matrix・tensor MU両loop)。thr=0で`while(RelChange > 0)`が初回
から成立し、指定num.iterを全て実行(thr>0の挙動は不変:
max(thr*10, xmin)=thr*10)。unit test: matrix(RecError長=num.iter+1)
とtensor(NumIter=num.iter)の両方。

## 8. Unit tests

`test_isoTensor_lbfgs.R`(新規): default=lbfgs確認、返り値構造、
次元、非負性、seed再現性、objective減少、θでのgradient FD一致、
pack/unpack roundtrip+次元整合、user init、zero含みinitの仕様
(1e-10へ引き上げ・可動)、convergence metadata、optimizer.control
検証(不正名error・mu時warning)、hybrid(refined ≤ MU final、構造、
非負性)、cross-optimizer(lbfgs/hybrid ≤ mu objective、relB<0.05)、
A不変、warning系がlbfgs経路でも発火。既存MU testsは全て維持
(MU固有testはoptimizer="mu"明示、names更新)+thr=0 test+
matrix mode optimizer無視warning test。

## 9. R CMD check

`R CMD check --no-manual`: **Status: OK(0 errors / 0 warnings /
0 notes)**。vignette(tensor節込み)・examples・全testsを含む。

## 10. Backward compatibility

- matrix mode: 完全不変(返り値・数値・names、regression test付き)
- tensor mode: `optimizer="mu"`で従来コードpathを完全維持(削除なし。
  単調性が証明された唯一のoptimizerとしてbaseline/デバッグ用に保持)
- 破壊的変更は2点のみ(0.99.x未リリースのため許容と判断):
  (i) tensor modeのdefault optimizerがMU→lbfgs(返り値の
  RelChange/NumIter/ConvergedはMU系のみに存在)、
  (ii) 全tensor返り値に`optimizer`/`objective` fieldを追加
  (names完全一致に依存するコードは要修正)。
  isoTensor-experiments側は`optimizer="mu"`の明示追加で従来と同一。
- thr=0の挙動変更は「1回も回らない」→「num.iter回実行」で、
  experiments側は既に1e-300で回避済みのため影響なし

## 11. 残課題

1. multi-start package API(n.start; 当面experiments側の責務)
2. rank selection、KL/beta-divergence(algorithm引数の実体化)、
   spillover対策、A uncertainty(いずれも未着手のまま)
3. hybridのRecErrorは「MU履歴+refined 1点」— L-BFGS段の反復履歴は
   optim()が返さないため取得不能(偽値は作らない方針)
4. lbfgsのまれなoptim例外(<1%)はerror/fallbackで処理するが、
   自動リトライ(別seed)は未実装
5. MU c-loop未最適化(正しさ優先、従来どおり)

## 状態

optimizer移植のみ完了、ここで停止。Sciege再現・GTEx/ROSMAP・
real-data・rank selectionには進んでいない。
