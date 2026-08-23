# Identifiability検証報告 (2026-08-23)

対象: 「R_N ≥ R_I では符号制約なしにgenericに復元不能」claimの厳密検証。
実装変更なし(package API不変)。変更: `notes/tensor-model.md`(identifiability
節を全面改訂)、`notes/tensor-simulation.md`(grid修正)、
`dev/identifiability_check.R`(開発用検証スクリプト、build対象外)。

## 1. R_N ≥ R_I claimは正しかったか

**正しかった。しかも従来の分類(LOCAL/GENERIC EVIDENCE)から
PROVEDに格上げでき、さらに強い形で成立する。**

証明(構成的): x_nᵀ = u_nᵀ K_n Wᵀ, K_n = G ×₂ a_n。R_N ≥ R_I なら
**任意の**genericな代替core G' に対し K'_n = G'×₂a_n はrank R_I
(行空間 = R^{R_I})なので、各nで u'_nᵀK'_n = u_nᵀK_n が解ける。
(U', G', W) は X を厳密再現し、generic G' では B' ≠ B。
countingではなく明示的な解の構成であり、C ≥ 2、a_n ≠ 0、W列フルランク
以外に仮定は不要。**特に A のvariationには一切依存しない**(どれだけ
Aが変動しても救えない)。

強化(Theorem A′): 真値が内点(U > 0, G > 0)なら G' = G + εE、
u'_n → u_n (ε→0) により**非負制約下でも**連続族が存在する。
数値実証: (R_N,R_I,C,N)=(2,2,3,20)、ε=1.25e-2 で全factor正、
‖X'−X‖/‖X‖ = 1.4e-16、‖B'−B‖/‖B‖ = 1.9e-2。

## 2. 旧derivationの論理の飛躍

- 旧記述は「X = HWᵀ の分解を所与として」という2段階の迂回をしており、
  W段の解決を仮定しているように読めた。実際には同じWを使う直接構成で
  足り、仮定は不要だった(結論は変わらず、証明が簡潔・強くなった)。
- 旧記述は「非負性が救うかはsimulation課題」と全面的に開いていたが、
  内点真値については**救わないことが証明できる**(制約が非活性のため)。
  未解決なのはsparse/boundary真値の場合のみ。
- 旧countingの閾値 N* は、Jacobian解析で**局所deficiencyまで厳密一致**
  する必要条件であることを確認(下記5)。

## 3. Latent tensor B のidentifiabilityについて確実に言えること

3水準(A: factor / B: latent tensor / C: downstream)を区別した上で、
B水準について(notes/tensor-model.md §4, §9):

- **Theorem A [PROVED]**: R_N ≥ R_I, C ≥ 2 ⇒ B非識別(実数値factor、
  構成的、Aに依らない)。
- **Theorem A′ [PROVED, generic内点真値]**: 上記は strictly positive
  な真値なら非負制約下でも成立。⇒ dense真値では R_N < R_I が必要。
- **Theorem B [PROVED]**: rank(A) = C は必要条件。rank(A) = r < C なら
  (C−r)R_N R_I 次元の不可視方向が存在(同一行 r=1 が最悪)。
  十分条件ではない(Theorem Aより)。
- **Theorem C [NECESSARY CONDITION]**: R_N < R_I のとき
  N(R_I−R_N) ≥ R_N(C·R_I−R_N)(W既知条件付きのfiber-dimension論法;
  Case 2の非識別はCase 4に遺伝するのでfull problemにも必要)。
- **Case 1 [PROVED]**: U, W既知なら線形逆問題で、識別 ⟺
  rank(W)=R_I かつ rank(U⊙_r A)=R_N C(⊙_r = 行方向Khatri–Rao)。
- **Case 3 [PROVED, 十分方向]**: U既知なら rank(U⊙_r A)=R_N C で
  Bは線形に識別可能(R_N vs R_I 無関係)。U を知ることは問題を解決し、
  W を知ることは解決しない — 困難は U と A が index n を共有する
  縮約に集中している。
- 手計算例(R_N=R_I=1): B_nci = u_n g_c w_i で任意の g' に対し
  u'_n = u_n(a_nᵀg)/(a_nᵀg') が同じXを与えB'≠B — 最小スケールで確認。

数値(exact Jacobian、defic = rank J_B − rank J_X):
Theorem A系列は構成族の次元 R_N C R_I + N(R_N−R_I) − R_N² と
**全ケース厳密一致**(1,6,4,8,13)。旧gridのdefault
(N=20,C=3,R_N=R_I=2)もdefic=8で非識別。局所フルランク ≠ 大域一意
(defic=0はLOCAL/GENERIC EVIDENCEに留まる)。

## 4. A に必要そうな条件

- rank(A) = C [PROVED必要]、行間variation [PROVED必要(同一行はr=1)]。
- ただし**full rank + 高variationでも不十分**(Theorem Aが任意のAで成立)。
- 実効的に効くのは A 単独ではなく {u_n ⊗ a_n} の張り方
  (rank(U⊙_r A) ≤ min(N, R_N·rank(A)))。
- 列相関大: 識別可能でも悪条件化 → noise増幅(識別性と適切性は別;
  EMPIRICAL)。rare cell type: 列が厳密0なら該当sliceは完全に不可視
  [PROVED極限]; 近似0では誤差増幅 [EMPIRICAL]。

## 5. Rankに関して確実に言えること

- dense(内点)真値では **R_N < R_I が必要** [PROVED]。
- さらに N ≥ N* = ⌈R_N(C·R_I−R_N)/(R_I−R_N)⌉ が必要型条件
  [NECESSARY CONDITION]。数値では閾値未満のdeficが予測式
  R_N(C·R_I−R_N) − N(R_I−R_N) と厳密一致((1,2,C2): N*=3、
  (2,3,C2): N*=8、(2,3,C3): N*=14)、閾値以上でdefic=0。
- 十分条件は未証明: R_N < R_I + N ≥ N* + generic での大域一意性は
  EMPIRICAL HYPOTHESIS(multi-startで検証)。
- sparse真値での非負性による回復(R_N ≥ R_I を含む)は
  EMPIRICAL HYPOTHESIS。

## 6. Simulation gridの修正

「主に R_N < R_I に限定」という制約を撤廃し、
- R_N <, =, > R_I の3 regimeすべてをsmall gridで試す
- dense真値 × R_N ≥ R_I は**negative control**として活用
  (X残差≈0でB誤差大/multi-start分散が出なければpipelineのバグ)
- **factor sparsityを新しいgrid軸に追加**(非負性が効き得る唯一の機構。
  generatorへのsparsity引数はgrid開始時に追加予定、当面は手動構成)
- R_N < R_I では N を N* 跨ぎでsweepし、大域的・経験的な境界が
  局所理論の境界を追うか確認

## 7. MU実装前に残る理論上の問題

1. R_N < R_I regimeの**大域**一意性(局所証拠のみ。multi-start
   simulationが代替)
2. sparse/boundary真値での非負性の効果(separability型条件の要否)
3. W段(X = HWᵀ)のNMF一意性条件が実データで満たされる見込み
4. 識別可能regime内でのwell-posedness(κ(A)、列相関、rare cell type、
   noise、A誤差の定量的影響)
5. アルゴリズム問題(MUが識別された解に収束するか)は識別性とは別で、
   fitting実装後の課題

## 状態

Phase 2 fitting/MUは未実装のまま停止。レビュー対象:
`notes/tensor-model.md`(identifiability全面改訂版)、
`notes/tensor-simulation.md`(grid修正)、`dev/identifiability_check.R`。
