# 研究計画：refined Kauffman 多項式の large N duality（第 3 の研究、2026-10-14 作成）

規則は README.md。文献は literature.md、2018 年の試行の整理は past_attempt.md。

## 0. 問い
- Q1 unrefined の Kauffman の large N duality（Mariño 0904.1088）の構造は何か。自分のコードで再現できるか。【K0、済】
- Q2 2018 年の試行はどこで、なぜ破れたか。【K1、済（past_attempt.md）】
- Q3 refined では unknot の時点で何が起きるか。AS の SO(2N) refined unknot（D_N Macdonald 次元）と U(N) の composite の Macdonald 次元（どちらも閉じた式）から、refined の orientifold の開弦の構造（½ と Σ_{d 奇数} の refined 版、基底、極の打ち消し）を決められるか。【K2】
- Q4 その構造の下で、結び目（trefoil, T(2,5), 4_1）の refined Kauffman 不変量（AS の SO(2N) refined CS、または NRZ/GW の Kauffman ホモロジー）と refined composite 不変量（Cherednik–Elliot の composite DAHA 超多項式、または U(N) refined CS の composite 彩色）から抽出した \hat g_ρ, \hat h_ρ は整数係数の Laurent 多項式か。正値性は何に課せるか。【K3〜K5】

## 1. 段階

### K0：unrefined の再現【済、2026-10-14】
- Mariño の構造（g_R, h_R, M_{RS}、予想）を literature.md §1 に整理。scripts/rk_unrefined_check.py で trefoil の \hat g_[2], \hat g_[1,1]（Mariño (4.8)）と 4_1 の \hat h, \hat g（NRZ 付録 C）を全係数で再現。誤植 2 点（Mariño (4.8) の z^9、NRZ (C.9)）を特定。
- 残り：tests/ に移す（test_rk_unrefined.py）。

### K1：2018 年の試行の監査【済（草稿）、2026-10-14】
- past_attempt.md。構造上のミス (A)〜(F)、実装上の不整合 (a)〜(j)。

### K2：refined unknot から構造を決める（本計算に入る前の停止点）
入力（すべて閉じた式、記号計算）：
- Kauffman 側：AS の SO(2N) refined unknot dim_{q,t}^{D_N}(λ)（scripts/rk_refined_unknot_poles.py の eval_DN）。|λ| ≤ 3（必要なら 4）。
- composite 側：U(N) の composite 表現 [λ,μ]_N の A 型 Macdonald 次元（(q,t)-hook の式、N を記号として扱えないので N = 4..8 で計算して a = t^N の有理式に直す）。
- 変数の同定：両者を共通の Kähler 変数で表す（unrefined では Mariño の ν。refined では AS の Q = t^{2N−1} と U(N) の a = t^{N_U} の関係を unknot の [1] で決める）。
検討する仮説（優先順）：
1. 基底：Kauffman 側・composite 側とも KN の t-brane 基底 P_{λ^T}(−x;t,q)（2018 年の最終形）か、composite の片側を q-brane にするか。
2. 展開：Σ_d と Σ_{d 奇数} の refined 版。AS §7 の閉弦の構造（oriented 部分が (q, t²) の形、偶数 d の引き算に (1 + q^d)/(1 + t^d)）を開弦に持ち上げた形を第一候補にする。
3. 極の打ち消し：D_N 次元の分母 (1 − q t^{N−1}) 等と composite の分母 (1 − q t^{N_U−1}) 等が、どの同定と組み合わせで打ち消すかを 2 箱で調べる。
判定：unknot の \hat g, \hat h が（unrefined の \hat g_[1] = 1、\hat h_[2] = \hat h_[1,1] = −z^{−1}、他 0 を refined 化した形で）Laurent 多項式になる組み合わせがあるか。
停止点：K2 の結果（できた／できない、どの仮説が残るか）を報告して確認を求める。

### K2'：K2 の否定的な結論を受けた再定式化（2026-10-14 追加、deep research の報告に基づく）
- 判定基準を第二量子化に移す：AS (3.13) に倣い、Kauffman 側の母関数（unknot）が Π(1 − q^• t^• Q^•)^{整数} の形（向き付け部分は t の偶数冪、固定面の部分は t の昇冪なし）に書けるかを判定する。½ と Σ_{d 奇数} は使わない。
- 基底：A 型 Macdonald の代わりに lifted Koornwinder \tilde K_λ(x;q,t,T;t̄)（T = Kähler 変数）と Rains の lifted Cauchy 恒等式を使う。D_N の分母が C⁺ 因子（T の極）であること【計算】を利用して、分母を基底と計量に吸収できるかを調べる。
- Kähler 変数のずれ：交差帽の数ごとに Q を量子シフトする（Avetisyan–Mkrtchyan 2205.12832 の閉弦の三重化の開弦版）【仮説】。
- 境界パラメータ：(±1, ±q^{1/2})（AS）と (±1, ±t^{1/2})（LRW）、(±q^{1/2}, ±t^{1/2})（LRW Conjecture 1.4）を unknot の 2〜3 箱で比べる。
- 停止点：unknot の 2 箱で積表示が得られるか否かが決まったとき。
- K2'(a)（2026-10-14 深夜、実施済み）：probe 側の Askey–Wilson 係数を U = t^N の閉じた形にし、virtual 正規化 y^N z^N Z で Hopf 絡み目の log Z を Q(q,t,U) の上で計算、Möbius 反転した f_{ab} の分母と整数性を見る。
- K2'(b)（2026-10-14、実施済み）：成分 1 の probe を無限個にしても S 行列の対称性で 1 変数の \hat K_{(k)}(z) と冪和に縮約できる（多変数 Koornwinder 不要）。|R| ≤ 4、b ≤ 2 で Schur 基底の整数性を確認（report.md §4.2）。Mariño の \hat g_{[2],[1]} との係数ごとの対応は未実施。
- K3'（次、trefoil）：AS 1202.4456 (6.15), (6.17) の T(2,2m+1) の [1] を入力に、f_1 = −W_{[1]} の構造（refined の h/g 分解：½h^{ref}_[1] = A 型 refined CS の [1] に A = U²/t を入れたもの、g^{ref}_[1] = f_1 − ½h^{ref}_[1] が Laurent 多項式か＝refined Rudolph 性）を調べる。framing は refined Rosso–Jones（p_2 の Koornwinder 展開）で AS の式と合わせる。

### K3：結び目の refined 入力の整備
- refined Kauffman：AS の SO(2N) refined CS の S, T 行列（D_N Macdonald 多項式、小さい N で数値的に）から trefoil の [1], [2], [1,1]。NRZ の対称表現の Kauffman ホモロジー（trefoil, 4_1）との関係（AS は「トーラス結び目では変数変換で結べない」）を記録。
- refined composite：Cherednik–Elliot の composite DAHA 超多項式（trefoil の ([1],[1]) は本文に明示）。unreduced 化に必要な次元は K2 の結果に従う。
- refined HOMFLY：既存（src/rcsknots）。

### K4：結び目での検証
- trefoil（2 箱）で \hat g, \hat h の整数性。破れるなら、どの入力の組み合わせでも破れるかを分類する。
- 正値性は [1] から（CLAUDE.md の規則どおり、G 展開の前後）。

### K5：報告
- report.md に結果と導出、novelty を literature.md §3 に。

## 2. 停止点
- K2 の後（refined unknot の構造）。
- research_plan の各段階で 1 件 30 分を超える計算が必要になったとき（CLAUDE.md の準備の規則）。

## 3. 状態
| 段階 | 状態 |
|---|---|
| K0 | 済（2026-10-14） |
| K1 | 済・草稿（2026-10-14） |
| K2 | 済・否定的（2026-10-14、progress.md） |
| K2' | 済・条件付き（2026-10-14、report.md §4–5） |
| K2'(a) | 済・肯定的（2026-10-14、report.md §4.1）：Hopf 絡み目（probe 各 1 個、a,b ≤ 4）で f_{ab} の整数性。次：probe 2 個、trefoil |
| K2'(b) | 済・肯定的（2026-10-14、report.md §4.2）：probe 無限個、Schur 基底で |R| ≤ 4、b ≤ 2 の整数性。Macdonald 基底では破れる |
| K3〜K5 | 未着手 |
