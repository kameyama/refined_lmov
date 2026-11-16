# yui-kr の計算結果（Nakagane–Sano, arXiv:2111.00388）

- 取得元：https://github.com/taketo1024/yui-kr コミット 0a07fbf19b0ea0deedb17df4ebd9b1f91136d877（2024-11-26）、MIT ライセンス（LICENSE）。
- 取得日：2026-10-05。
- homology.json：`results/*.json` をまとめたもの。{結び目名: [[i, j, k, 階数], ...]}。reduced HOMFLY-PT ホモロジー。
  yui-kr の Poincaré 多項式は t^{(k−j)/2} q^i a^j（`src/lib/result.rs` の `poincare_poly`）。
- braid.csv：入力のブレイド（KnotInfo）。homfly.csv：KnotInfo の HOMFLY-PT（yui-kr の検証用）。この 2 つは KnotInfo のデータで、KnotInfo が再配布の条件を示していないので git では追跡しない（2026-11-16 から）。`python3 scripts/framed_knot/fetch_knot_table.py` で、上のコミットの `data/braid.csv`・`data/homfly.csv` を取得する（sha256 を照合する）。使った成果を発表するときは KnotInfo を引用する：C. Livingston and A. H. Moore, KnotInfo: Table of Knot Invariants, knotinfo.org（閲覧日を添える）。
- sage_homfly.txt：braid.csv のブレイドから Sage で計算した HOMFLY-PT（scripts/framed_knot/knot_table_homfly.sage の出力。このリポジトリで計算したもの）。
