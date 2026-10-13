# yui-kr の計算結果（Nakagane–Sano, arXiv:2111.00388）

- 取得元：https://github.com/taketo1024/yui-kr コミット 0a07fbf19b0ea0deedb17df4ebd9b1f91136d877（2024-11-26）、MIT ライセンス（LICENSE）。
- 取得日：2026-10-05。
- homology.json：`results/*.json` をまとめたもの。{結び目名: [[i, j, k, 階数], ...]}。reduced HOMFLY-PT ホモロジー。
  yui-kr の Poincaré 多項式は t^{(k−j)/2} q^i a^j（`src/lib/result.rs` の `poincare_poly`）。
- braid.csv：入力のブレイド（KnotInfo）。homfly.csv：KnotInfo の HOMFLY-PT（yui-kr の検証用）。
