# refined CS invariants of iterated torus knots

Kameyama–Nawata "Refined large N duality for knots"（arXiv:1703.05408）の refined BPS 不変量の計算を、
トーラス結び目からその**ケーブル**（iterated torus 結び目）に拡張する。
対象は |λ| ≤ 3 の全ヤング図形 [1], [2], [1,1], [3], [2,1], [1,1,1]。

詳細な仕様は [CLAUDE_TASK_refined_cable_knots.md](CLAUDE_TASK_refined_cable_knots.md)、
結果は [refined_cables/REPORT.md](refined_cables/REPORT.md)、作業経過は [refined_cables/LOG.md](refined_cables/LOG.md) を参照。

## セットアップ

```bash
just setup
```

```bash
just fetch-refs
```

## テスト

```bash
.venv/bin/python -m pytest refined_cables/tests -q
```
