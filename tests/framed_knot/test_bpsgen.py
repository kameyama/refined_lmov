"""任意の箱数の \\hat f（rcsknots.bpsgen、scripts/framed_knot/highbox_torus.py）の機能テスト。"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp

from rcsknots.bps import M_matrix, T as ST
from rcsknots.bpsgen import B_class, EpsTable
from rcsknots.characters import chi
from rcsknots.field import V, P31, make_conv
from rcsknots.gn import GN
from rcsknots.partitions import partitions

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts/framed_knot"))
import highbox_torus as H  # noqa: E402


def test_B_class_diagonalizes_M():
    """M(t) は B(t) との Kronecker 積なので、指標の基底で対角：Σ_μ M_{μρ} χ_μ(κ) = B(κ) χ_ρ(κ)。"""
    for n in (2, 3, 4):
        M = M_matrix(n, ST)
        parts = partitions(n)
        for kappa in parts:
            for j, rho in enumerate(parts):
                lhs = sum(M[i, j] * chi(mu, kappa) for i, mu in enumerate(parts))
                rhs = B_class(kappa, ST) * chi(rho, kappa)
                assert sp.simplify(lhs - rhs) == 0


def test_eps_table_matches_evaluate_apoly():
    p, k = P31, 4
    rng = np.random.default_rng(3)
    G = GN(V(rng.integers(2, 10**6, k), p), V(rng.integers(2, 10**6, k), p), V(np.full(k, 3), p), make_conv(p, k))
    vec = G.op(2, 2, -3, G.op(1, 2, -3, {(): G.one}))
    a, b = EpsTable(G).evaluate(vec), G.evaluate_apoly(vec)
    assert len(a) == len(b) and all(x == y for x, y in zip(a, b))


def test_reproduces_saved_hat_f_up_to_3_boxes():
    """新しい経路（GN × §13 の正規化式 → 評価点ごとの plethystic log → 補間）が保存済みの \\hat f（第 1・3 期）を全係数で再現する。"""
    for m, n in ((2, 3), (3, 4)):
        for N in (1, 2, 3):
            hats, _ = H.run(m, n, N, workers=1, log=lambda s: None)
            for rho, h in hats.items():
                d = json.loads(Path(f"data/framed_knot/bps/T{m}_{n}_{''.join(map(str, rho))}.json").read_text())
                assert h == {tuple(map(int, k.split(","))): v for k, v in d["hat_f"].items()}


def test_trefoil_4_boxes_positive():
    hats, _ = H.run(2, 3, 4, workers=1, log=lambda s: None)
    for rho, h in hats.items():
        r, Nh = H.judge(h, 4)
        assert r["before_neg"] == 0 and r["after_neg"] == 0
