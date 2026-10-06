"""open_problems D1：正規化 (N1)(N2) で決めた単項式 M_□ の閉じた式。

M_□ = (−q/t)^{b−1}·(a (q/t)^{1/2})^{δ}。b はブレイドの本数（Π p_i）、δ は反復ごとに δ ↦ p δ + (p−1)(q−1)/2（トーラス結び目 T(m,n) は δ = (m−1)(n−1)/2）。
正のブレイド表示では δ = (w − b + 1)/2（w は writhe）で、Markov 跡による正規化の refined 版と読める。保存済みの全トーラス結び目と全ケーブルで照合する。
"""
import json
from fractions import Fraction
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent


def steps_of(d):
    if "torus" in d:
        return [tuple(d["torus"])]
    if "steps" in d:
        return [tuple(s) for s in d["steps"]]
    return [tuple(d["cable"]["companion"]), tuple(d["cable"]["pattern"])]


def predicted(steps):
    b, delta = 1, Fraction(0)
    for p, q in steps:
        b, delta = b * p, p * delta + Fraction((p - 1) * (q - 1), 2)
    # M = sign · a^δ · (q/t)^{δ/2 + b − 1}
    return (-1) ** (b - 1), delta, delta / 2 + b - 1


FILES = sorted(ROOT.glob("data/torus/*_1.json")) + sorted(ROOT.glob("data/cables/K*_1.json"))


@pytest.mark.parametrize("path", FILES, ids=[f.stem for f in FILES])
def test_box_monomial_formula(path):
    d = json.loads(path.read_text())
    m = d["monomial_applied"]
    sign, a_exp, qt_exp = predicted(steps_of(d))
    assert m["sign"] == sign
    assert Fraction(m["a"]) == a_exp
    assert Fraction(m["q"]) == qt_exp and Fraction(m["t"]) == -qt_exp
