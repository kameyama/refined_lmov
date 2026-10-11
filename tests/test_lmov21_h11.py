"""H11 の構成（scripts/lmov21_h11.py）の機能テスト：窓の規則 R ＋ framing の正値性 ＋ BPS 状態の総数の最小化。"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "scripts"))

import lmov21_h11 as H  # noqa: E402


def quiet(_):
    pass


def test_face_t25_unique_with_R_and_negative_framing():
    """T(2,5)：R ＋ framing −2..0 の spin・genus・既約表示の正値性 ＋ 総数最小で、面上の整数解は正解 1 個（約 15 秒）。"""
    res = H.face("T2_5", [-2, -1, 0], use_R=True, genus=True, irr=True, seeds=2, log=quiet)
    assert res["min_total"] == 0 and res["face_bounded"]
    assert res["face_count"] == 1 and res["face_count_contains_true"]


def test_face_t25_framing_zero_is_not_unique():
    """framing 0 だけ（R、spin のみ）では面上に複数の整数解が残る（framing が効くことの対照）。"""
    res = H.face("T2_5", [0], use_R=True, seeds=2, log=quiet)
    assert res["min_total"] == 0 and res["face_bounded"]
    assert res["face_count"] != 1


def test_predict_k1_reproduces_eha():
    """ツイスト結び目 K_1（三葉結び目）：GNSSS と X0 の入力から、実行可能解に自由度は無く、EHA の [2,1] に一致する。"""
    res = H.predict_twist(1, [-1, 0], log=quiet)
    assert res["face_count"] == 1 and res["matches_known"] == [True]


@pytest.mark.results
def test_face_t27_unique():
    """T(2,7)：R ＋ framing −6..0 ＋ 全ての正値性で一意（数分）。"""
    res = H.face("T2_7", list(range(-6, 1)), use_R=True, genus=True, irr=True, seeds=2, log=quiet)
    assert res["face_count"] == 1 and res["face_count_contains_true"]
