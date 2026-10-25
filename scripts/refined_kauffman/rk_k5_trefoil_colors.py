r"""(b)：trefoil の refined 不変量（KN の U(N) refined CS、a 一般、framing 0 の (N1)+(N2) 正規化）を、
SO(6) = SU(4) の Λ^k V に必要な彩色について作る：Λ¹V = [1,1]、Λ²V = [2,1,1]、Λ³V = [2] ⊕ [2,2,2]（SO(4) = SU(2)² は [1], [2]）。
既存の scripts/stage2_torus.py（GN を有限体上で評価して厳密に再構成、t = q を framing 0 の HOMFLY に合わせる）をそのまま使い、
保存先だけ data/refined_kauffman/torus に変える。"""
import sys
from pathlib import Path
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman"); sys.path.insert(1, "scripts/framed_knot")
import stage2_torus as S

S.OUT = Path("data/refined_kauffman/torus")
S.COLORS.update({"211": (2, 1, 1), "222": (2, 2, 2), "2": (2,), "11": (1, 1), "1": (1,), "2111": (2, 1, 1, 1), "33": (3, 3), "321": (3, 2, 1), "31": (3, 1), "3111": (3, 1, 1, 1), "2211": (2, 2, 1, 1), "1111": (1, 1, 1, 1), "2222": (2, 2, 2, 2), "11111": (1, 1, 1, 1, 1), "3": (3,), "21": (2, 1), "111": (1, 1, 1), "22": (2, 2)})
import os
KNOT = os.environ.get("KNOT", "T2_3")      # 例：KNOT=T2_5
for c in sys.argv[1:] or ["1", "2", "11", "211", "222"]:
    norm, meta = S.compute(KNOT, c)
    print(c, meta["monomial_applied"], meta["seconds"], "s", len(norm), "項", flush=True)
