"""トーラス結び目の箱数を 4 から順に増やし、1 件の計算時間が 10 分を超えない範囲で止める（ユーザーの指示、2026-10-10）。

次の箱数の時間は、直前の 2 件の比（無ければ 5 倍、4 倍未満なら 4 倍）で見積もる。見積もりが LIMIT を超えたら、その結び目はそこで止める。
1 件が実際に LIMIT を超えた場合も記録して止める。使い方：python scripts/highbox_batch.py [NMAX]
ログは logs/highbox_torus.log（各件）と logs/highbox_batch.log（打ち切りの判断）。
"""
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import highbox_torus as H  # noqa: E402

LIMIT = 600
KNOTS = [(2, 3), (2, 5), (3, 4), (2, 7), (3, 5), (2, 9), (2, 11), (3, 7), (4, 5)]


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 9
    with open("logs/highbox_batch.log", "a") as lg:
        def log(s):
            print(s, flush=True)
            lg.write(s + "\n")
            lg.flush()
        log(f"== batch start {time.strftime('%Y-%m-%d %H:%M:%S')} NMAX={nmax} LIMIT={LIMIT}s")
        for m, n in KNOTS:
            times = []
            for N in range(4, nmax + 1):
                if times:
                    ratio = times[-1] / times[-2] if len(times) > 1 and times[-2] > 5 else 5
                    est = times[-1] * max(ratio, 4)          # 比が小さく出ても 4 倍以上と見る（T(2,9) の 8 箱で 2.5 倍の見積もりが外れた）
                    if est > LIMIT:
                        log(f"T{m}_{n}: N={N} の見積もり {est:.0f}s > {LIMIT}s のため N={N - 1} で止める")
                        break
                sec = H.compute(m, n, N)
                times.append(sec)
                log(f"T{m}_{n} N={N}: {sec:.0f}s")
                if sec > LIMIT:
                    log(f"T{m}_{n}: N={N} が {LIMIT}s を超えたのでここで止める")
                    break
        log(f"batch done {time.strftime('%Y-%m-%d %H:%M:%S')}")


if __name__ == "__main__":
    main()
