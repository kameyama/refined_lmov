"""結び目表のブレイドと HOMFLY-PT（KnotInfo のデータ）を取得する。

KnotInfo（C. Livingston and A. H. Moore, KnotInfo: Table of Knot Invariants, knotinfo.org）は
データの再配布の条件を示していないので、このリポジトリには置かない。yui-kr のリポジトリ
（https://github.com/taketo1024/yui-kr、MIT）が入力として同梱しているものを、コミットを固定して取得する。

実行：python3 scripts/framed_knot/fetch_knot_table.py（raw.githubusercontent.com に接続する）
出力：data/framed_knot/reference/yui_kr/braid.csv, homfly.csv（どちらも git では追跡しない）
このデータを使った成果を発表するときは、KnotInfo を引用する（書式は NOTICE.md）。
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

COMMIT = "0a07fbf19b0ea0deedb17df4ebd9b1f91136d877"  # 2024-11-26
BASE = f"https://raw.githubusercontent.com/taketo1024/yui-kr/{COMMIT}/data/"
SHA256 = {
    "braid.csv": "400cf8b53742146915d3f814aa6b6df3323d9e48c07c77f3ef7eae81647134eb",
    "homfly.csv": "a3580a255cb2655697feb160f31f52561edd7709d3d8c064d8b30cf054dbf81e",
}
DEST = Path(__file__).resolve().parents[2] / "data/framed_knot/reference/yui_kr"


def main(dest=DEST):
    for name, expected in SHA256.items():
        data = urllib.request.urlopen(BASE + name).read()
        got = hashlib.sha256(data).hexdigest()
        if got != expected:
            sys.exit(f"{name}: sha256 が記録と違う（{got}）。取得元が変わっていないか確かめる。")
        (Path(dest) / name).write_bytes(data)
        print(f"{name}: {len(data)} バイト")


if __name__ == "__main__":
    main(*sys.argv[1:2])
