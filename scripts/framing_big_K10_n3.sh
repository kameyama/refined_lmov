#!/bin/bash
# K10 の 3 箱の framed refined LMOV の判定（framing −2..2、x = 0 は framing 0 の既存の判定 data/lmov/K10_3 と同じなので除く）。1 件約 1 時間、親プロセス約 9 GB。
# 実行：nohup bash scripts/framing_big_K10_n3.sh > logs/framing_big_K10_n3.log 2>&1 &。完了の目印は「K10 n=3 終了」の行。保存済みの framing は飛ばす。
cd "$(dirname "$0")/.." || exit 1
export PYTHONPATH=src:scripts
echo "=== K10 n=3 開始 $(date '+%F %T')"
.venv/bin/python -u scripts/framing_big.py --workers 4 K10 3 -1,-2,2 2>&1 | grep --line-buffered -v axis
echo "=== K10 n=3 終了 $(date '+%F %T')"
