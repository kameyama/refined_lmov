#!/bin/bash
# framed refined LMOV の判定を (結び目, n, framing の範囲) ごとにまとめて実行する（scripts/framed_knot/framing_big.py、保存済みの単位は飛ばす）。
# 使い方：bash scripts/framed_knot/framing_big_batch.sh WORKERS 結び目 n XMIN XMAX [結び目 n XMIN XMAX ...]
# 完了の目印は「framing big batch done」の行。
cd "$(dirname "$0")/.." || exit 1
export PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov
W=$1; shift
while [ $# -ge 4 ]; do
  /usr/bin/time -l .venv/bin/python scripts/framed_knot/framing_big.py --workers "$W" "$1" "$2" "$(echo $(seq -- "$3" "$4") | tr " " ",")" 2>&1 \
    | awk '/axis/ {next} /peak memory|real/ {print; next} /^ *[0-9]+  +[a-z]/ {next} {print}'  # time -l の出力は実時間とピークメモリだけ残す
  shift 4
done
echo "=== framing big batch done $(date '+%F %T')"
