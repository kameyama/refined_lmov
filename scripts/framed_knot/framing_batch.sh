#!/bin/bash
# framed knot の refined LMOV の判定をまとめて実行する（scripts/framed_knot/framing_check.py を framing ごとに別プロセスで、1 件 600 秒まで）。
# 実行：bash scripts/framed_knot/framing_batch.sh > logs/framed_knot/framing_batch.log 2>&1。完了の目印は「framing batch done」の行。
R=${PYTHON:-/Users/masaya/projects/ai_research/refined_cs/.venv/bin/python}
export PYTHONPATH=src:scripts/framed_knot:scripts/dgr_lmov
run() {  # 結び目 XMIN XMAX NMAX
  for x in $(seq $2 $3); do
    s=$(date +%s)
    timeout 600 $R scripts/framed_knot/framing_check.py $1 $x $x $4 2>&1 || echo "$1 framing $x：失敗または時間切れ（終了コード $?）"
    echo "--- $1 framing $x：$(( $(date +%s) - s )) 秒"
  done
}
run unknot -6 6 3
run T2_3 -6 6 3
run 4_1 -4 4 3
run T2_5 -4 4 3
run T2_7 -3 3 3
run T3_4 -3 3 3
for p in -6 -5 -4 -3 -2 -1 1 2 3 4 5 6; do run twist:$p -4 4 2; done
run K1 -2 2 3
run K3 -2 2 3
echo "=== framing batch done $(date '+%F %T')"
