#!/bin/bash
# framing を記号として扱う判定（scripts/framing_newton.py）：2 箱を起点のずらし 15 まで、間の framing の直接判定つきで。
# 完了の目印は「framing newton done」の行。3 箱は 6 根のニュートン形では証明書が見つからないので含めない（logs/framing_newton_first.log）。
R=${PYTHON:-/Users/masaya/projects/ai_research/refined_cs/.venv/bin/python}
export PYTHONPATH=src:scripts
for k in unknot T2_3 4_1 T2_5 T2_7 T3_4 K1 K3; do echo "=== $k 2 箱 $(date '+%T')"; timeout 7200 $R scripts/framing_newton.py $k 2 15; done
echo "=== framing newton done $(date '+%F %T')"
