#!/bin/sh
# 続き 141：K1 の D_N の 2 箱を N ごとに（基底の値を 1 個ずつ直列に、その後 Z を組み立てる）。使い方：run_cable_lp_K1.sh N...（既定 11 12 13）
cd "$(dirname "$0")/../.."
export PYTHONPATH=src:scripts/refined_kauffman
NS=${*:-"11 12 13"}
for N in $NS; do
  for nu in 4 31 22; do
    .venv/bin/python scripts/refined_kauffman/rk_c_cable_lp.py dnV $N $nu || exit 1
  done
  .venv/bin/python scripts/refined_kauffman/rk_c_cable_lp.py dn $N $N K1 || exit 1
done
echo "ALL DONE"
