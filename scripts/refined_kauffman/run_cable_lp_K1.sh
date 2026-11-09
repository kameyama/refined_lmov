#!/bin/sh
# 続き 141：K1 の D_N の 2 箱を N = 11..13 で（基底の値を 1 個ずつ直列に、その後 Z を組み立てる）。
cd "$(dirname "$0")/../.."
export PYTHONPATH=src:scripts/refined_kauffman
for N in 11 12 13; do
  for nu in 4 31 22; do
    .venv/bin/python scripts/refined_kauffman/rk_c_cable_lp.py dnV $N $nu || exit 1
  done
  .venv/bin/python scripts/refined_kauffman/rk_c_cable_lp.py dn $N $N K1 || exit 1
done
echo "ALL DONE"
