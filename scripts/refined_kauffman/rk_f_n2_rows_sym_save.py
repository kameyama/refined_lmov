"""rk_f_n2_rows_sym.build(kmax) の記号版ブロックを文字列にして保存する（サンドボックス外で実行）。
実行：PYTHONPATH=src:scripts/refined_kauffman python scripts/refined_kauffman/rk_f_n2_rows_sym_save.py 20"""
import sys, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
import rk_f_n2_rows_sym as RS

kmax = int(sys.argv[1]) if len(sys.argv) > 1 else 14
t0 = time.time()
N, x11, adams_coef, inv_den = RS.build(kmax, verbose=True)
out = {"N": {(lab, k): str(v.as_expr()) for lab in RS.LAMS for k, v in N[lab].c.items()},
       "x11": {k: str(v.as_expr()) for k, v in x11.c.items()},
       "adams_coef": {lab: str(v) for lab, v in adams_coef.items()},
       "inv_den": {k: str(v.as_expr()) for k, v in inv_den.c.items()}}
f = f"data/refined_kauffman/daha/n2_rows_sym_m1_k{kmax}.pkl"
pickle.dump(out, open(f, "wb"))
print(f"保存：{f}（{time.time() - t0:.0f} 秒、ブロック {sorted(out['N'])}）", flush=True)
