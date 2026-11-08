r"""rk_c_3box_largeN の組み立て（数分）をやり直し、g だけを sympy の文字列で data/refined_kauffman/daha/largeN_3box_m{m}_g.pkl に保存する
（f, f^cov を含む全体の保存は m = 4 で 6 時間を超えたため。progress.md 続き 139）。形式は largeN_3box_m{m}.pkl の "g" と同じ。
使い方：python rk_c_3box_g_only.py m"""
import sys, os, pickle, time
sys.path.insert(0, "src"); sys.path.insert(0, "scripts/refined_kauffman")
from fractions import Fraction
import rk_c_3box_largeN as T3

m = int(sys.argv[1])
t0 = time.time()
T3.MK = m
F = T3.plog(T3.so_Z(T3.so_W(), T3.K_matrix()))
Wo, pair = T3.cover_parts()
Fc = T3.plog(T3.cover_Z(Wo, pair))
G = {k: F.get(k, T3.ZERO) - T3.conv(Fraction(1, 2)) * Fc.get(k, T3.ZERO) for k in set(F) | set(Fc)}
print(f"組み立て {time.time()-t0:.0f}s", flush=True)
pickle.dump({"g": {str(k): str(v.to_sympy()) for k, v in G.items()}}, open(os.path.join(T3.OUT, f"largeN_3box_m{m}_g.pkl"), "wb"))
print(f"保存 {time.time()-t0:.0f}s", flush=True)
