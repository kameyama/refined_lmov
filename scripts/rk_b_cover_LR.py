r"""拡張版の S2：被覆側 Σ N N f^{cov}_{R,S}（|R|, |S| ≤ 4、Y' = Y、Z' = Z）を求めて保存する（kauffman_hopf2.cover_f_LR）。
出力：data/refined_kauffman/hopf_fcovLR_A4.pkl。完了の目印：ログの最後の行「保存 完了」。"""
import sys, time, pickle
sys.path.insert(0, "src")
from rcsknots.kauffman_hopf2 import cover_f_LR

if __name__ == "__main__":
    t0 = time.time()
    f = cover_f_LR(4, 4)
    print(f"cover_f_LR(4,4): {time.time()-t0:.0f}s、鍵 {len(f)}", flush=True)
    pickle.dump({k: str(v) for k, v in f.items()}, open("data/refined_kauffman/hopf_fcovLR_A4.pkl", "wb"))
    print("保存 完了", flush=True)
