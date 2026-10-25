"""framed knot の refined LMOV：ユーザーの以前の計算（Dropbox/Mathematica/marino-vafa のノートブック）の出力と、本研究のコード（scripts/framed_knot/framing_check.py）の結果を全係数で比べる。

ノートブックは scripts/framed_knot/nb2txt.py でテキストにしたもの（既定は $TMPDIR/mv/。MV_DIR で変更）。結果は logs/framed_knot/framing_notebook_check.log。
使い方：MV_DIR=<テキストのフォルダ> python scripts/framed_knot/framing_notebook_check.py
"""
import sys, re; sys.path.insert(0, "scripts/framed_knot")
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
from framing_check import run
from rcsknots.bps import to_bold
from rcsknots.extract import A_SHIFT
A,T,Q,g=sp.symbols("A T Q g")
TR=standard_transformations+(implicit_multiplication_application,)
def mparse(s): return sp.expand(parse_expr(s.replace("^","**").replace("\\IndentingNewLine",""),local_dict={"A":A,"T":T,"Q":Q,"g":g},transformations=TR))
def outputs(path):
    txt=open(path).read().split("=== [")
    pairs=[]; last=None
    for blk in txt:
        if blk.startswith("Input]"): last=blk[len("Input]"):].strip()
        elif blk.startswith("Output]") and last: pairs.append((last, blk[len("Output]"):].strip())); last=None
    return pairs
import os
D=os.environ.get("MV_DIR", os.path.join(os.environ.get("TMPDIR", "/tmp"), "mv")) + "/"
fr=outputs(D+"framed_results.txt"); v2=outputs(D+"framed refined Oogri-Vafa ver2.txt")
def find(pairs,key):
    hits=[o for i,o in pairs if key in i]; return hits
cases=[ # (結び目, x, rho, ノートブックの入力の目印, ファイル, 種類, 符号)
 ("T2_3",2,(1,),"hatf[{1},2]/a^(1/2)",fr,"G",1),
 ("T2_3",2,(2,),"-hatf[{2},2]",fr,"G",-1),
 ("T2_3",2,(1,1),"-hatf[{1,1},2]",fr,"G",-1),
 ("unknot",3,(3,),"test3//knotparam//gfactor",fr,"G",1),
 ("unknot",3,(2,1),"hatf21*a^(1/2)/.x->3",fr,"G",1),
 ("unknot",3,(1,1,1),"test111=a^(1/2)*hatf111",fr,"B",1),
 ("T2_3",2,(3,),"hatf3*a^(1/2)/.x->2",v2,"G",1),
 ("4_1",0,(2,),"-hatf2/.x->0",v2,"G",-1),
 ("4_1",2,(2,),"-hatf2/.x->2",v2,"G",-1),
]
cache={}
for name,x,rho,key,pairs,kind,sg in cases:
    outs=find(pairs,key)
    if not outs: print(name,x,rho,"ノートブックに出力なし"); continue
    user=mparse(outs[0])
    k=(name,x,sum(rho))
    if k not in cache: cache[k]=run(name,[x],nmax=sum(rho),log=lambda s:None)
    msg,N,hat=cache[k][(x,rho)]
    if kind=="G":
        P=sum(c*A**(2*b)*g**gg*T**F for (gg,b,F),c in N.items())
    else:
        bold=to_bold(hat,A_SHIFT["nontorus"][sum(rho)])
        P=sum(c*A**ea*Q**eq*T**et for (ea,eq,et),c in bold.items())
    ok=sp.expand(P-user)==0; okm=sp.expand(P+user)==0
    print(f"{name} framing {x} ρ={list(rho)} [{kind}]：本研究 {msg}／ノートブックと {'一致' if ok else ('符号違いで一致' if okm else '不一致')}",flush=True)
