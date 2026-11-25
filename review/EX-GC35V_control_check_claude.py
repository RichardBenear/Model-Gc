"""Claude check of the GC35V published overlapping-rule control under four graph conventions (vertex labels x input anchors). Run: python review/EX-GC35V_control_check_claude.py"""
import sys,time,pathlib; sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"src"))
import exgc35v as V
from itertools import combinations
orig=V.graph
def unl(state,seed,types,anchor=True):
    g=orig(state,seed,types,anchor)
    for n,d in g.nodes(data=True):
        if n[0]=="v": d["label"]="vtx"
    return g
rule=((("C",0,1),("C",2,1)),(("C",0,3),("C",1,3),("C",2,3)))
seed=(("C",0,3),("C",1,3),("C",2,3))
st=V.make(seed)
for mode in ["labelled","unlabelled"]:
    V.graph = orig if mode=="labelled" else unl
    for anchor in [True,False]:
        t=time.time()
        res=[V.join(rule,st,a,b,depth=3,anchor=anchor) for a,b in combinations(V.choices(rule,st[0]),2) if set(a[0])&set(b[0])]
        print(mode,anchor,res,round(time.time()-t,1),flush=True)
