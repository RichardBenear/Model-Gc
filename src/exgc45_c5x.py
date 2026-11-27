import sys, json; sys.path.insert(0, "src")
import numpy as np, exgc45 as X
C = json.loads(open("results/EX-GC45/C.json").read())["result"]
h = 1e-5; out = {}
for key, c in C.items():
    m0, g = c["m0"], c["g"]
    for Lc in (6,):
        pr = [X.pair_dense_full(Lc, 0.1, m0, g, c["U"], (3,), t=(1+s*h,)*3) for s in (1,-1)]
        fu = [X.pair_dense_full(Lc, 0.1, m0, g, c["U"], (3,4,5,6), t=(1+s*h,)*3) for s in (1,-1)]
        neg = lambda a: np.roll(np.flip(a, axis=(0,1,2)), 1, axis=(0,1,2))
        ed = lambda t: float((lambda e: (e[...,3]+neg(e[...,3])).min())(X.eig(Lc, 0.1, m0=m0, g=g, t=(t,)*3)[0]))
        Qe = (ed(1+h)-ed(1-h))/(2*h)
        dp = (pr[0]-pr[1])/(2*h)-Qe; df = (fu[0]-fu[1])/(2*h)-Qe
        Ep = X.pair_dense_full(Lc, 0.1, m0, g, c["U"], (3,)); Ef = X.pair_dense_full(Lc, 0.1, m0, g, c["U"], (3,4,5,6))
        out[key] = dict(L=Lc, dQa_proj=dp, dQa_4band=df, diff=df-dp, rel=abs(df-dp)/abs(dp), diff_over_Eb=(df-dp)/c["Eb"], E_proj=Ep, E_4band=Ef)
        print(key, out[key], flush=True)
json.dump(out, open("results/EX-GC45/C5_extended.json","w"), indent=1)
