"""EX-GC34E POST HOC (not registered): extend the re-slicing-direction scan from w = 3k to wh = 2.5 (lattice cutoff).
-> results/EX-GC34E/posthoc_cutoff.json"""
import sys; from pathlib import Path; ROOT=Path(__file__).resolve().parents[1]; sys.path.insert(0,str(ROOT/'src')); import json; OUTL=[]
import numpy as np, scipy.linalg as sl, mpmath as mp
import exgc34e as E
mp.mp.dps=50; st=E.tet_static()
def sv(q,w):
    M=E.full_kernel(st,q[None],w)[0]; G=E.gauge_modes(q[None],w)[0]; Q=sl.null_space(G.conj().T)
    s=np.linalg.svd(Q.conj().T@M@Q,compute_uv=False); return s[::-1]/s[0]
for kh,dv in ((0.1,(0.8,0.5,0.33)),(0.05,(1,1,1)),(0.2,(1,1,0))):
    q=kh*np.array(dv)/np.linalg.norm(dv); k=kh/E.H
    ws=np.geomspace(3,2.5/E.H/k,600)*k
    s=np.array([sv(q,w)[:2] for w in ws])
    i=np.argmin(s[:,0]); OUTL.append(dict(kh=kh,dir=list(dv),s0_min=float(s[i,0]),wh_at_min=float(ws[i]*E.H),s0_lo=float(s[0,0]),s0_hi=float(s[-1,0]),s1_min=float(s[:,1].min()))); print(kh,dv,'min s0 %.2e at wh=%.3f'%(s[i,0],ws[i]*E.H), 's0 at ends %.2e %.2e'%(s[0,0],s[-1,0]), 'min s1 %.2e at wh=%.3f'%(s[:,1].min(), ws[np.argmin(s[:,1])]*E.H))
(ROOT/'results/EX-GC34E/posthoc_cutoff.json').write_text(json.dumps(OUTL,indent=1)+'\n')
