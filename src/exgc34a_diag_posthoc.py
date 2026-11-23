import sys; sys.path.insert(0,'/home/user/model-gc/src')
import numpy as np, json, exgc34a as S, exgc34a_rays as R
from multiprocessing import Pool
# 1) wrap diagnosis: energy fraction in outer layers of final fields
def frac(tag, mk):
    z=np.load(f'/home/user/model-gc/local_raw/EX-GC34A/{tag}.npz'); f=z['f']
    e=(f**2) if f.ndim==3 else (f[1:]**2).sum(0) if f.shape[0]==4 else (f**2).sum(0)
    return float(e[mk(e.shape)].sum()/e.sum())
def edge_m1(sh): m=np.zeros(sh,bool); m[:10]=m[-10:]=True; return m   # within 2 units of x faces
def edge_m2(sh): m=np.zeros(sh,bool); m[:,-10:]=True; m[:,:10]=True; m[:,:,:10]=True; m[:,:,-10:]=True; return m
out={'M1_B0_x_F1_edge_frac':frac('M1_B0_x_F1',edge_m1),'M2_B2_F1_edge_frac':frac('M2_B2_F1',edge_m2)}
# 2) ray late-window angle for massless in bumps
def job(b):
    bg=R.Bg('bump',gppn={'B2':0.,'B3':1.,'B4':0.}[b],nordstrom=(b=='B4'),dt_lf=0.05)
    t=np.arange(0,36.0001,0.5); c=R.rays(bg,[-18,6,0],[3,0,0],2.0,0.0,0.2,36.0,0.05,nq=4,times=t)
    m=t>=28; v=[np.polyfit(t[m],c[m,a],1)[0] for a in range(2)]; return b,float(np.arctan2(v[1],v[0]))
with Pool(3) as P: out['ray_theta_late_massless']=dict(P.map(job,['B2','B3','B4']))
print(json.dumps(out,indent=1))
