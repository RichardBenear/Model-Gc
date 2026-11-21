"""Claude review check for EX-GC26B: Lorentz boost with the original EX-GC26A NumPy solver (h=0.2, box 5).
Usage: python review/EX-GC26B_check/claude_boost_check.py 100 0.3   (or 111). Needs the gitignored 26A state."""
import os,sys,json
os.environ['OMP_NUM_THREADS']='1';os.environ['OPENBLAS_NUM_THREADS']='1'
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','src'))
import numpy as np
from scipy.ndimage import map_coordinates
import exgc26a_skyrmion as S
axis=[int(c) for c in sys.argv[1]];u=float(sys.argv[2]);T=0.8
G=S.Grid(0.2,5);p0=np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','..','results','EX-GC26A','relaxed_h0.2_L5.npy'));E0=G.energy(p0)
n=np.array(axis,float);n/=np.linalg.norm(n);gam=1/np.sqrt(1-u*u);X=np.array(G.X)
def boosted(t):
    xl=np.einsum('i,ixyz->xyz',n,X);xr=X+n[:,None,None,None]*((gam-1)*xl-gam*u*t)
    c=xr/G.h+(G.n-1)/2
    q=np.array([map_coordinates(p0[k],c,order=3,mode='nearest') for k in range(4)])
    return q/np.sqrt((q*q).sum(0))
eps=1e-3
p=G.fix(boosted(0.0));v=S.tangent(p,(boosted(eps)-boosted(-eps))/(2*eps))*G.inner
def dens(phi,v):
    H=0.
    for s in (1,-1):
        dp,D=G.strain(phi,s);tr=D[0,0]+D[1,1]+D[2,2]
        K=(v*v).sum(0)*(1+tr)-sum((v*dp[i]).sum(0)**2 for i in range(3));e=tr+.5*(tr*tr-(D*D).sum((0,1)))
        H=H+.5*(K+e)*G.mask[s]
    return H
dt=0.01;rows=[]
for i in range(81):
    if i%10==0:
        Hd=dens(p,v);H,P=S.hamiltonian(G,p,v);rows.append(dict(t=i*dt,H=H,P=P,C=[float((Hd*x).sum()/Hd.sum()) for x in X]))
    if i==80:break
    k1v=S.accel(G,p,v);k1p=v
    k2v=S.accel(G,G.fix(p+.5*dt*k1p),v+.5*dt*k1v);k2p=v+.5*dt*k1v
    k3v=S.accel(G,G.fix(p+.5*dt*k2p),v+.5*dt*k2v);k3p=v+.5*dt*k2v
    k4v=S.accel(G,G.fix(p+dt*k3p),v+dt*k3v);k4p=v+dt*k3v
    p=G.fix(p+dt/6*(k1p+2*k2p+2*k3p+k4p));v=S.tangent(p,v+dt/6*(k1v+2*k2v+2*k3v+k4v))*G.inner
t=np.array([r['t'] for r in rows]);C=np.array([r['C'] for r in rows])@n;P=np.array([r['P'] for r in rows]);H=np.array([r['H'] for r in rows])
vc=np.polyfit(t,C,1)[0];Pn=P@n;Ptr=np.linalg.norm(P-Pn[:,None]*n,axis=1).max()/E0
g=1/np.sqrt(1-vc*vc)
out=dict(axis=axis,u=u,E0=E0,v_centroid=vc,P_over_H=float((Pn/H).mean()),ratio=float(vc/(Pn/H).mean()),M_over_E0=float(Pn.mean()/(g*vc*E0)),P_over_gamma_u_E0=float(Pn.mean()/(gam*u*E0)),H_over_gamma_E0=float(H.mean()/(gam*E0)),dH=float(abs(H/H[0]-1).max()),transverse_P_max=float(Ptr))
print(json.dumps(out));json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),f'chk_{sys.argv[1]}_{u}.json'),'w'))
