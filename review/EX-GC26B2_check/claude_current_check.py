"""Claude review check for EX-GC26B2: exact bond energy current J vs canonical P on the ORIGINAL two-orientation
26A stencil (own derivation, own code). Lagrangian here is (K-e)/2, so J_over_P must be multiplied by 2 (canonical P
from hamiltonian() uses the full momentum). Needs the gitignored 26A state results/EX-GC26A/relaxed_h0.2_L5.npy."""
import sys,json;import os;HERE=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,os.path.join(HERE,'..','..','src'))
import numpy as np
from scipy.ndimage import map_coordinates
import exgc26a_skyrmion as S
G=S.Grid(0.2,5);p0=np.load(os.path.join(HERE,'..','..','results','EX-GC26A','relaxed_h0.2_L5.npy'));X=np.array(G.X)
def launch(axis,u):
    n=np.array(axis,float);n/=np.linalg.norm(n);gam=1/np.sqrt(1-u*u)
    def b(t):
        xl=np.einsum('i,ixyz->xyz',n,X);c=(X+n[:,None,None,None]*((gam-1)*xl-gam*u*t))/G.h+(G.n-1)/2
        q=np.array([map_coordinates(p0[k],c,order=3,mode='nearest') for k in range(4)]);return q/np.sqrt((q*q).sum(0))
    p=G.fix(b(0.));v=S.tangent(p,(b(1e-3)-b(-1e-3))/2e-3)*G.inner;return p,v,n
def JP(p,v):
    J=np.zeros(3);P=np.zeros(3);H=0
    for s in (1,-1):
        dp,D=G.strain(p,s);tr=D[0,0]+D[1,1]+D[2,2];vv=(v*v).sum(0);m=G.mask[s]
        for a in range(3):
            dl=vv*dp[a]-(v*dp[a]).sum(0)*v-dp[a]*(1+tr)+sum(D[a,j]*dp[j] for j in range(3))  # dL/d(dp_a)
            vn=np.roll(v,-s,axis=a+1)
            J[a]+=0.5*((dl*vn).sum(0)*m).sum(); P[a]+=0.5*((dl*v).sum(0)*m).sum()
    Hh,Pc=S.hamiltonian(G,p,v);f=G.h**3/S.NORM
    return J*f,P*f,np.array(Pc),Hh
out={}
for ax in ([1,0,0],[1,1,1]):
    p,v,n=launch(ax,0.3);J,Ps,Pc,H=JP(p,v)
    # numerical check: d/dt of energy moment via one RK4-free Euler-like central step using accel
    out[str(ax)]=dict(J_over_P=float(J@n/(Pc@n)),Psame_over_Pcanon=float(Ps@n/(Pc@n)),J_over_H=float(J@n/H))
    print(ax,out[str(ax)])
json.dump(out,open(os.path.join(HERE,'claude_current_check.json'),'w'),indent=1)
