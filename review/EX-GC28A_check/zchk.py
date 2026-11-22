# Independent recomputation of the EX-GC28A canonical radial Z (separate code path; shares numpy/scipy only).
# Usage: python zchk.py 0   and   python zchk.py 0.17839644285107545   (about 1-2 min each). Output in zchk_output.txt.
# Profile: own energy minimisation of the hedgehog on a fine grid (L-BFGS); effective potential from
# numerical second variation (finite differences of the energy density), not the analytic formula.
import numpy as np, sys
from scipy.optimize import minimize
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh_tridiagonal
C=float(sys.argv[1]); R=24.0; A=1.0; B=1.0
n=4800; r=np.linspace(0,R,n+1); dr=r[1]-r[0]
def E(fi):
    f=np.r_[np.pi,fi,0.]; t=.5*(f[1:]+f[:-1]); z=np.diff(f)/dr; rm=.5*(r[1:]+r[:-1]); s=np.sin(t)
    return dr*np.sum((A*rm**2+2*B*s*s)*z*z+2*A*s*s+B*s**4/rm**2+C*rm**2*(1-np.cos(t)))/(3*np.pi)
def G(fi):
    f=np.r_[np.pi,fi,0.]; t=.5*(f[1:]+f[:-1]); z=np.diff(f)/dr; rm=.5*(r[1:]+r[:-1]); s=np.sin(t); c=np.cos(t)
    dt=(2*B*2*s*c)*z*z+2*A*2*s*c+B*4*s**3*c/rm**2+C*rm**2*s; dz=2*(A*rm**2+2*B*s*s)*z
    g=.5*dt*dr; gz=dz; gl=g-gz; gr=g+gz  # derivative wrt left node f_i: .5*dt*dr - dz ; right node f_{i+1}: .5*dt*dr + dz
    tot=np.zeros(n+1); tot[:-1]+=gl; tot[1:]+=gr; return tot[1:-1]/(3*np.pi)
f0=2*np.arctan(np.divide(1.0,r,out=np.full_like(r,1e8),where=r>0)**2)  # crude start
res=minimize(E,f0[1:-1],jac=G,method='L-BFGS-B',options=dict(maxcor=50,maxiter=200000,maxfun=400000,gtol=1e-13,ftol=1e-16))
f=np.r_[np.pi,res.x,0.]; print('classical E',E(res.x),'grad',np.abs(G(res.x)).max())
fs=CubicSpline(r,f); 
# second variation numerically: density e(r,f,fp) = K(f,r) fp^2 + U(f,r); canonical potential
def K(r,f): return A*r*r+2*B*np.sin(f)**2
def U(r,f): return 2*A*np.sin(f)**2+B*np.sin(f)**4/r**2+C*r*r*(1-np.cos(f))
def V(rr):
    fr=fs(rr); fp=fs(rr,1); e=2e-3
    Uff=(U(rr,fr+e)-2*U(rr,fr)+U(rr,fr-e))/e**2; Kf=(K(rr,fr+e)-K(rr,fr-e))/(2*e); Kff=(K(rr,fr+e)-2*K(rr,fr)+K(rr,fr-e))/e**2
    # second variation of density: K dfp^2 + [Kf fp] d(df^2)/dr... -> potential term q = (Uff + Kff fp^2)/2 - d/dr(Kf fp)
    q=lambda x:None
    w=.5*(Uff+Kff*fp**2)
    d=2e-3; dKf=lambda x:((K(x,fs(x)+e)-K(x,fs(x)-e))/(2*e))*fs(x,1)
    w-= (dKf(rr+d)-dKf(rr-d))/(2*d)
    sq=lambda x:np.sqrt(K(x,fs(x)))
    sqpp=(sq(rr+d)-2*sq(rr)+sq(rr-d))/d**2
    return w/K(rr,fr)+sqpp/sq(rr)
for h in (.1,.05,.025):
    m=round(R/h); x=np.arange(1,m)*h; Vd=V(x); Vv=2/x**2+C/(2*A)
    o=np.full(m-2,-1/h**2); w=np.sqrt(eigh_tridiagonal(2/h**2+Vd,o,eigvals_only=True)); vv,Q=eigh_tridiagonal(2/h**2+Vv,o)
    wv=np.sqrt(vv); born=.25*np.sum((Q**2*(Vd-Vv)[:,None]).sum(0)/wv)
    print(h,'Z',.5*np.sum(w-wv)-born)
