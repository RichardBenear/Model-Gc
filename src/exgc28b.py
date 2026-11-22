"""Connected vacuum interaction of disjoint, positive scalar scatterers (not gravity)."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import json, time, resource, platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import cholesky, solve_triangular, eigvalsh, eigh_tridiagonal
from scipy.spatial.distance import cdist
from scipy.special import roots_legendre
from scipy.integrate import quad

OUT=Path(__file__).resolve().parents[1]/'results/EX-GC28B'
START=time.process_time()

def shape(dim,h,sigma,Q):
    z=np.arange(-int(1/h),int(1/h)+1)*h
    x=np.stack(np.meshgrid(*([z]*dim),indexing='ij'),axis=-1).reshape(-1,dim)
    r=np.linalg.norm(x,axis=1); x=x[r<1-1e-12]; r=r[r<1-1e-12]
    v=np.exp(-r*r/(2*sigma*sigma)); v*=Q/v.sum()
    return x,v

def greens(r,k,dim,h):
    if dim==1:return np.exp(-k*r)/(2*k)
    # Equal-volume sphere average of Yukawa kernel at self point.
    a=h*(3/(4*np.pi))**(1/3); u=k*a
    diag=3/(8*np.pi*a) if u<1e-5 else 3*(-np.expm1(-u)-u*np.exp(-u))/(4*np.pi*a*u*u)
    return np.where(r>0,np.exp(-k*r)/(4*np.pi*np.where(r>0,r,1)),diag)

def setup(dim,h,Q1,Q2,s1,s2,d):
    x,v=shape(dim,h,s1,Q1); y,w=shape(dim,h,s2,Q2); shift=np.zeros(dim);shift[0]=d
    return cdist(x,x),cdist(y,y),cdist(x,y+shift),np.sqrt(v),np.sqrt(w)

def logs(r1,r2,r12,v,w,k,dim,h):
    A=np.eye(len(v))+v[:,None]*greens(r1,k,dim,h)*v
    C=np.eye(len(w))+w[:,None]*greens(r2,k,dim,h)*w
    X=v[:,None]*greens(r12,k,dim,h)*w
    L=cholesky(A,lower=True); M=cholesky(C,lower=True)
    B=solve_triangular(M,solve_triangular(L,X,lower=True).T,lower=True).T
    # Stable far-tail: logdet(I-BB^T) via eigenvalues or trace expansion.
    P=B@B.T; t=np.trace(P)
    if t<1e-4: return -t-.5*np.sum(P*P)-np.trace(P@P@P)/3
    ev=eigvalsh(P);assert ev[-1]<1,(k,ev[-1]);return np.log1p(-np.maximum(ev,0)).sum()

def energy(dim,h,d,Q1,Q2,m=0,s1=.35,s2=.35,n=48,end=32,square=False):
    r1,r2,r12,v,w=setup(dim,h,Q1,Q2,s1,s2,d)
    nodes,weights=roots_legendre(n); t=(nodes+1)*end/2; weights*=end/2
    if square:
        u=t/np.sqrt(end);weights*=2*u/np.sqrt(end);t=u*u
    l=np.array([logs(r1,r2,r12,v,w,np.hypot(u/d,m),dim,h) for u in t]);E=float(weights@l/(2*np.pi*d))
    result=dict(dim=dim,h=h,d=d,Q1=Q1,Q2=Q2,m=m,s1=s1,s2=s2,n=n,end=end,energy=E,n1=len(v),n2=len(w),square=square)
    if dim==3:
        def a(r,v):return float(v@np.linalg.solve(np.eye(len(v))+v[:,None]*greens(r,0,3,h)*v,v)/(4*np.pi))
        a1,a2=a(r1,v),a(r2,w);result.update(a1=a1,a2=a2,asymptotic=-a1*a2/(4*np.pi*d**3),born_point=-Q1*Q2/(64*np.pi**3*d**3))
    return result

def validate():
    # Single point potentials: exact reflection formula vs whitened determinant.
    d=7.; q1=.8;q2=1.7;m=.15
    def exact(x):
        k=np.hypot(x,m);return np.log1p(-q1*q2*np.exp(-2*k*d)/((2*k+q1)*(2*k+q2)))
    ref=quad(exact,0,np.inf,epsabs=1e-13)[0]/(2*np.pi)
    test=quad(lambda x:logs(np.zeros((1,1)),np.zeros((1,1)),np.full((1,1),d),np.array([q1**.5]),np.array([q2**.5]),np.hypot(x,m),1,.1),0,np.inf,epsabs=1e-13)[0]/(2*np.pi)
    assert abs(ref-test)<1e-11
    # Independent finite box lattice spectrum vs exact same finite-box matrix determinant integral.
    h=.2;x=np.arange(-10+h,10,h);N=len(x);mass=.3;dv=np.full(N,2/h**2+mass**2);o=np.full(N-1,-1/h**2)
    v1=np.where(abs(x+2)<1,np.exp(-(x+2)**2/.25),0);v2=np.where(abs(x-2)<1,.7*np.exp(-(x-2)**2/.25),0)
    def w(v):return np.sqrt(eigh_tridiagonal(dv+v,o,eigvals_only=True))
    direct=.5*np.sum(w(v1+v2)-w(v1)-w(v2)+w(0))
    D=np.diag(dv)+np.diag(o,1)+np.diag(o,-1); lam,U=np.linalg.eigh(D);i=np.flatnonzero(v1);j=np.flatnonzero(v2)
    def integ(z):
        G=(U/(lam+z*z))@U.T;v=np.sqrt(v1[i]);u=np.sqrt(v2[j]);L=cholesky(np.eye(len(i))+v[:,None]*G[np.ix_(i,i)]*v,lower=True);M=cholesky(np.eye(len(j))+u[:,None]*G[np.ix_(j,j)]*u,lower=True);B=solve_triangular(M,solve_triangular(L,v[:,None]*G[np.ix_(i,j)]*u,lower=True).T,lower=True);return np.log1p(-eigvalsh(B@B.T)).sum()
    det=quad(integ,0,np.inf,epsabs=1e-12)[0]/(2*np.pi);assert abs(direct-det)<1e-10
    null=energy(3,.4,4,0,1)['energy'];assert null==0
    return dict(delta_reference=ref,delta_matrix=test,lattice_spectral=float(direct),lattice_det=det,null=null)

def main():
    resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));OUT.mkdir(parents=True,exist_ok=True); checks=validate();print(checks,flush=True);rows=[]
    def run(dim,h,d,q,p,m=0,s1=.35,s2=.35,n=48,end=32):
        r=energy(dim,h,d,q,p,m,s1,s2,n,end);rows.append(r);print(dim,h,d,q,p,m,r['energy'],flush=True)
        (OUT/'records.json').write_text(json.dumps(dict(validation=checks,rows=rows,complete=False),indent=1)+'\n')
    for dim,h in [(1,.05),(3,.25)]:
        for q,p,s1,s2 in [(.05,.05,.35,.35),(1,1,.35,.35),(10,10,.35,.35),(1,10,.35,.35),(1,1,.35,.5)]:
            for m in (0,.15):
                for d in (4,8,16,32,64):run(dim,h,d,q,p,m,s1,s2)
        for hh in ([.1] if dim==1 else [.4,.2]):
            for q in (.05,1,10):
                for d in (4,16,64):run(dim,hh,d,q,q)
        for d in (4,64):
            for q in (.05,10):run(dim,h,d,q,q,n=96,end=48)
    (OUT/'records.json').write_text(json.dumps(dict(validation=checks,rows=rows,complete=True,cpu_s=time.process_time()-START,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),indent=1)+'\n');print('DONE',time.process_time()-START,flush=True)

def infrared():
    rows=[]
    def run(h,d,q,p,s1=.35,s2=.35,n=96):
        r=energy(1,h,d,q,p,0,s1,s2,n,48,True);rows.append(r);print(h,d,q,p,n,r['energy'],flush=True)
    for q,p,s1,s2 in [(.05,.05,.35,.35),(1,1,.35,.35),(10,10,.35,.35),(1,10,.35,.35),(1,1,.35,.5)]:
        for d in (4,8,16,32,64):run(.05,d,q,p,s1,s2)
    for q in (.05,1,10):
        for d in (4,16,64):run(.1,d,q,q)
        for d in (128,256,512,1024):run(.05,d,q,q)
        for d in (4,64,1024):run(.05,d,q,q,n=192)
    (OUT/'infrared.json').write_text(json.dumps(dict(rows=rows,cpu_s=time.process_time()-START,complete=True),indent=1)+'\n')

if __name__=='__main__':
    import sys
    infrared() if 'infrared' in sys.argv else main()
