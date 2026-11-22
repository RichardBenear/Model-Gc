"""EX-GC28A exploratory harmonic spectra; radial Skyrme is NOT full 3D."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[key]='1'
import json,time,resource,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import eigh_tridiagonal,solve_banded,eigh
from scipy.optimize import root
from scipy.integrate import solve_bvp

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC28A';START=time.process_time()
HS=(.2,.1,.05,.025);C0=.17839644285107545;ETA=.02
def save(tag,r):
    OUT.mkdir(parents=True,exist_ok=True);p=OUT/(tag+'.json');q=p.with_suffix('.tmp');q.write_text(json.dumps(r,indent=1,allow_nan=False)+'\n');q.replace(p)
def matrix(d,o):return np.diag(d)+np.diag(o,1)+np.diag(o,-1)
def spectrum(d,o,dv,ov,collective=None):
    vals,vec=eigh_tridiagonal(d,o);vv,V=eigh_tridiagonal(dv,ov);assert min(vv)>0
    rawvals=vals.copy();overlap=None
    if collective is not None:
        q,mu=collective;idx=int(np.argmax(abs(vec.T@q)));overlap=float(abs(vec[:,idx]@q));assert idx==0 and overlap>.99 and abs(vals[0])<.02*mu**2;vals[0]=0.
    assert min(vals)>=0,(vals[:4],overlap)
    w=np.sqrt(vals);wv=np.sqrt(vv);delta_d=d-dv;delta_o=o-ov
    trace_diag=np.sum(V**2*delta_d[:,None],axis=0)+2*np.sum(V[:-1]*V[1:]*delta_o[:,None],axis=0)
    raw=float(.5*np.sum(w-wv));born=float(.25*np.sum(trace_diag/wv));sub=raw-born
    # Independent vacuum-eigenbasis trace of Dpert, and frequency-bin convergence diagnostics.
    cuts=[float(wv[-1]*q) for q in (.25,.5,.75,1.)];partial=[]
    for cut in cuts:
        keep=wv<=cut;partial.append(dict(omega_cut=cut,modes=int(keep.sum()),raw=float(.5*np.sum((w-wv)[keep])),born=float(.25*np.sum((trace_diag/wv)[keep]))))
    return dict(raw=raw,born=born,sub=sub,lambda_min=float(rawvals[0]),collective_overlap=overlap,frequencies=w.tolist(),vacuum_frequencies=wv.tolist(),born_mode_terms=(.25*trace_diag/wv).tolist(),partial=partial)

def kink(h,R,eta=0.):
    mu=np.exp(eta);n=round(R/h);x=np.linspace(-R,R,2*n+1);h=x[1]-x[0];m=2*mu;xp=x[n:];f=np.tanh(mu*xp);f[0]=0;f[-1]=1
    def fun(y):
        z=np.r_[0.,y,1.];return (2*z[1:-1]-z[:-2]-z[2:])/h**2+2*mu**2*z[1:-1]*(z[1:-1]**2-1)
    def jac(y):return matrix(2/h**2+2*mu**2*(3*y*y-1),np.full(len(y)-1,-1/h**2))
    sol=root(fun,f[1:-1],jac=jac,tol=1e-10);f=np.r_[0.,sol.x,1.];p=np.r_[-f[:0:-1],f];res=float(max(abs(fun(sol.x))));assert res<1e-8
    E=float(np.sum(.5*np.diff(p)**2/h)+h*np.sum(.5*mu**2*(p[1:-1]**2-1)**2));d=2/h**2+2*mu**2*(3*p[1:-1]**2-1);o=np.full(len(d)-1,-1/h**2);dv=np.full_like(d,2/h**2+m*m);q=np.gradient(p,h)[1:-1];q/=np.linalg.norm(q)
    r=spectrum(d,o,dv,o,(q,mu));k=np.arange(1,len(d)+1);exact=np.sqrt(m*m+4/h**2*np.sin(np.pi*k/(2*(len(d)+1)))**2);r.update(h=h,R=R,eta=eta,mu=mu,E=E,stationary_residual=res,shape_omega2=r['frequencies'][1]**2,reference=mu*(1/(2*np.sqrt(3))-3/np.pi),vacuum_error=float(max(abs(np.array(r['vacuum_frequencies'])-exact))),cpu_total=time.process_time()-START);return r

def radial_parts(f,h,A,B,C):
    n=len(f)-1;r=(np.arange(n)+.5)*h;t=.5*(f[:-1]+f[1:]);z=np.diff(f)/h;s=np.sin(t);s2=np.sin(2*t);c2=np.cos(2*t);K=A*r*r+2*B*s*s;Kp=2*B*s2;Kpp=4*B*c2
    E2=h*np.sum(A*(r*r*z*z+2*s*s))/(3*np.pi);E4=h*np.sum(B*(2*s*s*z*z+s**4/r**2))/(3*np.pi);EV=h*np.sum(C*r*r*(1-np.cos(t)))/(3*np.pi)
    Lf=Kp*z*z+2*A*s2+2*B*s*s*s2/r**2+C*r*r*s;Ld=2*K*z;Lff=Kpp*z*z+4*A*c2+2*B*(s2*s2+2*s*s*c2)/r**2+C*r*r*np.cos(t);Lfd=2*Kp*z;Ldd=2*K
    F=h/(3*np.pi);gl=F*(.5*Lf-Ld/h);gr=F*(.5*Lf+Ld/h);g=gr[:-1]+gl[1:];dl=F*(.25*Lff-Lfd/h+Ldd/h**2);dr=F*(.25*Lff+Lfd/h+Ldd/h**2);d=dr[:-1]+dl[1:];o=F*(.25*Lff-Ldd/h**2)[1:-1];return np.array([E2,E4,EV]),g,d,o
def relax(h,R,A,B,C):
    n=round(R/h);h=R/n;r=np.arange(n+1)*h;f=2*np.arctan(np.divide(.85,r,out=np.full_like(r,1e8),where=r>0)**2);f[0]=np.pi;f[-1]=0
    for it in range(100):
        p,g,d,o=radial_parts(f,h,A,B,C);err=float(max(abs(g)));E=p.sum()
        if err<2e-11:break
        # Newton direction with diagonal shift only if Hessian not positive definite.
        lmin=float(eigh_tridiagonal(d,o,select='i',select_range=(0,0),eigvals_only=True)[0]);dd=d+max(0.,1e-8-lmin);ab=np.zeros((3,len(d)));ab[0,1:]=o;ab[1]=dd;ab[2,:-1]=o;step=solve_banded((1,1),ab,-g);slope=float(g@step);assert slope<0
        for j in range(35):
            q=f.copy();q[1:-1]+=2.**(-j)*step;new=radial_parts(q,h,A,B,C)[0].sum()
            if new<=E+1e-4*2.**(-j)*slope+1e-15:f=q;break
        else:raise RuntimeError('line search failed')
    p,g,d,o=radial_parts(f,h,A,B,C);assert max(abs(g))<2e-8,(h,R,A,B,C,max(abs(g)),it);return f,p,g,d,o,it
def radial(h,R,Cbase,eta=0.,kind='consistent'):
    A=np.exp(2*eta);B=1.;C=Cbase*np.exp(4*eta) if kind=='consistent' else Cbase;n=round(R/h);h=R/n;f,p,g,d,o,it=relax(h,R,A,B,C);r=np.arange(n+1)*h;M=2*h*(A*r[1:-1]**2+2*B*np.sin(f[1:-1])**2)/(3*np.pi);_,_,dv,ov=radial_parts(np.zeros_like(f),h,A,B,C);Mv=2*h*A*r[1:-1]**2/(3*np.pi);D=d/M;O=o/np.sqrt(M[:-1]*M[1:]);Dv=dv/Mv;Ov=ov/np.sqrt(Mv[:-1]*Mv[1:]);z=spectrum(D,O,Dv,Ov);E=float(p.sum());z.update(h=h,R=R,Cbase=Cbase,C=C,A=A,B=B,eta=eta,kind=kind,E=E,parts=p.tolist(),virial=float((p[0]-p[1]+3*p[2])/E),s_classical=float(2*p[0]+(4*p[2] if kind=='consistent' else 0)),stationary_residual=float(max(abs(g))),iterations=it,profile=f.tolist(),mass_metric_min=float(M.min()),cpu_total=time.process_time()-START);return z

def validate():
    h=.2;R=4;A=1.13;B=.91;C=.18;f,p,g,d,o,it=relax(h,R,A,B,C);rng=np.random.default_rng(28);q=rng.normal(size=len(g));q/=np.linalg.norm(q);eps=1e-5;pp=f.copy();pm=f.copy();pp[1:-1]+=eps*q;pm[1:-1]-=eps*q;ep,gp,*_=radial_parts(pp,h,A,B,C);em,gm,*_=radial_parts(pm,h,A,B,C);err=float(np.max(abs((gp-gm)/(2*eps)-matrix(d,o)@q)));gerr=float(abs((ep.sum()-em.sum())/(2*eps)-g@q));r=np.arange(1,len(f)-1)*h;M=2*h*(A*r*r+2*B*np.sin(f[1:-1])**2)/(3*np.pi);v=eigh(matrix(d,o),np.diag(M),eigvals_only=True);w=eigh_tridiagonal(d/M,o/np.sqrt(M[:-1]*M[1:]),eigvals_only=True);dense=float(max(abs(v-w)));save('validation',dict(hessian_error=err,gradient_error=gerr,generalized_eigenvalue_error=dense,min_lambda=float(min(w)),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__));assert err<1e-7 and gerr<1e-7 and dense<1e-9
def main():
    resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));validate();k=[];s=[]
    for R in (12,20):
        for h in HS:
            row=kink(h,R);k.append(row);save('kink',k);print('kink',h,R,row['raw'],row['sub'],flush=True)
    for h in (.1,.05):
        for eta in (-ETA,ETA):k.append(kink(h,20,eta));save('kink',k)
    for C in (0.,C0):
        for R in (12,24):
            for h in HS:
                row=radial(h,R,C);s.append(row);save('radial',s);print('radial',C,h,R,row['E'],row['raw'],row['sub'],flush=True)
        for h in (.1,.05):
            for kind in ('consistent','independent'):
                for eta in (-ETA,ETA):s.append(radial(h,24,C,eta,kind));save('radial',s)
            for eta in (-ETA,ETA):
                lam=np.exp(-eta);row=radial(h*lam,24*lam,C,eta);row['comoving']=True;row['nominal_h']=h;s.append(row);save('radial',s)
    save('completion',dict(cpu_s=time.process_time()-START,complete=True));print('DONE',time.process_time()-START,flush=True)
def canonical(h,R,Cbase,eta=0.,kind='consistent'):
    A=np.exp(2*eta);B=1.;C=Cbase*np.exp(4*eta) if kind=='consistent' else Cbase
    def terms(r,f,z):
        s=np.sin(f);s2=np.sin(2*f);c2=np.cos(2*f);K=A*r*r+2*B*s*s;Kf=2*B*s2;Kff=4*B*c2
        Uf=2*A*s2+2*B*s*s*s2/r**2+C*r*r*s
        Uff=4*A*c2+2*B*(s2*s2+2*s*s*c2)/r**2+C*r*r*np.cos(f)
        zz=(Uf-Kf*z*z-4*A*r*z)/(2*K)
        return K,Kf,Kff,Uff,zz
    def rhs(r,y):return np.vstack((y[1],terms(r,*y)[-1]))
    r0=1e-4;rr=np.linspace(r0,R,1000);f,p,*_=relax(.025,R,A,B,C);rg=np.linspace(0,R,len(f));guess=np.interp(rr,rg,f)
    sol=solve_bvp(rhs,lambda ya,yb:np.array([ya[0]-np.pi-r0*ya[1],yb[0]]),rr,np.vstack((guess,np.gradient(guess,rr))),tol=1e-7,max_nodes=30000)
    assert sol.success,sol.message
    n=round(R/h);h=R/n;r=np.arange(1,n)*h;f,z=sol.sol(r);K,Kf,Kff,Uff,zz=terms(r,f,z);K1=2*A*r+Kf*z;K2=2*A+Kff*z*z+Kf*zz
    eff=.5*K2/K-.25*(K1/K)**2+(Uff-Kff*z*z-2*Kf*zz)/(2*K);vac=2/r**2+C/(2*A);o=np.full(n-2,-1/h**2)
    out=spectrum(2/h**2+eff,o,2/h**2+vac,o);out.update(h=h,R=R,Cbase=Cbase,C=C,A=A,B=B,eta=eta,kind=kind,profile=f.tolist(),delta_potential=(eff-vac).tolist(),bvp_residual=float(max(sol.rms_residuals)),cpu_total=time.process_time()-START);return out

def canonical_run():
    resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));rows=[]
    for C in (0.,C0):
        for R in (12,24):
            for h in HS:
                row=canonical(h,R,C);rows.append(row);save('canonical',rows);print(C,h,R,row['raw'],row['sub'],flush=True)
        for h in (.1,.05):
            for kind in ('consistent','independent'):
                for eta in (-ETA,ETA):rows.append(canonical(h,24,C,eta,kind));save('canonical',rows)
            for eta in (-ETA,ETA):
                lam=np.exp(-eta);row=canonical(h*lam,24*lam,C,eta);row.update(comoving=True,nominal_h=h);rows.append(row);save('canonical',rows)
    save('canonical_completion',dict(cpu_s=time.process_time()-START,complete=True))

if __name__=='__main__':
    import sys
    canonical_run() if 'canonical' in sys.argv else main()
