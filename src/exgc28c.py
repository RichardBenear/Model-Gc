"""One-loop induced auxiliary-coordinate stiffness and monopole tests, EXPLORATORY."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,time,resource,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import exp1,roots_legendre,erf
from scipy.integrate import quad
from scipy.optimize import brentq

OUT=Path(__file__).resolve().parents[1]/'results/EX-GC28C'; C=1/(16*np.pi**2)
x,w=roots_legendre(64);x=(x+1)/2;w=w/2;y=x*(1-x)
t,wt=roots_legendre(24);t=(t+1)/2;wt=wt/2
def I(a,L):return C*(L*L*np.exp(-a/L**2)-a*exp1(a/L**2))
def B(a,b,L,z=0):return C*np.dot(w,exp1((x*a+(1-x)*b+y*z)/L**2))
def diff(a,b,L,z):
    # B(0)-B(z) without cancellation, by integrating E1 derivative.
    d=x*a+(1-x)*b; v=d[:,None]+y[:,None]*z*t
    return C*z*np.dot(w*y,(np.exp(-v/L**2)/v)@wt)
def slope(a,b,L,z=0):
    v=x*a+(1-x)*b+y*z;return C*np.dot(w*y,np.exp(-v/L**2)/v)

def model(A,L,a=1,b=4):
    if A=='scalar':b=a;f=1.
    else:f=(a-b)**2
    B0=B(a,b,L);Z=f*slope(a,b,L);return a,b,f,B0,Z

def response(a,b,L,f,Z,mu,r,source='unit'):
    delta=Z*mu*mu
    def K(p):return delta+f*diff(a,b,L,p*p)
    def J(p):return np.exp(-p*p/2)*(1 if source=='unit' else B(a,b,L,p*p)/B(a,b,L))
    if mu==0:
        # Subtract exactly the Gaussian-smeared Coulomb pole and integrate smooth remainder.
        pole=erf(r/np.sqrt(2))/(4*np.pi*Z*r)
        def integrand(p):
            if p<1e-7:return 0.
            return p*(J(p)/K(p)-np.exp(-p*p/2)/(Z*p*p))
    else:
        pole=0.
        def integrand(p):return p*J(p)/K(p)
    value,error=quad(integrand,0,12,weight='sin',wvar=r,epsabs=1e-10,epsrel=1e-10,limit=200)
    return float(pole+value/(2*np.pi**2*r)),float(error/(2*np.pi**2*r))

def validate():
    rows=[]
    for L in (2,5,10,20):
        for a,b in ((1,1),(1,4)):
            ref=(I(a,L)-I(b,L))/(b-a) if a!=b else -(I(a+1e-4,L)-I(a-1e-4,L))/(2e-4)
            be=abs(B(a,b,L)-ref);fd=(B(a,b,L,-1e-4)-B(a,b,L,1e-4))/(2e-4);se=abs(fd-slope(a,b,L));assert be<1e-9 and se<1e-10
            # Constant rotation preserves spectrum; tadpole/seagull and bubble cancel.
            d=a-b;theta=.234;R=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]]);e=np.linalg.eigvalsh(R@np.diag([a,b])@R.T);rot=float(max(abs(e-[a,b])));seagull=d*(I(b,L)-I(a,L));bubble=d*d*B(a,b,L);assert abs(seagull-bubble)<1e-12 and rot<1e-12
            rows.append(dict(L=L,a=a,b=b,bubble_identity_error=be,slope_error=se,rotation_eigen_error=rot,seagull=seagull,bubble=bubble))
    return rows

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));OUT.mkdir(parents=True,exist_ok=True);checks=validate();rows=[];sources=[];curves=[]
    for L in (2,5,10,20):
        for A in ('scalar','angle'):
            a,b,f,B0,Z=model(A,L);assert Z>0
            ucrit=f*B0 if A=='scalar' else 0.
            quartic=1-C*np.exp(-a/L**2)*(1+a/L**2)/(6*a*a) if A=='scalar' else None
            for mu in (0,.05,.2):
                pole=0. if mu==0 else brentq(lambda m2:Z*mu*mu+f*diff(a,b,L,-m2),0,mu*mu*1.1)
                residue=1/(f*slope(a,b,L,-pole));mass=np.sqrt(pole)
                row=dict(arm=A,L=L,a=a,b=b,B0=B0,Z=Z,ucrit=ucrit,untuned_K0=-f*B0 if A=='scalar' else 0.,mu=mu,delta=Z*mu*mu,pole_mass=mass,pole_residue=residue,quartic_local=quartic,dispersion=[[k,np.sqrt(k*k+pole)] for k in (0,.05,.2,1)])
                rows.append(row)
                for p in (.001,.01,.1,.5,1,2):curves.append(dict(arm=A,L=L,mu=mu,p=p,K=Z*mu*mu+f*diff(a,b,L,p*p),K_over_Zp2=(Z*mu*mu+f*diff(a,b,L,p*p))/(Z*p*p)))
                for r in (4,8,16,32,64,128):
                    value,err=response(a,b,L,f,Z,mu,r);expected=residue*np.exp(pole/2-mass*r)/(4*np.pi*r)
                    if A=='scalar':secondary,serr=response(a,b,L,f,Z,mu,r,'loop')
                    else:secondary=serr=None
                    row.setdefault('response',[]).append(dict(r=r,unit=value,error=err,pole_asymptotic=expected,ratio=value/expected,secondary_normalized=secondary,secondary_error=serr))
                print(A,L,mu,Z,mass,row['response'][2]['ratio'],flush=True)
            if A=='scalar':
                for lp,lm in ((1,0),(0,1),(1,1),(2,0)):
                    E=(lp+lm)*I(a,L)/2;J=(lp-lm)*B0/2;sources.append(dict(arm=A,L=L,lp=lp,lm=lm,E=E,J=J,J_over_E=J/E))
            else:
                for off in (0,.5,-.5):
                    E=(I(a,L)+I(b,L))/2;J=(a-b)*off*B0;sources.append(dict(arm=A,L=L,diag=[1,1],off=off,E=E,J=J,J_over_E=J/E,explicit_pinned_anisotropy=bool(off)))
    degenerate=model('angle',10,1,1)[-1];assert degenerate==0
    out=dict(validation=checks,degenerate_Z=degenerate,rows=rows,sources=sources,kernel=curves,complete=True,cpu_s=time.process_time()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (OUT/'records.json').write_text(json.dumps(out,indent=1)+'\n')

if __name__=='__main__':main()
