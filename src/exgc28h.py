"""EX-GC28H: symmetric scalar loop and fixed-core screening / fifth force."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,time,resource,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.linalg import solve_banded
from scipy.optimize import brentq
from exgc28c import I,B,diff,slope
from exgc28def import dress,C0
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC28H';NORM=3*np.pi

def loops():
    rows=[];sources=[]
    for a in ((1.,1.),(1.,4.)):
        for L in (2.,5.,10.,20.):
            tad=sum(g*I(g,L)/2 for g in a);u=sum(g*g*B(g,g,L)/2 for g in a);Z=sum(g*g*slope(g,g,L)/2 for g in a)
            cubic=sum(g**3*np.exp(-g/L**2)/(32*np.pi**2*g) for g in a)
            def t(s):return sum(g*I(g*(1+s),L)/2 for g in a)
            h=1e-3;fd=(t(h)-2*t(0)+t(-h))/h**2;assert abs(fd/cubic-1)<1e-5
            for mu in (0.,.05,.2):
                pole=brentq(lambda z:Z*mu*mu+sum(g*g*diff(g,g,L,-z)/2 for g in a),0,1.1*mu*mu) if mu else 0.
                rows.append(dict(a=list(a),L=L,tadpole=tad,ucrit=u,untuned_K0=-u,Z=Z,cubic=cubic,cubic_fd=fd,mu=mu,delta=Z*mu*mu,pole_mass=np.sqrt(pole)))
            for lam in ((1,0),(0,1),(1,1),(2,0)):
                E=sum(v*I(g,L)/2 for g,v in zip(a,lam));J=sum(g*v*B(g,g,L)/2 for g,v in zip(a,lam));sources.append(dict(a=list(a),L=L,lam=list(lam),E=E,J=J,J_over_E=J/E))
    return rows,sources

def fifth():
    p=json.loads((ROOT/'results/EX-GC28G/records.json').read_text());x=[r for r in p['rows'] if r['success'] and r['C']>0 and r['G']==.03 and r['R']==24 and r['tol']==2e-7];rows=[]
    for r in x:
        qm=r['secondary_Q']/r['M_ADM'];rows.append(dict(q=r['q'],alpha=r['alpha'],Q=r['secondary_Q'],M=r['M_ADM'],Q_over_M=qm,ratios=[dict(G=G,ratio=NORM*qm*qm/G) for G in (.03,1e-6,1e-12,1e-30)]))
    bounds=[]
    for G in (.03,1e-6,1e-12,1e-30):
        for eps in (1e-5,1e-13):
            for r in rows:
                ratio=NORM*r['Q_over_M']**2/G
                mr=0. if ratio<=eps else brentq(lambda t:np.log1p(t)-t-np.log(eps/ratio),0,200)
                bounds.append(dict(G=G,epsilon=eps,max_Q_over_M=np.sqrt(eps*G/NORM),q=r['q'],alpha=r['alpha'],unsuppressed_ratio=ratio,mr_required=mr))
    return rows,bounds

def screening(C,q,alpha,mu,R,h,base,s):
    n=round(R/h);r=np.arange(1,n+1)*h;f=s.sol(r)[0];rho=(2*np.sin(f/2)**2)**2;V=mu*mu+alpha*rho
    d=2/h**2+V;d[-1]=1/h**2+mu/h+.5*V[-1];o=np.full(n-1,-1/h**2);rhs=q*r*rho;rhs[-1]*=.5
    ab=np.zeros((3,n));ab[0,1:]=o;ab[1]=d;ab[2,:-1]=o;u=solve_banded((1,1),ab,rhs);w=np.ones(n)*h;w[-1]*=.5
    # Integrate field energy with exterior tail using u=r phi. Integration by parts cancels -u_R²/R.
    Efield=(np.sum(np.diff(np.r_[0,u])**2)/h+mu*mu*np.dot(w,u*u)+mu*u[-1]**2)/(2*NORM)
    Epol=.5*alpha*np.dot(w,rho*u*u)/NORM;Eint=-q*np.dot(w,r*rho*u)/NORM;E=base['Etotal']+Efield+Epol+Eint
    Q=u[-1]*np.exp(mu*R)/NORM;S=rho*(q*r-alpha*u);weight=np.sinh(mu*r)/mu if mu else r;moment=np.dot(w,weight*S)/NORM
    # Exact discrete summation-by-parts with homogeneous lattice growing solution.
    kap=2*np.arcsinh(mu*h/2)/h;test=np.sinh(kap*r)*h/np.sinh(kap*h) if mu else r
    # Sum full discrete source against test; boundary term at last row.
    discrete=np.dot(w,test*S)/NORM
    boundary=u[-1]*(test[-1]*(1/h+mu+.5*h*mu*mu)-test[-2]/h)/NORM
    return dict(C=C,q=q,alpha=alpha,mu=mu,R=R,h=h,Q=Q,Q_moment=moment,Q_discrete=discrete,Q_discrete_boundary=boundary,energy=E,Ecore=base['Etotal'],Efield=Efield,Epol=Epol,Einteraction=Eint,ledger=2*(Efield+Epol)+Eint,Q_over_E=Q/E,phi_center=u[0]/h,ratio_G03=NORM*(Q/E)**2/.03)

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));OUT.mkdir(parents=True,exist_ok=True);l,ls=loops();f,fb=fifth();rows=[]
    for C in (0.,C0):
        for R in (24,48):
            base,s=dress(C,0,0,R);assert s is not None
            for q in (.25,.5):
                for alpha in (0,1,5,25,100,1000,10000,1000000):
                    for mu in (0.,.05,.2):
                        if R==48 and (q!=.25 or alpha not in (0,100,1000000)):continue
                        for h in (.1,.05,.025):
                            if h!=.05 and (q!=.25 or alpha not in (0,100,1000000) or R!=24):continue
                            rows.append(screening(C,q,alpha,mu,R,h,base,s))
    data=dict(loop=l,loop_sources=ls,fifth=f,fifth_bounds=fb,screening=rows,cpu_s=time.process_time()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (OUT/'records.json').write_text(json.dumps(data,indent=1,allow_nan=False)+'\n');print('DONE',len(rows),time.process_time()-start)

if __name__=='__main__':main()
