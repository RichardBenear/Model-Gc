"""EX-GC30B3: finite-scale one-loop velocity locking, analytical controls."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
import json,time,platform,resource
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.linalg import expm

D=np.array([0,1,5,10,20,30,40]);TS=D*np.log(10);PI=np.pi
def yuk_beta(y):
    f,b,g=y
    bf=g*g*(f-b)/(6*PI**2*b*(f+b)**2)
    bb=g*g*(b*b-f*f)/(8*PI**2*b*f**3)
    bg=g**3*(3*b*f*f+2*b*b*f+b**3+4*f**3)/(8*PI**2*b*f**3*(f+b)**2)
    return np.array([bf,bb,bg])
def alpha(y):
    f,b,g=y;return g*g/(16*PI**2*((f+b)/2)**3)
def lin_yuk(eta,a,t,fixed=False):return eta*(np.exp(-14*a*t/3) if fixed else (1+10*a*t)**(-7/15))
def mismatch(c):return (max(c)-min(c))/np.mean(c)
def armY():
    rows=[]; ratioerr=0.;jacerr=0.
    for f,b in ((1,.8),(1,1.2),(1,.99),(1,1.01),(1,1)):
        for a in (.001,.01,.03,.1):
            # normalize coupling by common-reference c0=1, not by a species-specific clock
            y0=np.array([f,b,np.sqrt(16*PI**2*a)])
            sols=[]
            for tol in (1e-9,1e-12):
                sol=solve_ivp(lambda t,y:-yuk_beta(y),(0,TS[-1]),y0,t_eval=TS,rtol=tol,atol=tol*.01)
                assert sol.success and np.all(sol.y>0);sols.append(sol.y.T)
            z=sols[-1];err=np.max(abs(z-sols[0]));assert err<1e-7
            for q in z:
                ff,bb,gg=q;r=bb/ff;bf,bbb,_=yuk_beta(q)
                br=(bbb-r*bf)/ff;eq19=gg*gg*(r-1)*(8*r+6*(1+r)**3)/(48*PI**2*bb*ff**2*(1+r)**2)
                ratioerr=max(ratioerr,abs(br-eq19));assert br*(r-1)>=-1e-15
            eta=z[:,1]/z[:,0]-1;linear=lin_yuk(b/f-1,a,TS)
            assert np.all(np.diff(abs(eta))<=1e-10)
            if b==f:assert max(abs(z[:,0]-f))<1e-12 and max(abs(z[:,1]-b))<1e-12 and max(abs(z[:,2]**2/(16*PI**2)-(a/(1+10*a*TS))))<1e-9
            rows.append(dict(f0=f,b0=b,alpha0=a,values=z,alpha=[alpha(q) for q in z],eta=eta,linear_eta=linear,mismatch=[mismatch(q[:2]) for q in z],refinement=err))
    # finite-difference velocity Jacobian at equality, holding g constant
    for a in (.001,.01,.03,.1):
        y=np.array([1.,1.,np.sqrt(16*PI**2*a)]);h=1e-5;J=np.column_stack([(yuk_beta(y+h*np.eye(3)[i])[:2]-yuk_beta(y-h*np.eye(3)[i])[:2])/(2*h) for i in (0,1)])
        eig=np.sort(np.linalg.eigvals(J));jacerr=max(jacerr,max(abs(eig-np.array([0,14*a/3]))))
    assert ratioerr<1e-12 and jacerr<1e-9
    return rows,dict(ratio_error=ratioerr,jacobian_error=jacerr)

def gauge_matrix(q):
    q=np.asarray(q);L=np.zeros((len(q)+1,len(q)+1));r=8*q*q/3
    for i,p in enumerate(r,1):L[0,i]=p;L[i,0]=2*p;L[0,0]-=p;L[i,i]-=2*p
    return L,float(r.sum())
def evolve_g(x,a,q,t,fixed=False):
    L,b=gauge_matrix(q);s=a*t if fixed or b==0 else np.log1p(b*a*t)/b
    z=expm(L*s)@x;aa=a if fixed else a/(1+b*a*t)
    return z,aa
def armG():
    rows=[];initial=np.array([0.,-.1,.1]);wt=np.array([2.,1.,1.]);maxerr=0.
    for q in ((1,1),(1,.3),(1,.1),(1,.01),(1,0)):
        L,b=gauge_matrix(q);assert np.max(abs(wt@L))<1e-12
        eig=np.sort(np.linalg.eigvals(L).real);assert max(eig)<1e-12
        for a in (.001,.01,.03):
            for mode in ('massless','threshold'):
                out=[];alph=[];threshold=10*np.log(10)
                for t in TS:
                    if mode=='threshold' and t>threshold:
                        xx,aa=evolve_g(initial,a,q,threshold);xx,aa=evolve_g(xx,aa,(q[0],0),t-threshold)
                    else:xx,aa=evolve_g(initial,a,q,t)
                    out.append(xx);alph.append(aa)
                out=np.array(out);assert max(abs(out@wt-initial@wt))<1e-12
                if q[-1]==0:assert max(abs(out[:,2]-initial[2]))<1e-12
                if mode=='threshold':assert max(abs(out[3:,2]-out[3,2]))<1e-12
                if mode=='massless':
                    sol=solve_ivp(lambda t,y:np.r_[y[-1]*L@y[:-1],-b*y[-1]**2],(0,TS[-1]),np.r_[initial,a],t_eval=TS,rtol=1e-12,atol=1e-14)
                    assert sol.success;err=max(abs(sol.y.T-np.column_stack([out,alph])).ravel());maxerr=max(maxerr,err);assert err<1e-9
                rows.append(dict(q=q,alpha0=a,mode=mode,values=out,alpha=alph,mismatch=[mismatch(1+z) for z in out],eigenvalues=eig))
        out=np.array([evolve_g(initial,.03,q,t,fixed=True)[0] for t in TS]);rows.append(dict(q=q,alpha0=.03,mode='constant_alpha',values=out,alpha=np.full(len(TS),.03),mismatch=[mismatch(1+z) for z in out],eigenvalues=eig))
    return rows,maxerr

def armT():
    # Linearized operator mixing restricted to constant symmetric traceless sources;
    # same species-space map for each Lorentz tensor component (not a nonlinear tensor beta function).
    M=np.array([[-2/3,2/3],[4,-4.]])
    # covariant tensor with Minkowski trace zero, order t,x,y,z
    basis=[np.diag([3.,1,1,1]),np.diag([0.,1,-1,0]),np.array([[0,0,0,0],[0,0,1,0],[0,1,0,0],[0,0,0,0]]),np.array([[0,1,0,0],[1,0,0,0],[0,0,0,0],[0,0,0,0]])]
    metric=np.diag([-1,1,1,1]);rows=[]
    for a in (.001,.01,.03):
        for name,B in zip(('isotropic_traceless','shear_xy_diag','shear_xy_offdiag','tilt_tx'),basis):
            assert abs(np.sum(metric*B))<1e-12;B=B/np.linalg.norm(B);initial=np.array([-.005,.005]);vals=[]
            for t in TS:
                s=np.log1p(10*a*t)/10;coeff=expm(M*s)@initial;assert abs((coeff[1]-coeff[0])-lin_yuk(.01,a,t))<1e-12;assert np.max(abs(expm(M*s)@np.ones(2)-1))<1e-12
                vals.append(dict(coeff=coeff,tensors=coeff[:,None,None]*B,relative=np.linalg.norm((coeff[1]-coeff[0])*B)))
            rows.append(dict(name=name,alpha0=a,values=vals))
    return rows

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));Y,Ychecks=armY();G,Gcheck=armG();T=armT();burden=[]
    for a in (.001,.01,.03,.1):
        for target in (.01,1e-6,1e-15):
            running=np.expm1((15/7)*np.log(.1/target))/(10*a*np.log(10));fixed=np.log(.1/target)/(14*a/3*np.log(10));burden.append(dict(alpha0=a,target=target,decades_running=running,decades_fixed=fixed))
    out=Path(__file__).resolve().parents[1]/'results/EX-GC30B3';out.mkdir(parents=True,exist_ok=True)
    meta=dict(decades=D,cpu_s=time.process_time()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    def conv(z):
        if isinstance(z,np.ndarray):return z.tolist()
        if isinstance(z,np.generic):return z.item()
        raise TypeError(type(z))
    for name,data in [('yukawa',dict(rows=Y,checks=Ychecks,burden=burden)),('gauge',dict(rows=G,check=Gcheck)),('tensor',dict(rows=T))]:
        (out/(name+'.json')).write_text(json.dumps(dict(**meta,**data),default=conv,separators=(',',':'),allow_nan=False)+'\n')
    print('CPU',meta['cpu_s'],'checks',Ychecks,Gcheck)
    for r in Y:
        if r['b0'] in (.8,1.2): print('Y',r['b0'],r['alpha0'],r['mismatch'][-1])
    for r in G:
        if r['alpha0']==.03:print('G',r['q'],r['mode'],r['mismatch'][-1])

if __name__=='__main__':main()
