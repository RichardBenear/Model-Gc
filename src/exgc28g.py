"""EX-GC28G: static minimally coupled Einstein/Skyrme/secondary scalar, spherical."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,time,resource,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.integrate import solve_bvp,solve_ivp,simpson,cumulative_trapezoid
from exgc28def import dress,C0

OUT=Path(__file__).resolve().parents[1]/'results/EX-GC28G'
C=3*np.pi

def exterior(R,m,d,v,vp,G,details=False):
    # x=1/r, delta unnormalized; scalar flux J constant.
    J=np.exp(d)*(1-2*G*m/R)*R*R*vp
    def rhs(x,y):
        m,d,v=y;N=1-2*G*m*x;f=J*np.exp(-d)/N
        return [-.5*N*f*f/C,-G*x*f*f/C,-f]
    sol=solve_ivp(rhs,(1/R,0),(m,d,v),rtol=2e-10,atol=2e-12,dense_output=details)
    if not sol.success:raise RuntimeError(sol.message)
    return sol.y[:,-1],J,sol

def run(Cmass,q,alpha,G,R=24,tol=2e-7,seed=None):
    r0=1e-4;r=np.linspace(r0,R,900)
    if seed is None:
        base,s=dress(Cmass,q,alpha,R,tol);assert s is not None
        f,fp,v,vp=s.sol(r);ss=np.sin(f);rho=(1-np.cos(f))**2
        e=(r*r+2*ss*ss)*fp*fp+.5*r*r*vp*vp+2*ss*ss+ss**4/r**2+Cmass*r*r*(1-np.cos(f))+r*r*rho*(-q*v+.5*alpha*v*v)
        m=cumulative_trapezoid(e,r,initial=0)/C;d=G*cumulative_trapezoid((2*(r*r+2*ss*ss)*fp*fp+r*r*vp*vp)/r,r,initial=0)/C;d-=d[-1];y0=np.vstack((f,fp,v,vp,m,d))
    else:
        oldR=seed.x[-1];y0=seed.sol(np.minimum(r,oldR))
        if R>oldR:
            yb=seed.y[:,-1];_,J,ext=exterior(oldR,yb[4],yb[5],yb[2],yb[3],G,True);mask=r>oldR;rr=r[mask];mm,dd,vv=ext.sol(1/rr);y0[:,mask]=np.vstack((np.zeros_like(rr),np.zeros_like(rr),vv,J*np.exp(-dd)/(rr**2*(1-2*G*mm/rr)),mm,dd));y0[5]-=y0[5,-1]
    def rhs(r,y):
        f,fp,v,vp,m,d=y;ss=np.sin(f);s2=np.sin(2*f);rho=(1-np.cos(f))**2;drho=2*(1-np.cos(f))*ss;K=r*r+2*ss*ss;Kf=2*s2;N=1-2*G*m/r
        U=2*ss*ss+ss**4/r**2+Cmass*r*r*(1-np.cos(f))+r*r*rho*(-q*v+.5*alpha*v*v)
        Uf=2*s2+2*ss*ss*s2/r**2+Cmass*r*r*ss+r*r*drho*(-q*v+.5*alpha*v*v)
        mp=(N*K*fp*fp+.5*N*r*r*vp*vp+U)/C;dp=G*(2*K*fp*fp+r*r*vp*vp)/(C*r);Np=2*G*(m/r**2-mp/r)
        fpp=(Uf-N*Kf*fp*fp-4*N*r*fp-2*K*(Np+N*dp)*fp)/(2*N*K)
        vpp=rho*(alpha*v-q)/N-(2/r+Np/N+dp)*vp
        return np.vstack((fp,fpp,vp,vpp,mp,dp))
    def bc(ya,yb):
        inf,_,_=exterior(R,yb[4],yb[5],yb[2],yb[3],G)
        return np.array([ya[0]-np.pi-r0*ya[1],yb[0],ya[3],inf[2],ya[4],yb[5]])
    sol=solve_bvp(rhs,bc,r,y0,tol=tol,max_nodes=30000)
    result=dict(C=Cmass,q=q,alpha=alpha,G=G,R=R,tol=tol,success=bool(sol.success),message=sol.message)
    if not sol.success:return result,None
    rr=np.linspace(r0,R,16001);f,fp,v,vp,m,d=sol.sol(rr);N=1-2*G*m/rr;inf,J,ext=exterior(R,m[-1],d[-1],v[-1],vp[-1],G,True);M,di,vi=inf;red=np.exp(d-di)
    s=np.sin(f);rho=(1-np.cos(f))**2;A=N*rr**2*fp*fp;B=s*s;P=2*N*B*fp*fp;T=B*B/rr**2;F=.5*N*rr**2*vp*vp;V=Cmass*rr**2*(1-np.cos(f));I=-q*rr**2*rho*v;Pol=.5*alpha*rr**2*rho*v*v;U=V+I+Pol;dens=A+2*B+P+T+F+U
    pr=2*A+2*P+2*F-dens;pt=2*B+P+2*T-dens
    integ=lambda x:float(simpson(x,x=rr)/C)
    Kmass=integ(red*(dens+pr+2*pt));tail=M-m[-1];Q=integ(red*rr**2*rho*(q-alpha*v));Qflux=-J*np.exp(-di)/C
    # Independent lapse readout at large radius; two radii extrapolate leading 1/r contamination.
    def lapse_mass(r):
        mm,dd,_=ext.sol(1/r);return -np.expm1(2*(dd-di)+np.log1p(-2*G*mm/r))*r/(2*G) if G else M
    lm1,lm2=lapse_mass(128*R),lapse_mass(256*R);lm=2*lm2-lm1
    # Boundary pressure produces finite wall error; exact vacuum scalar tail has zero Tolman density.
    boundary_komar=red[-1]*(m[-1]+rr[-1]*pr[-1]/C)
    result.update(M_ADM=float(M),M_Komar=Kmass,M_lapse=float(lm),M_boundary_komar=float(boundary_komar),M_integral=integ(dens)+float(tail),E2=integ(A+2*B),E4=integ(P+T),EV=integ(V),Efield=integ(F)+float(tail),Einteraction=integ(I),Epolarization=integ(Pol),secondary_Q=Q,secondary_flux=float(Qflux),Q_over_M=Q/M,stress_integral=integ(red*(pr+2*pt)),core_only=integ(A+2*B+P+T+V),tail=float(tail),min_N=float(N.min()),max_residual=float(max(sol.rms_residuals)),phi_infinity=float(vi),phi_center=float(v[0]),translation_flat=(2*integ(A+2*B)+4*integ(P+T)+2*(integ(F)+float(tail)))/3 if G==0 else None,radial_monotone=bool(np.max(fp)<1e-6))
    return result,sol

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));OUT.mkdir(parents=True,exist_ok=True);rows=[]
    def checkpoint(complete=False):
        data=dict(rows=rows,complete=complete,cpu_s=time.process_time()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
        tmp=OUT/'records.tmp';tmp.write_text(json.dumps(data,indent=1,allow_nan=False)+'\n');tmp.replace(OUT/'records.json')
    for cm in (C0,0.):
        seeds={}
        for q,a in [(0,0),(0,5),(.25,0),(-.25,0),(.25,1),(-.25,1),(.25,5),(-.25,5),(.25,-.25),(.5,1)]:
            seed=None
            for G in (0.,.002,.01,.03):
                r,seed=run(cm,q,a,G,seed=seed);rows.append(r);checkpoint();print(cm,q,a,G,r.get('M_ADM'),r.get('M_Komar'),r.get('secondary_Q'),flush=True)
                if not r['success'] or r['min_N']<=0:break
                seeds[q,a,G]=seed
        for q,a,G in [(0,0,.03),(.25,0,.03),(.25,5,.03),(.25,-.25,.03),(.25,1,.01)]:
            for R,tol in [(48,2e-7),(24,5e-8)]:r,_=run(cm,q,a,G,R,tol,seed=seeds.get((q,a,G)));rows.append(r);checkpoint();print('refine',cm,q,a,G,R,r['success'],flush=True)
    checkpoint(True)

if __name__=='__main__':main()
