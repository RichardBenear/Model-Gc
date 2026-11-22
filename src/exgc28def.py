"""EX-GC28D/E/F: heat kernel, dressed Skyrmion, spherical response shells."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,time,resource,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.special import roots_legendre,exp1
from scipy.integrate import solve_bvp,simpson,quad
from scipy.linalg import eigh_tridiagonal
from exgc28a import relax
ROOT=Path(__file__).resolve().parents[1]; C0=.17839644285107545

def save(stage,data,start):
    p=ROOT/f'results/EX-GC28{stage}';p.mkdir(parents=True,exist_ok=True)
    data.update(cpu_s=time.process_time()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    (p/'records.json').write_text(json.dumps(data,indent=1,allow_nan=False)+'\n')

def heat(t,r,kind,xi=0):
    u=t/r**2;l=np.arange(int(np.ceil(np.sqrt(45/u)))+5,dtype=float)
    if kind=='scalar':lam=l*(l+3)+12*xi;g=(l+1)*(l+2)*(2*l+3)/6
    else:lam=(l+2)**2;g=4*(l+1)*(l+2)*(l+3)/3
    return 6*u*u*np.dot(g,np.exp(-u*lam))

def ward():
    # Linearized Einstein tensor in Euclidean flat space, Fourier convention.
    p=np.array([.7,-.3,.4,1.]);p2=p@p; eye=np.eye(4)
    def G(h):return .5*(p2*h-np.outer(p,h@p)-np.outer(h@p,p)+np.outer(p,p)*np.trace(h)+eye*((p@h@p)-p2*np.trace(h)))
    rng=np.random.default_rng(28);v=rng.normal(size=4);g=np.outer(p,v)+np.outer(v,p);a=rng.normal(size=(4,4));a=(a+a.T)/2;P=eye-np.outer(p,p)/p2;h=P@a@P;h-=P*np.trace(h)/3
    return dict(gauge_null=float(np.linalg.norm(G(g))),ward=float(np.linalg.norm(p@G(a))),tt_ratio=float(np.sum(h*G(h))/np.sum(h*h)),tt_expected=p2/2,nonmetric_mass_gauge_norm=float(np.linalg.norm(G(g)+.1*g)))

def D():
    start=time.process_time();rows=[];fits=[];n,w=roots_legendre(96)
    for kind,xi in [('scalar',0),('scalar',1/6),('scalar',1/3),('dirac',0)]:
        a0=1 if kind=='scalar' else 4;a1=1/6-xi if kind=='scalar' else -1/3;sg=1 if kind=='scalar' else -1
        u=np.array([.001,.002,.004,.008]);v=np.array([heat(t,1,kind,xi) for t in u]);co=np.polynomial.polynomial.polyfit(u,v,3);fits.append(dict(kind=kind,xi=xi,a0=float(co[0]),a1=float(co[1]/12),expected_a0=a0,expected_a1=a1,heat=v.tolist()))
        for L in (2,5,10,20):
            for M in (1,2):
                lo=np.log(1/L**2);hi=np.log(40/M**2);t=np.exp((n+1)*(hi-lo)/2+lo);ww=w*(hi-lo)/2
                I1=L**2*np.exp(-M*M/L**2)-M*M*exp1(M*M/L**2)
                I0=quad(lambda z:np.exp(-M*M*z)/z**3,1/L**2,np.inf,epsabs=1e-9)[0]
                ref=sg*a1*I1/(16*np.pi**2);c0=-sg*a0*I0/(32*np.pi**2);vals=[]
                for r in (8,16,32):
                    R=12/r**2;hh=np.array([heat(tt,r,kind,xi) for tt in t]);val=sg*np.dot(ww,np.exp(-M*M*t)*(hh-a0)/t**2)/(16*np.pi**2*R);vals.append(float(val))
                ext=(4*vals[2]-vals[1])/3
                rows.append(dict(kind=kind,xi=xi,L=L,M=M,Mpl2_reference=ref,Mpl2_sphere=vals,Mpl2_extrap=ext,c0=c0))
    comp=[]
    for L in (2,5,10,20):
        for mf in (1,2):
            s=next(r for r in rows if r['kind']=='scalar' and r['xi']==0 and r['L']==L and r['M']==1);f=next(r for r in rows if r['kind']=='dirac' and r['L']==L and r['M']==mf)
            comp.append(dict(L=L,fermion_mass=mf,Mpl2=4*s['Mpl2_reference']+f['Mpl2_reference'],c0=4*s['c0']+f['c0']))
    save('D',dict(heat_fits=fits,rows=rows,mixtures=comp,ward=ward()),start)

def shell(r,a,b,z,beta):
    if b<=a:return np.ones_like(r),np.zeros_like(r),np.zeros_like(r)
    t=np.pi*(r-a)/(b-a);mask=(r>a)&(r<b);s=np.where(mask,np.sin(t)**2,0);ds=np.where(mask,np.pi/(b-a)*np.sin(2*t),0)
    return 1+(z-1)*s,(z-1)*ds,beta*s

def scalar_eigen(R,a,b,z,beta,profile=None,alpha=0):
    # u=r*delta_phi, Dirichlet at origin and finite outer wall; independent fixed-profile diagnostic.
    n=round(R/.05);h=R/n;r=np.arange(1,n)*h;Z,Zp,B=shell(r,a,b,z,beta);links=shell((np.arange(n)+.5)*h,a,b,z,beta)[0]
    U=B+Zp/r
    if profile is not None:U=U+alpha*(1-np.cos(profile.sol(r)[0]))**2
    d=(links[:-1]+links[1:])/h**2+U;o=-links[1:-1]/h**2
    return float(eigh_tridiagonal(d,o,select='i',select_range=(0,0),eigvals_only=True)[0])

def dress(C,q,alpha,R=24,tol=2e-7,a=0,b=0,zshell=1,beta=0,previous=None):
    r0=1e-4;r=np.linspace(r0,R,1000);f,*_=relax(.05,R,1,1,C);rg=np.linspace(0,R,len(f));guess=np.interp(r,rg,f)
    y0=np.vstack((guess,np.gradient(guess,r),q*.4/np.sqrt(1+r*r),-q*.4*r/(1+r*r)**1.5)) if previous is None else previous.sol(r)
    def rhs(r,y):
        f,fp,v,vp=y;s=np.sin(f);s2=np.sin(2*f);rho=(1-np.cos(f))**2;drho=2*(1-np.cos(f))*s;K=r*r+2*s*s;Kf=2*s2;Z,Zp,B=shell(r,a,b,zshell,beta)
        Uf=2*s2+2*s*s*s2/r**2+C*r*r*s+r*r*drho*(-q*v+.5*alpha*v*v)
        fpp=(Uf-Kf*fp*fp-4*r*fp)/(2*K);vpp=((alpha*rho+B)*v-q*rho)/Z-(2/r+Zp/Z)*vp
        return np.vstack((fp,fpp,vp,vpp))
    def bc(ya,yb):return np.array([ya[0]-np.pi-r0*ya[1],yb[0],ya[3],yb[3]+yb[2]/R])
    sol=solve_bvp(rhs,bc,r,y0,tol=tol,max_nodes=30000)
    if not sol.success:return dict(C=C,q=q,alpha=alpha,R=R,tol=tol,a=a,b=b,zshell=zshell,beta=beta,success=False,message=sol.message),None
    rr=np.linspace(r0,R,16001);f,fp,v,vp=sol.sol(rr);s=np.sin(f);rho=(1-np.cos(f))**2;Z,Zp,B=shell(rr,a,b,zshell,beta)
    integ=lambda x:float(simpson(x,x=rr)/(3*np.pi));E2=integ(rr**2*fp**2+2*s*s);E4=integ(2*s*s*fp*fp+s**4/rr**2);EV=integ(C*rr**2*(1-np.cos(f)));Ef=integ(.5*Z*rr**2*vp*vp)+R*v[-1]**2/(6*np.pi);Eint=integ(-q*rr**2*rho*v);Epol=integ(.5*alpha*rr**2*rho*v*v);Eshell=integ(.5*B*rr**2*v*v);E=E2+E4+EV+Ef+Eint+Epol+Eshell;Q0=integ(q*rr**2*rho);Qcore=integ(rr**2*rho*(q-alpha*v));Qshell=integ(-rr**2*B*v);Q=Qcore+Qshell;A=R*v[-1];MT=(2*E2+4*E4+2*Ef)/3
    # Sufficient positivity bound for fixed radial scalar operator (infinite exterior), not coupled Hessian.
    anti=float(simpson(rr*np.maximum(-(alpha*rho+B),0),x=rr)/min(1,zshell));
    out=dict(C=C,q=q,alpha=alpha,R=R,tol=tol,a=a,b=b,zshell=zshell,beta=beta,success=True,E2=E2,E4=E4,EV=EV,Efield=Ef,Einteraction=Eint,Epolarization=Epol,Eshell=Eshell,Etotal=E,dressed_core=E2+E4+EV,translation_mass=MT,virial=(E2+Ef-E4+3*(EV+Eint+Epol+Eshell))/E,Qbare=Q0,Qcore=Qcore,Qshell=Qshell,Q=Q,amplitude=A,flux_error=A-3*np.pi*Q,source_per_energy=Q/E,phi_center=float(v[0]),f_monotone=bool(np.max(fp)<1e-6),residual=float(max(sol.rms_residuals)),anti_bound=anti,scalar_positive_sufficient=anti<1,field_ledger=2*(Ef+Epol+Eshell)+Eint)
    return out,sol

def fixed(base,sol,R,a,b,z,beta):
    r0=1e-4;r=np.linspace(r0,R,1000);C=base['C'];q=base['q'];alpha=base['alpha'];RB=base['R']
    def orig(rr):
        y=sol.sol(np.minimum(rr,RB));return np.where(rr<=RB,(1-np.cos(y[0]))**2*(q-alpha*y[2]),0.)
    def rhs(r,y):
        Z,Zp,B=shell(r,a,b,z,beta);return np.vstack((y[1],(B*y[0]-orig(r))/Z-(2/r+Zp/Z)*y[1]))
    ys=sol.sol(np.minimum(r,RB))[2:];ys[:,r>RB]=np.vstack((base['amplitude']/r[r>RB],-base['amplitude']/r[r>RB]**2));s=solve_bvp(rhs,lambda u,v:np.array([u[1],v[1]+v[0]/R]),r,ys,tol=2e-7,max_nodes=30000)
    if not s.success:return dict(success=False,message=s.message,a=a,b=b,zshell=z,beta=beta,R=R)
    rr=np.linspace(r0,R,16001);v,vp=s.sol(rr);Z,_,B=shell(rr,a,b,z,beta);integ=lambda x:float(simpson(x,x=rr)/(3*np.pi));Qc=integ(rr**2*orig(rr));Qs=integ(-rr**2*B*v);A=R*v[-1];Ef=integ(.5*Z*rr**2*vp*vp)+R*v[-1]**2/(6*np.pi);Ep=integ(.5*B*rr**2*v*v);Ei=integ(-rr**2*orig(rr)*v);anti=float(simpson(rr*np.maximum(-B,0),x=rr)/min(1,z));return dict(success=True,a=a,b=b,zshell=z,beta=beta,R=R,Qcore=Qc,Qshell=Qs,Q=Qc+Qs,amplitude=A,flux_error=A-3*np.pi*(Qc+Qs),phi_center=float(v[0]),Efield=Ef,Eshell=Ep,Einteraction=Ei,field_ledger=2*(Ef+Ep)+Ei,anti_bound=anti,scalar_positive_sufficient=anti<1)

def E():
    start=time.process_time();rows=[]
    for C in (0,C0):
        for q,alpha in [(0,0),(0,5),(.25,0),(-.25,0),(.5,0),(.25,1),(-.25,1),(.5,1),(.25,5),(-.25,5),(.5,5),(.25,-.25),(-.25,-.25),(.5,-.25)]:
            row,_=dress(C,q,alpha);rows.append(row);print('E',C,q,alpha,row.get('Etotal'),row.get('Q'),flush=True)
        for R,tol in [(12,2e-7),(24,5e-8)]:
            for q,alpha in [(.25,0),(.25,5),(.25,-.25)]:row,_=dress(C,q,alpha,R,tol);rows.append(row)
    save('E',dict(rows=rows),start)

def F():
    start=time.process_time();base,s=dress(C0,.25,1);assert s is not None;rows=[]
    for a,b in ((2,4),(4,8)):
        for z,beta in ((1,0),(.1,0),(10,0),(1,.1),(1,1),(1,5),(1,-.01),(1,-.03)):
            for R in ([24,48] if (z,beta) in ((.1,0),(1,1),(1,-.03)) else [24]):
                r=fixed(base,s,R,a,b,z,beta);r['kind']='fixed';r['scalar_min_eigen']=scalar_eigen(R,a,b,z,beta);rows.append(r);r,ss=dress(C0,.25,1,R,a=a,b=b,zshell=z,beta=beta);r['kind']='relaxed';r['scalar_min_eigen']=scalar_eigen(R,a,b,z,beta,ss,1) if ss is not None else None;rows.append(r);print('F',a,b,z,beta,R,rows[-2].get('Q'),r.get('Q'),flush=True)
    save('F',dict(baseline=base,rows=rows),start)

if __name__=='__main__':
    import sys
    resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));{'D':D,'E':E,'F':F}[sys.argv[1]]()
