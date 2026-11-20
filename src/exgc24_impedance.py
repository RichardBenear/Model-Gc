"""EX-GC24A, prescribed-sink impedance and independent wave integration."""
import os
for var in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[var]='1'
import json, time, platform
from pathlib import Path
import numpy as np
from scipy.fft import fftn, ifftn
from scipy.ndimage import map_coordinates

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/EX-GC24'
BASE=dict(S=5.32,sigma=.8,c=2.,tau=20.,L=32.,h=.5)

def grid(S,sigma,c,tau,L,h):
    n=round(L/h); x=(np.arange(n)-n//2)*h
    s=np.exp(-(x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2)/(2*sigma**2));s*=S/(s.sum()*h**3)
    k=2*np.pi*np.fft.fftfreq(n,h); kz=k[None,None,:]
    q=k[:,None,None]**2+k[None,:,None]**2+kz**2
    sh=fftn(s)*h**3;sh[0,0,0]=0
    return s,sh,kz,q

def spectral(p):
    s,sh,kz,q=grid(**p); c,tau,L=p['c'],p['tau'],p['L']; a=c*c*q; a[0,0,0]=1
    w=2*abs(sh)**2*kz**2/L**3
    m0=float(np.sum(w*c*c/a**2)); b0=m0/tau
    mlow=float(np.sum(w*c*c/a**2*(1-1/(tau*tau*a))))
    rows=[]
    for om in (.001,.003,.01,.03,.1,.3):
        # Stable subtraction, exp(-i omega t); Z=-F/V, m=-Im Z/omega.
        f=np.sum(w*c*c*(om**2+1j*om/tau)/(a*(a-om**2-1j*om/tau)))
        z=f/(1j*om); b=float(z.real); m=float(-z.imag/om)
        rows.append(dict(omega=om,b=b,m=m,ratio=m/b))
    return dict(**p,m_lossless=m0,b_zero=b0,m_zero_finite_damping=mlow,m_infinite=p['S']**2/(6*np.pi**1.5*c*c*p['sigma']),rows=rows)

def steady(p,v):
    s,sh,kz,q=grid(**p); c,tau,L,h=p['c'],p['tau'],p['L'],p['h']; om=v*kz
    den=c*c*q-om**2-1j*om/tau;den[0,0,0]=1
    g=1+ifftn(-(1/tau-1j*om)*sh/den).real/h**3
    j=ifftn(1j*c*c*kz*sh/den).real/h**3
    assert g.min()>0
    wt=s/(s+.001*s.max()); f=float(2*np.sum(s*j/g)*h**3); fl=float(2*np.sum(s*j)*h**3)
    bare=-p['S']*v; bare_local=-float(np.sum(s*wt)*h**3)*v
    return dict(**p,v=v,Gmin=float(g.min()),wake=f,wake_G1=fl,bare=bare,bare_local=bare_local,
        total=f+bare,total_local=float(2*np.sum(s*wt*j/g)*h**3)+bare_local,
        applied_flux=[dict(J=J,force=float(2*np.sum(s*(j+J)/g)*h**3)+bare) for J in (-.02,0,.02)])

def weighted_harmonic(p,om):
    """Linearize the actual 1/G recoil weighting, not only its weak-field proxy."""
    s,sh,kz,q=grid(**p);c,tau,h=p['c'],p['tau'],p['h'];a=c*c*q;a[0,0,0]=1
    g0=1+ifftn(-sh/(tau*a)).real/h**3;j0=ifftn(1j*c*c*kz*sh/a).real/h**3
    ds=-1j*kz*sh;den=a-om**2-1j*om/tau
    dg=ifftn(-(1/tau-1j*om)*ds/den)/h**3;dj=ifftn(1j*c*c*kz*ds/den)/h**3
    dsx=ifftn(ds)/h**3
    force=2*np.sum(dsx*j0/g0+s*dj/g0-s*j0*dg/g0**2)*h**3
    z=force/(1j*om)
    return dict(**p,omega=om,b=float(z.real),m=float(-z.imag/om),Gmin=float(g0.min()))

def integrate(dtmax,A):
    p=dict(BASE,L=12.,h=.75,tau=2.);s,sh,kz,q=grid(**p);c,tau,L=p['c'],p['tau'],p['L'];q[0,0,0]=1
    # Integrate original G and div J equations; Jz reconstructed from longitudinal flux.
    g=-sh/(c*c*q*tau); d=-sh.copy(); q[0,0,0]=0
    omega=.5; period=2*np.pi/omega; nper=int(np.ceil(period/dtmax));dt=period/nper
    def source(t):return sh*np.exp(-1j*kz*A*np.cos(omega*t))
    def rhs(g,d,t):return -d-source(t),c*c*q*g-d/tau
    tlist=[]; forces=[]
    for i in range(6*nper+1):
        t=i*dt
        if i>=4*nper and i<6*nper:
            j=-1j*kz*d/np.where(q>0,q,1)
            f=2*np.sum(np.conj(source(t))*j).real/L**3
            tlist.append(t);forces.append(f)
        if i==6*nper:break
        a,b=rhs(g,d,t);e,f=rhs(g+dt*a/2,d+dt*b/2,t+dt/2);j,k=rhs(g+dt*e/2,d+dt*f/2,t+dt/2);l,m=rhs(g+dt*j,d+dt*k,t+dt)
        g+=dt*(a+2*e+2*j+l)/6;d+=dt*(b+2*f+2*k+m)/6
    t=np.asarray(tlist); F=np.asarray(forces); v=-A*omega*np.sin(omega*t);acc=-A*omega**2*np.cos(omega*t)
    b,m=np.linalg.lstsq(np.stack([-v,-acc],axis=1),F,rcond=None)[0]
    a=c*c*q;a[0,0,0]=1;w=2*abs(sh)**2*kz**2/L**3
    f=np.sum(w*c*c*(omega**2+1j*omega/tau)/(a*(a-omega**2-1j*omega/tau)));z=f/(1j*omega)
    return dict(**p,omega=omega,A=A,dt=dt,b=float(b),m=float(m),b_spectral=float(z.real),m_spectral=float(-z.imag/omega),relative_b_error=float(abs(b/z.real-1)),relative_m_error=float(abs(m/(-z.imag/omega)-1)),fit_rms=float(np.sqrt(np.mean((F+b*v+m*acc)**2))))

def centroid():
    rows=[]
    for K in (.2,.01):
        path=ROOT/'local_raw/EX-GC23'/f'settle_K{K:g}_L6.npz'
        with np.load(path) as dat: W=dat['u'][2]; V=dat['V'][2]
        n=W.shape[0]; h=.25
        for factor in (1,2):
            hh=h/factor; newn=(n-1)*factor+1
            if factor==1: w=W; vel=V
            else:
                idx=np.arange(newn)/factor; coords=np.meshgrid(idx,idx,idx,indexing='ij');w=map_coordinates(W,coords,order=3);vel=map_coordinates(V,coords,order=3)
            z=(np.arange(newn)-(newn-1)/2)*hh; adv0=-np.gradient(w*vel,hh,axis=2);adv1=-np.gradient(w*(vel+.05),hh,axis=2)
            cen=float(np.sum(w*z)/w.sum());speed=float(np.sum((adv1-adv0)*(z-cen))/w.sum())
            rows.append(dict(K=K,h=hh,n=newn,centroid_increment=speed,ratio=speed/.05))
    return rows

def main():
    start=time.process_time();OUT.mkdir(parents=True,exist_ok=True)
    ps=[BASE];sweep=dict(tau=[.25,2.,200.],S=[.266,2.66,10.64],sigma=[.4,1.6],c=[1.,4.],L=[16.,64.]);ps += [dict(BASE,**{k:v}) for k,vs in sweep.items() for v in vs];ps.append(dict(BASE,L=16.,h=.25))
    spec=[spectral(p) for p in ps]; print('spectral done',flush=True)
    st=[steady(dict(BASE,S=S),v) for S in (5.32,.266) for v in (0.,.01,-.01,.025,-.025,.05,-.05)]
    td=[integrate(.02,.02),integrate(.01,.02),integrate(.01,.01)];print('time checks',td,flush=True)
    ct=centroid();print('centroid',ct,flush=True)
    weighted=[weighted_harmonic(dict(BASE,S=S),om) for S in (5.32,.266) for om in (.003,.01,.03)]
    result=dict(spectral=spec,steady=st,time_domain=td,centroid=ct,weighted_harmonic=weighted,cpu_seconds=time.process_time()-start,platform=platform.platform(),numpy=np.__version__)
    (OUT/'diagnostics.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    assert max(max(r['relative_b_error'],r['relative_m_error']) for r in td)<.01

if __name__=='__main__':main()
