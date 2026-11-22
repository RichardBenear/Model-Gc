"""EX-GC29A: elastic particle form factors and radial mode response; no bulk UV claim."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import json,time,platform,resource
from pathlib import Path
import numpy as np
import scipy
from scipy.interpolate import CubicSpline
from scipy.integrate import simpson
from scipy.linalg import eigh_tridiagonal
from scipy.optimize import brentq
from exgc28a import canonical,C0
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC29A'

def modes(h,R,C,old,eta=0.):
    a=np.exp(2*eta);c=C*np.exp(4*eta);row=next((x for x in old if x['h']==h and x['R']==R and x['eta']==eta and x['Cbase']==C),None)
    reused=row is not None
    if row is None:row=canonical(h,R,C,eta)
    n=round(R/h);h=R/n;r=np.arange(1,n)*h;f=np.array(row['profile']);v=2/r**2+c/(2*a);d=2/h**2+v+np.array(row['delta_potential']);o=np.full(n-2,-1/h**2)
    ev,V=eigh_tridiagonal(d,o);evv,Vv=eigh_tridiagonal(2/h**2+v,o);assert min(ev)>0 and min(evv)>0
    om=np.sqrt(ev);err=float(np.max(abs(om-np.array(row['frequencies']))));assert err<1e-8
    return dict(h=h,R=R,C=C,eta=eta,A=a,Ceffective=c,reused=reused,spectrum_error=err,r=r,f=f,d=d,o=o,omega=om,V=V,vacomega=np.sqrt(evv),Vvac=Vv,min_lambda=float(ev[0]))

def form(m):
    R=m['R'];h=m['h'];a=m['A'];c=m['Ceffective'];r=np.linspace(0,R,round(R/.0025)+1);sp=CubicSpline(np.r_[0,m['r'],R],np.r_[np.pi,m['f'],0],bc_type='natural');f=sp(r);z=sp(r,1);s=np.sin(f);p=2*np.sin(f/2)**2
    safe=np.where(r>0,r,1.);weights={'baryon':-2/np.pi*s*s*z,'P1':r*r*p,'P2':r*r*p*p,'energy':a*(r*r*z*z+2*s*s)+2*s*s*z*z+s**4/safe**2+c*r*r*p}
    out=[]
    for name,w in weights.items():
        norm=float(simpson(w,x=r));rms=float(np.sqrt(simpson(w*r*r,x=r)/norm))
        def F(k):return float(simpson(w*np.sinc(k*r/np.pi),x=r)/norm)
        # Scan for first crossing, since form factors can have nodes.
        ks=np.linspace(0,12,241);vals=np.array([F(k)**2 for k in ks]);kcut=[]
        for t in (.5,.1):
            i=np.flatnonzero(vals<=t)[0];kcut.append(brentq(lambda k:F(k)**2-t,ks[i-1],ks[i]))
        out.append(dict(operator=name,norm=norm,rms=rms,k50=kcut[0],k90=kcut[1],k50rms=kcut[0]*rms,gaussian_k50=np.sqrt(3*np.log(2))/rms,values=[dict(k=k,F=F(k),gaussian=float(np.exp(-k*k*rms*rms/6)),resolved=h*k<=.5) for k in (0,.5,1,2,4,8,16,32,64)]))
    return out

def response(m):
    r=m['r'];h=m['h'];f=m['f'];K=m['A']*r*r+2*np.sin(f)**2;rows=[]
    for name,p in [('P1',np.sin(f)),('P2',2*(2*np.sin(f/2)**2)*np.sin(f))]:
        J=r*r*p/np.sqrt(6*np.pi*K);norm=h*np.dot(J,J);quad=h*(np.dot(m['d']*J,J)+2*np.dot(m['o']*J[:-1],J[1:]));
        for system,om,V in [('defect',m['omega'],m['V']),('vacuum',m['vacomega'],m['Vvac'])]:
            w0=h*(V.T@J)**2;assert abs(w0.sum()/norm-1)<1e-10
            if system=='defect':assert abs(np.dot(w0,om**2)/quad-1)<1e-8
            w1=np.empty_like(w0)
            for i in range(0,len(om),64):
                j=slice(i,i+64);drive=J[:,None]*np.sinc(r[:,None]*om[None,j]/np.pi);w1[j]=h*np.sum(V[:,j]*drive,axis=0)**2
            for drive,w in [('uniform',w0),('on_shell',w1)]:
                cuts=[]
                for cut in (2,4,8,16,32,64):
                    mask=om<=cut;ws=w[mask];freq=om[mask];mom=[];q90=[]
                    for pwr in (0,2,4):
                        wt=ws*freq**pwr;tot=float(wt.sum());mom.append(tot);q90.append(float(freq[min(np.searchsorted(np.cumsum(wt),.9*tot),len(freq)-1)]) if tot>0 else None)
                    cuts.append(dict(cut=cut,resolved=h*cut<=.5,moments=mom,q90=q90,modes=int(mask.sum())))
                rows.append(dict(operator=name,system=system,drive=drive,norm=float(norm),parseval=float(w0.sum()),quadratic_form=float(quad),uniform_frequency2=float(np.dot(w0,om**2)),cuts=cuts,total_unweighted=float(w.sum()),max_mode_frequency=float(om[-1])))
    return rows

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));OUT.mkdir(parents=True,exist_ok=True);path=ROOT/'results/EX-GC28A/canonical.json';old=json.loads(path.read_text()) if path.exists() else [];rows=[]
    def run(h,R,C,eta=0):
        m=modes(h,R,C,old,eta);out={k:v for k,v in m.items() if k not in ('r','f','d','o','omega','V','vacomega','Vvac')};out['form']=form(m);out['response']=response(m);rows.append(out);print('case',h,R,C,eta,[(x['operator'],round(x['k50'],4)) for x in out['form']],flush=True)
    for C in (0.,C0):
        for R in (12,24):
            for h in (.1,.05,.025,.0125):run(h,R,C)
        for eta in (-np.log(1.5),np.log(1.5)):
            lam=np.exp(-eta);run(.025*lam,24*lam,C,eta)
    meta=dict(cpu_s=time.process_time()-start,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__)
    for C,label in [(0.,'massless'),(C0,'massive')]:
        (OUT/f'{label}.json').write_text(json.dumps(dict(rows=[r for r in rows if r['C']==C],**meta),separators=(',',':'),allow_nan=False,default=lambda x:x.item() if isinstance(x,np.generic) else x)+'\n')

if __name__=='__main__':main()
