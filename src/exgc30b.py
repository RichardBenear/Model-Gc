"""EX-GC30B: local homogeneous two-component Bogoliubov spectra, no metric input."""
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[name]='1'
import json,time,platform,resource
from pathlib import Path
import numpy as np

I=np.eye(3); T=np.array([[1.15,.12,0],[.12,.85,.05],[0,.05,1.1]])
P=np.array([.2,-.1,.1]); E=np.eye(3)
dirs=np.array([*E,*[(E[i]+s*E[j])/np.sqrt(2) for i in range(3) for j in range(i+1,3) for s in (1,-1)],np.ones(3)/np.sqrt(3)])
design=np.array([[x*x,y*y,z*z,2*x*y,2*x*z,2*y*z] for x,y,z in dirs])
def sym(b):
    x,y,z,xy,xz,yz=b
    return np.array([[x,xy,xz],[xy,y,yz],[xz,yz,z]])

def cases():
    def c(name,kind='fixed',n=(1,1),g=(1,1),cross=.2,lam=.2,M=(I,I),p=(0,0,0)):
        return dict(name=name,kind=kind,n=np.array(n,dtype=float),g=np.array(g,dtype=float),cross=cross,lam=lam,M=np.array(M,dtype=float),p=np.array(p,dtype=float))
    yield c('decoupled','control',cross=0,lam=0)
    yield c('exchange_split','control',cross=.1,lam=0)
    yield c('baseline')
    for n in (.8,.9,.95,1.05,1.1,1.2): yield c(f'density_{n}',n=(n,n))
    for cross in (.18,.22):yield c(f'cross_{cross}','parameter',cross=cross)
    for lam in (.18,.22):yield c(f'lambda_{lam}','parameter',lam=lam)
    yield c('flow',p=P)
    yield c('anisotropy',M=(T,T))
    yield c('flow_anisotropy',M=(T,T),p=P)
    yield c('density_flow_anisotropy',n=(1.1,1.1),M=(T,T),p=P)
    yield c('relative_density','relative',n=(1.1,.9))
    yield c('relative_g','relative',g=(1.1,.9))
    yield c('relative_tensor','relative',M=(T,I))
    yield c('relative_tensor_flow','relative',M=(T,I),p=P)
    for n in (.8,1.2):yield c(f'retuned_{n}','retuned',n=(n,n),lam=.2*n)

def spectrum(c,k):
    n,g,L,G,M,p=[c[x] for x in ('n','g','lam','cross','M','p')]
    eps=np.einsum('i,aij,j->a',k,M,k)/2; shift=np.einsum('i,aij,j->a',p,M,k)
    off=G*np.sqrt(n.prod()); A=np.diag(eps+g*n+L*np.sqrt(n[::-1]/n)); A[0,1]=A[1,0]=off-L
    B=np.diag(g*n); B[0,1]=B[1,0]=off
    S=np.diag(shift); H=np.block([[A+S,B],[-B,-A+S]])
    ev,V=np.linalg.eig(H); order=np.argsort(ev.real); ev=ev[order]; V=V[:,order]
    if max(abs(ev.imag))>1e-8 or ev[2].real<=0 or ev[1].real>=0:raise ValueError(f'nonstable sampled spectrum: {c["name"]}, {k}, {ev}')
    norms=np.sum(abs(V[:2,2:])**2,axis=0)-np.sum(abs(V[2:,2:])**2,axis=0)
    if min(norms)<=0:raise ValueError('nonpositive Bogoliubov norm')
    return ev[2:].real, float(max(abs(ev.imag)))

def analytic(c,k):
    n=c['n'][0];g=c['g'][0];G=c['cross'];L=c['lam'];M=c['M'][0]; eps=k@M@k/2;shift=c['p']@M@k
    return np.sort(np.sqrt([eps*(eps+2*n*(g+G)),(eps+2*L)*(eps+2*L+2*n*(g-G))])+shift)

def symmetric(c):return np.allclose(c['n'][0],c['n'][1]) and np.allclose(c['g'][0],c['g'][1]) and np.allclose(c['M'][0],c['M'][1])

def run(c):
    n,g,L,G,M,p=[c[x] for x in ('n','g','lam','cross','M','p')]; assert min(n)>0 and min(np.linalg.eigvalsh(M).ravel())>0
    V= -g*n-G*n[::-1]+L*np.sqrt(n[::-1]/n)-np.einsum('i,aij,j->a',p,M,p)/2 # mu=0
    residual=V+g*n+G*n[::-1]-L*np.sqrt(n[::-1]/n)+np.einsum('i,aij,j->a',p,M,p)/2
    assert max(abs(residual))<1e-12
    fits=[]; raw=[]; err=0.; imag=0.
    for top in (.08,.04):
        ks=np.linspace(top/8,top,8); x=ks/top; b=[];d=[];gap=[];quart=[];res=[]
        for j,u in enumerate(dirs):
            both=[]
            for k in ks:
                pm=[]
                for s in (1,-1):
                    w,im=spectrum(c,s*k*u);imag=max(imag,im);pm.append(w)
                    if symmetric(c):err=max(err,float(max(abs(w*w-analytic(c,s*k*u)**2))))
                both.append(pm)
            both=np.array(both); even=both.mean(axis=1)**2; odd=(both[:,0]-both[:,1])/2
            X=np.column_stack([np.ones(8),x*x,x**4]); Y=np.column_stack([x,x**3]); a=np.linalg.lstsq(X,even,rcond=None)[0]; o=np.linalg.lstsq(Y,odd,rcond=None)[0]
            b.append(a[1]/top**2);d.append(o[0]/top);gap.append(a[0]);quart.append(a[2]/top**4);res.append(float(max(abs((X@a-even).ravel()))))
            raw.append(dict(top=top,direction=j,k=ks.tolist(),omega=both.tolist()))
        b=np.array(b);d=np.array(d);bs=np.linalg.lstsq(design[:-1],b[:-1],rcond=None)[0];ds=np.linalg.lstsq(dirs[:-1],d[:-1],rcond=None)[0]
        B=np.array([sym(bs[:,i]) for i in range(2)]); drift=ds.T
        scale=np.mean(np.linalg.norm(B,axis=(1,2)));speed=np.sqrt(np.mean(np.trace(B,axis1=1,axis2=2))/3)
        fits.append(dict(top=top,B=B,drift=drift,gap2=np.mean(gap,axis=0),quartic=np.array(quart),quadratic_mismatch=np.linalg.norm(B[0]-B[1])/scale,drift_mismatch=np.linalg.norm(drift[0]-drift[1])/speed,holdout_B=np.max(abs(design[-1]@bs-b[-1]))/(scale/np.sqrt(3)),holdout_drift=np.max(abs(dirs[-1]@ds-d[-1]))/speed,fit_residual=max(res)))
    a,b=fits;scale=np.mean(np.linalg.norm(b['B'],axis=(1,2)));speed=np.sqrt(np.mean(np.trace(b['B'],axis1=1,axis2=2))/3)
    uncertainty_B=np.max(np.linalg.norm(a['B']-b['B'],axis=(1,2)))/scale;uncertainty_drift=np.max(np.linalg.norm(a['drift']-b['drift'],axis=1))/speed
    analytic_B_error=None
    if symmetric(c):
        target=np.array([n[0]*(g[0]+G)*M[0],(n[0]*(g[0]-G)+2*L)*M[0]])
        # two gapless branches sorted by frequency, otherwise gapless first
        if L==0 and G>0:target=target[::-1]
        analytic_B_error=float(np.max(np.linalg.norm(target-b['B'],axis=(1,2)))/scale)
    hi=[]
    for j,u in enumerate(dirs):
        for k in (.2,.5,1.,2.):
            w,im=spectrum(c,k*u);imag=max(imag,im); quad=np.sqrt(np.maximum(0,b['gap2']+k*k*np.einsum('i,aij,j->a',u,b['B'],u)))+k*(b['drift']@u)
            hi.append(dict(direction=j,k=k,omega=w,quadratic=quad))
    return dict(case=c,V=V,stationarity_error=max(abs(residual)),analytic_omega2_error=err if symmetric(c) else None,analytic_B_error=analytic_B_error,max_imag=imag,fits=fits,uncertainty_B=uncertainty_B,uncertainty_drift=uncertainty_drift,raw=raw,high_k=hi)

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));rows=[]
    for c in cases():
        try:row=run(c)
        except ValueError as e:row=dict(case=c,error=str(e))
        rows.append(row); print(c['name'],row.get('error',''),[(f['quadratic_mismatch'],f['drift_mismatch']) for f in row.get('fits',[])],flush=True)
    out=Path(__file__).resolve().parents[1]/'results/EX-GC30B';out.mkdir(parents=True,exist_ok=True)
    def conv(x):return x.tolist() if isinstance(x,np.ndarray) else x.item() if isinstance(x,np.generic) else None
    for i in range(0,len(rows),5):
        (out/f'cases_{i//5}.json').write_text(json.dumps(dict(rows=rows[i:i+5],python=platform.python_version(),numpy=np.__version__,cpu_s=time.process_time()-start,directions=dirs),default=conv,separators=(',',':'),allow_nan=False)+'\n')
    controls=[r for r in rows if symmetric(r['case'])]
    assert all('error' not in r and r['analytic_omega2_error']<1e-8 and r['analytic_B_error']<.002 for r in controls), 'analytic validation failed'
    print('CPU seconds',time.process_time()-start)

if __name__=='__main__':main()
