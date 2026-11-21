"""Independent EX-GC26A2 diagnostics. Run static first, then short boosts."""
import os
for var in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[var]='1'
import argparse, ctypes, hashlib, itertools, json, platform, subprocess, time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter, shift, fourier_shift
from scipy.interpolate import PchipInterpolator
import exgc26a_skyrmion as ref

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/EX-GC26A2'; RAW=ROOT/'local_raw/EX-GC26A2'; NORM=12*np.pi**2
PTR=np.ctypeslib.ndpointer(dtype=np.float64, flags='C_CONTIGUOUS'); LIB=None
def save(name,obj):
    OUT.mkdir(parents=True,exist_ok=True); p=OUT/(name+'.json');tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(obj,indent=1,allow_nan=False)+'\n');tmp.replace(p)
def load():
    global LIB
    RAW.mkdir(parents=True,exist_ok=True); p=RAW/'stencil.so'; src=Path(__file__).with_name('exgc26a2_stencil.cpp');
    if not p.exists() or p.stat().st_mtime<src.stat().st_mtime:
        tmp=p.with_name(f'stencil-{os.getpid()}.so');subprocess.run(['c++','-O3','-shared','-fPIC','-std=c++17',str(src),'-o',str(tmp)],check=True);tmp.replace(p)
    LIB=ctypes.CDLL(str(p));LIB.evaluate.argtypes=[ctypes.c_int,ctypes.c_double,ctypes.c_int,PTR,PTR,PTR];LIB.evaluate.restype=None;LIB.acceleration.argtypes=[ctypes.c_int,ctypes.c_double,PTR,PTR,PTR];LIB.acceleration.restype=None
def evaluate(G,p,mode):
    a=np.ascontiguousarray(np.moveaxis(p,0,-1));g=np.empty_like(a);e=np.zeros(2);LIB.evaluate(G.n,G.h,mode,a,g,e);return e,np.moveaxis(g,-1,0)
def acceleration(G,p,v):
    a=np.ascontiguousarray(np.moveaxis(p,0,-1));b=np.ascontiguousarray(np.moveaxis(v,0,-1));out=np.empty_like(a);LIB.acceleration(G.n,G.h,a,b,out);return np.moveaxis(out,-1,0)
def geometric(p):
    """Signed regular-value preimage count of normalized affine tetrahedra."""
    a=np.moveaxis(p,0,-1); vals=np.array([[-.91,.23,.27,.17],[-.77,-.43,.29,.33],[-.83,.32,-.37,.22]]);vals/=np.linalg.norm(vals,axis=1)[:,None]
    counts=np.zeros(3,int); min_dot=1.; bad=0; degenerate=0; near=1.; nt=0
    for perm in itertools.permutations(range(3)):
        steps=np.zeros((4,3),int)
        for i,j in enumerate(perm):steps[i+1]=steps[i];steps[i+1,j]+=1
        vertices=np.stack([a[x:x+a.shape[0]-1,y:y+a.shape[1]-1,z:z+a.shape[2]-1] for x,y,z in steps],axis=-1).reshape(-1,4,4)
        orientation=round(np.linalg.det(steps[1:]))
        for i in range(0,len(vertices),32768):
            v=vertices[i:i+32768];nt+=len(v); mid=v.sum(axis=-1);norm=np.linalg.norm(mid,axis=-1);mid/=np.maximum(norm[:,None],1e-300);dots=np.einsum('ni,nij->nj',mid,v);min_dot=min(min_dot,float(dots.min()));bad+=int(np.sum(dots.min(axis=1)<=0))
            d=np.linalg.det(v);ok=abs(d)>1e-13;degenerate+=int((~ok).sum());v=v[ok];d=d[ok]
            if not len(v):continue
            coeff=np.linalg.solve(v,np.broadcast_to(vals.T,(len(v),4,3)));hits=(coeff>0).all(axis=1)
            counts+=np.sum(hits*np.sign(d[:,None]).astype(int)*orientation,axis=0)
            if hits.any():near=min(near,float(np.min(coeff.transpose(0,2,1)[hits])))
    return dict(degrees=counts.tolist(),min_hemisphere_dot=min_dot,bad_tetrahedra=bad,degenerate_tetrahedra=degenerate,total_tetrahedra=nt,min_positive_preimage_coefficient=near)
def validate():
    G=ref.Grid(.25,1);rng=np.random.default_rng(2602);p=G.fix(G.vacuum()+.2*rng.normal(size=(4,G.n,G.n,G.n)));v=ref.tangent(p,rng.normal(size=p.shape))*G.inner;v/=np.linalg.norm(v)
    e,g=evaluate(G,p,0);old=np.array(G.energy_parts(p));oldg=G.grad(p)*G.h**3/NORM
    errors=dict(energy=float(np.max(abs(e-old))),gradient=float(np.max(abs(g-oldg))))
    for m in (0,1):
        e,g=evaluate(G,p,m);eps=1e-5;fd=(evaluate(G,G.fix(p+eps*v),m)[0].sum()-evaluate(G,G.fix(p-eps*v),m)[0].sum())/(2*eps);an=float(np.sum(g*v));errors[f'tangent_{m}']=abs(fd-an)/max(abs(an),1e-12)
    av=ref.tangent(p,rng.normal(scale=.05,size=p.shape))*G.inner;errors['acceleration_max_error']=float(abs(acceleration(G,p,av)-ref.accel(G,p,av)).max())
    chess=G.vacuum();chess[1]=.2*(-1.)**np.indices(chess.shape[1:]).sum(0);chess=G.fix(chess);errors['checkerboard_energy']=float(evaluate(G,chess,1)[0].sum());errors['vacuum_degree']=geometric(G.vacuum());save('validation',errors);print('validation',errors,flush=True);assert max(errors['energy'],errors['gradient'],errors['acceleration_max_error'])<1e-10 and max(errors['tangent_0'],errors['tangent_1'])<1e-5 and errors['checkerboard_energy']>0
def optimize(G,p,mode,maxiter=300):
    # Positive radial penalty fixes redundant norms; analytic chain rule through normalization.
    inner=G.inner; x0=p[:,inner].T.copy().ravel(); base=G.vacuum(); count=0;start=time.process_time()
    def fg(x):
        nonlocal count
        q=x.reshape(-1,4); n=np.linalg.norm(q,axis=1);u=q/n[:,None];base[:,inner]=u.T
        e,g=evaluate(G,base,mode);g=g[:,inner].T; g=(g-u*np.sum(u*g,axis=1)[:,None])/n[:,None]
        # keeps the arbitrary norms near one without altering the constrained minimum
        lam=1e-5; val=e.sum()+lam*np.sum((n-1)**2);g+=2*lam*(n-1)[:,None]*u;count+=1
        return val,g.ravel()
    r=minimize(fg,x0,method='L-BFGS-B',jac=True,options=dict(maxiter=maxiter,maxcor=8,ftol=1e-13,gtol=3e-8,maxls=30));q=r.x.reshape(-1,4);base[:,inner]=(q/np.linalg.norm(q,axis=1)[:,None]).T;e,g=evaluate(G,base,mode);force=ref.tangent(base,g)*inner
    return base.copy(),dict(iterations=int(r.nit),evals=count,success=bool(r.success),message=str(r.message),E2=float(e[0]),E4=float(e[1]),E=float(e.sum()),ratio=float(e[0]/e[1]),projected_gradient_max=float(abs(force).max()),cpu_s=time.process_time()-start)
def resample(G,p,off):
    x=np.empty_like(p)
    for i in range(4):x[i]=shift(p[i],np.array(off)/G.h,order=3,mode='constant',cval=float(i==0),prefilter=True)
    return G.fix(x)
def pinning(G,p,mode):
    rows=[]
    for frac in ((0,0,0),(.25,0,0),(.5,0,0),(.75,0,0),(1,0,0),(.5,.5,0),(.5,.5,.5)):
        f=resample(G,p,np.array(frac)*G.h);e,_=evaluate(G,f,mode);rows.append(dict(fraction=frac,E=float(e.sum())))
    E=rows[0]['E'];return dict(rows=rows,relative_range=(max(x['E'] for x in rows)-min(x['E'] for x in rows))/E,integer_shift_relative=(rows[4]['E']-E)/E)
def compact_pinning(h,L,mode):
    start=time.process_time();G=ref.Grid(h,L);p=np.load(path(h,L,mode));vac=G.vacuum();q=np.clip((G.r-(L-2.0))/1.0,0,1);t=1-(10*q**3-15*q**4+6*q**5);p=G.fix(vac+t*(p-vac));trans=np.fft.fftn(p-vac,axes=(1,2,3));rows=[]
    for frac in ((0,0,0),(.25,0,0),(.5,0,0),(.75,0,0),(1,0,0),(.5,.5,0),(.5,.5,.5),(-.5,0,0)):
        f=G.fix(vac+np.fft.ifftn(np.stack([fourier_shift(trans[i],frac) for i in range(4)]),axes=(1,2,3)).real);e,_=evaluate(G,f,mode);ec,_=evaluate(G,resample(G,p,np.array(frac)*h),mode);rows.append(dict(fraction=frac,E_fourier=float(e.sum()),E_cubic=float(ec.sum())))
    E=rows[0]['E_fourier'];out=dict(h=h,L=L,mode=mode,rows=rows,relative_range_fourier=(max(x['E_fourier'] for x in rows)-min(x['E_fourier'] for x in rows))/E,relative_range_cubic=(max(x['E_cubic'] for x in rows)-min(x['E_cubic'] for x in rows))/E,integer_shift_relative=(rows[4]['E_fourier']-E)/E,geometric_reflection=geometric(p[:,:,::-1,:]),cpu_s=time.process_time()-start);save(f'compact_pin_h{h:g}_L{L:g}_m{mode}',out);print('compact pin',h,L,mode,out['relative_range_fourier'],out['relative_range_cubic'],out['integer_shift_relative'],out['geometric_reflection']['degrees'],flush=True)
def path(h,L,mode):return RAW/f'field_h{h:g}_L{L:g}_m{mode}.npy'
def prepare(h,L,mode,kind='hedgehog',maxiter=300):
    clock=time.process_time();G=ref.Grid(h,L)
    if kind=='lopsided':
        rng=np.random.default_rng(26);p=G.hedgehog(ref.PROF,centre=(.3,-.2,.1),scale=(1.4,.8,1.));noise=np.array([gaussian_filter(rng.normal(size=p.shape[1:]),2.) for _ in range(4)]);p=G.fix(p+.15*noise/noise.std())
    elif mode and path(h,L,0).exists():p=np.load(path(h,L,0))
    else:p=G.hedgehog(ref.PROF)
    p,info=optimize(G,p,mode,maxiter);info.update(h=h,L=L,mode=mode,kind=kind,B_fd=G.degree(p),geometric=geometric(p),pinning=pinning(G,p,mode),actual_L=float(G.ax[-1]));dest=path(h,L,mode) if kind=='hedgehog' else RAW/'lopsided.npy';np.save(dest,p);info['sha256']=hashlib.sha256(dest.read_bytes()).hexdigest();info['cpu_total']=time.process_time()-clock;save(f'static_h{h:g}_L{L:g}_m{mode}_{kind}',info);print('static',h,L,mode,kind,{k:info[k] for k in ('E','ratio','B_fd','iterations','success','cpu_total')},info['geometric']['degrees'],'pin',info['pinning']['relative_range'],flush=True);return info
def radial_pin(h,L,mode):
    G=ref.Grid(h,L);r=json.loads((ROOT/'results/EX-GC26A/radial.json').read_text());rr=np.array(r['profile']['r']);ff=np.array(r['profile']['f']);interp=PchipInterpolator(np.r_[0.,rr],np.r_[np.pi,ff]);cut=rr[-1];k=r['k'];prof=lambda x:np.where(x<=cut,interp(np.minimum(x,cut)),k/np.maximum(x,1e-12)**2)
    rows=[]
    for frac in ((0,0,0),(.25,0,0),(.5,0,0),(.75,0,0),(1,0,0),(.5,.5,0),(.5,.5,.5)):
        p=G.hedgehog(prof,centre=np.array(frac)*h);e,_=evaluate(G,p,mode);rows.append(dict(fraction=frac,E=float(e.sum()),ratio=float(e[0]/e[1])))
    return dict(h=h,L=L,mode=mode,rows=rows,relative_range=(max(x['E'] for x in rows)-min(x['E'] for x in rows))/rows[0]['E'])
def boost(h,axis,dtfac,L=5):
    start=time.process_time(); G=ref.Grid(h,L);p=np.load(path(h,L,0));E0=G.energy(p);a=np.array(axis,dtype=float);a/=np.linalg.norm(a); v=ref.tangent(p,-.1*sum(a[i]*np.gradient(p,h,axis=i+1) for i in range(3)))*G.inner;dt=dtfac*h;T=1.5;n=int(round(T/dt));rows=[];sample=max(1,round(.15/dt))
    def measure(t,p,v):
        H,P=ref.hamiltonian(G,p,v);dens=np.zeros(p.shape[1:]);
        for s in (-1,1):
            dp,D=G.strain(p,s);tr=np.trace(D,axis1=0,axis2=1);dens+=.5*(tr+.5*(tr**2-(D*D).sum((0,1))))
        cen=[float(np.sum(x*dens)/np.sum(dens)) for x in G.X];return dict(t=t,H=H,P=P,B=G.degree(p),centre=cen)
    for i in range(n+1):
        if i%sample==0 or i==n:rows.append(measure(i*dt,p,v));print('boost',h,axis,dtfac,round(i*dt,3),rows[-1]['P'],flush=True)
        if i==n:break
        k1v=acceleration(G,p,v);k1p=v;k2p=v+dt*k1v/2;k2v=acceleration(G,G.fix(p+dt*k1p/2),k2p);k3p=v+dt*k2v/2;k3v=acceleration(G,G.fix(p+dt*k2p/2),k3p);k4p=v+dt*k3v;k4v=acceleration(G,G.fix(p+dt*k3p),k4p);p=G.fix(p+dt*(k1p+2*k2p+2*k3p+k4p)/6);v=ref.tangent(p,v+dt*(k1v+2*k2v+2*k3v+k4v)/6)*G.inner
    P0=np.array(rows[0]['P']);dP=max(np.linalg.norm(np.array(r['P'])-P0) for r in rows);rec=dict(h=h,L=L,axis=axis,dt=dt,E0=E0,rows=rows,max_dP_over_E0=float(dP/E0),max_dP_over_P0=float(dP/np.linalg.norm(P0)),max_dH_over_H0=max(abs(r['H']/rows[0]['H']-1) for r in rows),geometric_end=geometric(p),cpu_s=time.process_time()-start);tag=f'boost_h{h:g}_axis{"".join(map(str,axis))}_dt{dtfac:g}'+(f'_L{L:g}' if L!=5 else '');save(tag,rec);return rec
def main():
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['validate','static','radial','boost','lopsided','compact']);ap.add_argument('--h',type=float,default=.2);ap.add_argument('--L',type=float,default=5);ap.add_argument('--mode',type=int,default=0);ap.add_argument('--axis',nargs=3,type=int,default=[1,0,0]);ap.add_argument('--dtfac',type=float,default=.05);ap.add_argument('--maxiter',type=int,default=300);a=ap.parse_args();load()
    if a.phase=='validate':validate()
    elif a.phase=='radial':
        t=time.process_time();rec=radial_pin(a.h,a.L,a.mode);rec['cpu_s']=time.process_time()-t;save(f'radial_pin_h{a.h:g}_L{a.L:g}_m{a.mode}',rec);print(rec,flush=True)
    elif a.phase=='boost':boost(a.h,a.axis,a.dtfac,a.L)
    elif a.phase=='compact':compact_pinning(a.h,a.L,a.mode)
    else:prepare(a.h,a.L,a.mode,'lopsided' if a.phase=='lopsided' else 'hedgehog',a.maxiter)
if __name__=='__main__':main()
