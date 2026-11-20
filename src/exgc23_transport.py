"""EX-GC23 exploratory material transport. Run from root: python src/exgc23_transport.py --workers 4."""
import os
for var in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):
    os.environ[var]='1'
import argparse, ctypes, hashlib, json, platform, subprocess, tempfile, time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import numpy as np
import scipy
from scipy.ndimage import gaussian_filter
import exgc22_drag as g22
import exgc20_wave as wave
import ex008_longrange as e
import exgc09_lattice as g9

OUT=Path(e.ROOT)/'results/EX-GC23'
RAW=Path(e.ROOT)/'local_raw/EX-GC23'
PTR=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
LIB=None
def load():
    global LIB
    tmp=tempfile.TemporaryDirectory(prefix='exgc23-')
    lib=Path(tmp.name)/'kernel.so'
    subprocess.run(['c++','-O3','-shared','-fPIC',str(Path(__file__).with_name('exgc23_kernel.cpp')),'-o',str(lib)],check=True)
    LIB=ctypes.CDLL(str(lib)); LIB.chem.argtypes=[ctypes.c_int]+[ctypes.c_double]*7+[PTR]*3
    LIB.chem.restype=None; LIB.wave_step.argtypes=[ctypes.c_int]+[ctypes.c_double]*5+[PTR]*5; LIB.wave_step.restype=None
    return tmp

def dump(name,obj):
    OUT.mkdir(exist_ok=True,parents=True)
    dest=OUT/(name+'.json'); tmp=dest.with_suffix('.tmp'); tmp.write_text(json.dumps(obj,indent=1,allow_nan=False)+'\n'); tmp.replace(dest)

def integral(s,V,h):
    return np.einsum('ijk,aijk->a',s,V)*h**3

def material(V,h):
    return -np.array([sum(V[a]*np.gradient(v,h,axis=a) for a in range(3)) for v in V])

def conservative(s,V,h):
    # Sealed-face upwind mass flux; transport of momentum and density by same donor.
    dm=np.zeros_like(s); dp=np.zeros_like(V)
    for a in range(3):
        lo=[slice(None)]*3; hi=lo.copy();lo[a]=slice(None,-1);hi[a]=slice(1,None);lo=tuple(lo);hi=tuple(hi)
        vl=(a,*lo);vh=(a,*hi);ll=(slice(None),*lo);hh=(slice(None),*hi)
        vf=.5*(V[vl]+V[vh]); f=vf*np.where(vf>=0,s[lo],s[hi])/h
        p=f*np.where(vf[None]>=0,V[ll],V[hh]);dm[lo]-=f;dm[hi]+=f;dp[ll]-=p;dp[hh]+=p
    adv=np.divide(dp-V*dm,s[None],out=np.zeros_like(V),where=s[None]>1e-30)
    return adv,dm,dp

class Box(g22.Box):
    def __init__(self,K,L):
        super().__init__(.25,L,80,K);self.du=np.empty((3,*self.shape));self.gn=np.empty(self.shape)
    def chem_inplace(self,u,V):
        LIB.chem(self.n,self.h,self.K,e.A,e.B,e.DY,e.DW,e.C,u,V,self.du)
        self.du*=self.dt;u[:3]+=self.du
    def step(self,u,V,J,E):
        # G/s old for both chemistry and medium. Compute both before mutation.
        LIB.chem(self.n,self.h,self.K,e.A,e.B,e.DY,e.DW,e.C,u,V,self.du)
        LIB.wave_step(self.n,self.h,self.dt,self.K,20.,E,u,*J,self.gn)
        self.du*=self.dt;u[:3]+=self.du;u[3]=self.gn
    def terms(self,u,V,push,mode):
        s=self.K*u[3]*u[2]; sig=1.5/self.h; sm=gaussian_filter(s,sig,mode='constant'); sv=np.array([gaussian_filter(s*v,sig,mode='constant') for v in V]);
        R=s/(s+1e-3*s.max())*np.array(push)/20
        if mode=='cons':
            C=(sv-V*sm)/(sm.max()+1e-12); A,dm,dp=conservative(s,V,self.h)
        else:
            C=sv/(sm+1e-12)-V; A=material(V,self.h) if mode=='on' else np.zeros_like(V)
        return s,R,C,A
    def settle(self,T=30):
        _,_,u=g9.init(self.h,float(self.ax[-1]));V=np.zeros((3,*self.shape));u[3]=self.poisson(u[2],u[3]);rows=[]
        n=round(T/self.dt);every=round(1/self.dt)
        for i in range(n+1):
            if i%(self.kf*10)==0 or i==n: u[3]=self.poisson(u[2],u[3])
            if i%self.kf==0 or i==n:
                push=[2*j/u[3] for j in wave.centre_flux([-self.D*g for g in wave.grad_faces(u[3],self.h,self.bc(0))])];V=self.vel(u,V,push)
            if i%every==0 or i==n: rows.append(self.measure(i*self.dt,u,V))
            if i<n: self.chem_inplace(u,V)
        return u,V,rows
    def evolve(self,u0,V0,mode,kick,flux,T):
        u=u0.copy();V=V0.copy();V[2]+=kick;E=-flux/self.D;u[3]+=E*self.Z;J=[np.ascontiguousarray(-self.D*g) for g in wave.grad_faces(u[3],self.h,self.bc(E))]
        # Ledger in units of s*V; the omitted common assigned inertia multiplier is tau=20.
        p0=integral(self.K*u[2]*u[3],V,self.h);previous=p0.copy();last=previous.copy();ledger=np.zeros((4,3));rows=[]
        n=round(T/self.dt);every=round(.5/self.dt)
        for i in range(n+1):
            sample=i%every==0 or i==n;update=i%self.kf==0 and i<n
            if sample or update:
                s=self.K*u[2]*u[3];p=integral(s,V,self.h);ledger[3]+=p-last;last=p.copy()
            if sample:
                r=self.measure(i*self.dt,u,V);r.update(scen_z=float((s*self.Z).sum()/s.sum()),P=p.tolist(),ledger=ledger.tolist(),
                    residual=(p-p0-ledger.sum(axis=0)).tolist(),Gmin=float(u[3].min()),Vmax=float(abs(V).max()))
                rows.append(r)
                if not np.isfinite(u).all() or not np.isfinite(V).all() or r['Gmin']<=0 or r['Wmax']<.01 or u[2].min()<-1e-6:
                    return rows,'failed survival/finite/positivity'
            if i==n: break
            if update:
                push=[2*j/u[3] for j in wave.centre_flux(J)];s,R,C,A=self.terms(u,V,push,mode)
                for j,f in enumerate((R,C,A)):ledger[j]+=self.dtV*integral(s,f,self.h)
                V+=self.dtV*(R+C+A);last=integral(s,V,self.h)
            self.step(u,V,J,E)
        return rows,'complete'

def start(K,L):
    t=time.monotonic(); cpu=time.process_time(); b=Box(K,L);u,V,rows=b.settle();RAW.mkdir(exist_ok=True,parents=True)
    name=f'settle_K{K:g}_L{L:g}';path=RAW/(name+'.npz');np.savez_compressed(path,u=u,V=V);hash_=hashlib.sha256(path.read_bytes()).hexdigest()
    dump(name,dict(K=K,L=L,h=b.h,dt=b.dt,rows=rows,sha256=hash_,cpu_s=time.process_time()-cpu,wall_s=time.monotonic()-t))
    print(name,'done',rows[-1], 'cpu',round(time.process_time()-cpu,2),flush=True)
    return name

def run(case):
    name,K,L,mode,kick,flux,T=case;t=time.monotonic();cpu=time.process_time();b=Box(K,L)
    with np.load(RAW/f'settle_K{K:g}_L{L:g}.npz') as dat:u,V=dat['u'],dat['V']
    rows,status=b.evolve(u,V,mode,kick,flux,T)
    dump(name,dict(name=name,K=K,L=L,h=b.h,mode=mode,kick=kick,flux=flux,T=T,dt=b.dt,dtV=b.dtV,D=80,status=status,rows=rows,
        cpu_s=time.process_time()-cpu,wall_s=time.monotonic()-t,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__))
    print(name,status,rows[-1],'cpu',round(time.process_time()-cpu,2),flush=True);return name

CASES=[('strong_off',.2,6,'off',.05,0,20),('strong_on',.2,6,'on',.05,0,20),('weak_off',.01,6,'off',.05,0,20),('weak_on',.01,6,'on',.05,0,20),
       ('flux_K02',.2,6,'on',0,.02,20),('flux_K04',.4,6,'on',0,.02,20),('strong_null',.2,6,'on',0,0,20),('wall_on',.2,9,'on',.05,0,20),('strong_cons',.2,6,'cons',.05,0,20)]

def validate():
    b=g22.Box(.25,1,80,.2);rng=np.random.default_rng(23);u=rng.uniform(.8,1.2,(4,*b.shape));V=rng.normal(0,.05,(3,*b.shape));E=.002
    J=[np.ascontiguousarray(rng.normal(0,.02,g.shape)) for g in wave.grad_faces(u[3],b.h,b.bc(E))];out=np.empty((3,*b.shape));gn=np.empty(b.shape)
    LIB.chem(b.n,b.h,b.K,e.A,e.B,e.DY,e.DW,e.C,u,V,out);err1=float(abs(out-b.chem_numpy(u,V)).max())
    refG=u[3]+b.dt*(-wave.div_faces(J,b.h)-b.K*u[2]*u[3]);refJ=[j+b.dt*(-4*g-j/20) for j,g in zip(J,wave.grad_faces(u[3],b.h,b.bc(E)))];LIB.wave_step(b.n,b.h,b.dt,b.K,20.,E,u,*J,gn)
    err2=max(float(abs(gn-refG).max()),*(float(abs(j-r).max()) for j,r in zip(J,refJ)))
    x,y,z=np.meshgrid(b.ax,b.ax,b.ax,indexing='ij');v=np.array([.1*x,.2*y,.3*z+.05]);err3=float(abs(material(v,b.h)+v*np.array([.1,.2,.3])[:,None,None,None]).max())
    s=rng.uniform(.01,1,b.shape);adv,dm,dp=conservative(s,V,b.h);err4=float(max(abs(dm.sum()),abs(dp.sum(axis=(1,2,3))).max(),abs((s*adv+V*dm-dp)).max()));sm=gaussian_filter(s,2,mode='constant');sv=np.array([gaussian_filter(s*v,2,mode='constant') for v in V]);err5=float(abs(integral(s,(sv-V*sm)/sm.max(),b.h)).max())
    assert max(err1,err2,err3,err4,err5)<1e-10,(err1,err2,err3,err4,err5)
    dump('validation',dict(chem_max=err1,wave_max=err2,affine_material_max=err3,flux_balance_max=err4,cohesion_balance_max=err5));print('validation passed',err1,err2,err3,err4,err5,flush=True)

def main():
    p=argparse.ArgumentParser();p.add_argument('--workers',type=int,default=4);p.add_argument('--phase',choices=['validate','settle','run'],required=True);p.add_argument('--cases',nargs='*');a=p.parse_args();tmp=load()
    if a.phase=='validate':validate();return
    cases=[c for c in CASES if not a.cases or c[0] in a.cases];keys=sorted(set((c[1],c[2]) for c in cases))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        futures=[pool.submit(start,*k) for k in keys] if a.phase=='settle' else [pool.submit(run,c) for c in cases]
        for fut in as_completed(futures):fut.result()

if __name__=='__main__':main()
