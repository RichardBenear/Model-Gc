"""Mirror-symmetric Skyrme calibration; consistent static/dynamic stencil."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'):os.environ[k]='1'
import argparse,ctypes,hashlib,itertools,json,platform,subprocess,time,resource
from pathlib import Path
import numpy as np
import scipy
import exgc26a_skyrmion as ref
import exgc26a2_qualify as a2
import exgc26b_mass as b

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC26B2';RAW=ROOT/'local_raw/EX-GC26B2';NORM=ref.NORM;PTR=a2.PTR;LIB=None;POWER=1
def save(tag,r):
    OUT.mkdir(parents=True,exist_ok=True);p=OUT/(tag+'.json');q=p.with_suffix('.tmp');q.write_text(json.dumps(r,indent=1,allow_nan=False)+'\n');q.replace(p)
def load():
    global LIB
    RAW.mkdir(parents=True,exist_ok=True);p=RAW/'stencil8.so';src=Path(__file__).with_name('exgc26b2_stencil.cpp')
    if not p.exists() or p.stat().st_mtime<src.stat().st_mtime:
        q=p.with_name(f'stencil8-{os.getpid()}.so');subprocess.run(['c++','-O3','-march=native','-shared','-fPIC','-std=c++17',str(src),'-o',str(q)],check=True);q.replace(p)
    LIB=ctypes.CDLL(str(p));LIB.evaluate8.argtypes=[ctypes.c_int,ctypes.c_double,PTR,PTR,PTR];LIB.acceleration8.argtypes=LIB.acceleration8_coreforce.argtypes=[ctypes.c_int,ctypes.c_double,PTR,PTR,ctypes.c_double,PTR];LIB.evaluate8.restype=LIB.acceleration8.restype=LIB.acceleration8_coreforce.restype=None
def evaluate(G,p,mode=0):
    x=np.ascontiguousarray(np.moveaxis(p,0,-1));y=np.empty_like(x);e=np.zeros(2);LIB.evaluate8(G.n,G.h,x,y,e);return e,np.moveaxis(y,-1,0)
def accel(G,p,v,f=0.):
    x=np.ascontiguousarray(np.moveaxis(p,0,-1));y=np.ascontiguousarray(np.moveaxis(v,0,-1));z=np.empty_like(x);(LIB.acceleration8 if POWER==1 else LIB.acceleration8_coreforce)(G.n,G.h,x,y,f,z);return np.moveaxis(z,-1,0)
def orient(G,p,s):
    mask=np.zeros_like(G.inner);mask[tuple(slice(None,-1) if x==1 else slice(1,None) for x in s)]=True
    d=[(np.roll(p,-s[i],axis=i+1)-p)*s[i]/G.h*mask for i in range(3)];D=np.array([[b.dot(x,y) for y in d] for x in d]);return mask,d,D
def adj(G,q,i,s,mask):return s*(np.roll(q*mask,s,axis=i+1)-q*mask)/G.h
def numpy_check(G,p,v,f):
    """Independent direct NumPy Euler-Lagrange sum; validation only."""
    E=np.zeros(2);grad=np.zeros_like(p);R=np.zeros_like(p);M=np.zeros((4,4)+p.shape[1:]);vv=b.dot(v,v)
    for s in itertools.product((-1,1),repeat=3):
        m,d,D=orient(G,p,s);_,dv,_=orient(G,v,s);tr=sum(D[i,i] for i in range(3));e4=.5*(tr*tr-np.sum(D*D,axis=(0,1)));E+=np.array([tr.sum(),e4.sum()])/8;vd=[b.dot(v,x) for x in d];sd=sum(b.dot(x,y) for x,y in zip(d,dv));M+=((1+tr)*np.eye(4)[:,:,None,None,None]-sum(x[:,None]*x[None,:] for x in d))/8
        R-=(2*sd*v-sum(dv[i]*vd[i]+d[i]*b.dot(dv[i],v) for i in range(3)))/8
        for i in range(3):
            q=(1+tr)*d[i]-sum(D[i,j]*d[j] for j in range(3));grad+=adj(G,2*q,i,s[i],m)/8;R+=adj(G,vv*d[i]-vd[i]*v-q,i,s[i],m)/8
    R[0]-=.5*f*G.X[0]*POWER*(1-p[0])**(POWER-1);A=np.moveaxis(M,(0,1),(-2,-1));rhs=np.stack([np.moveaxis(R,0,-1),np.moveaxis(p,0,-1)],-1);sol=np.linalg.solve(A,rhs);a=np.moveaxis(sol[...,0],-1,0);c=np.moveaxis(sol[...,1],-1,0);return E*G.h**3/NORM,grad*G.h**3/NORM,(a+c*(-vv-b.dot(p,a))/b.dot(p,c))*G.inner
def fields(G,p,v):
    H=np.zeros(p.shape[1:]);P=np.zeros((3,)+H.shape);S=np.zeros((3,3)+H.shape);vv=b.dot(v,v)
    for s in itertools.product((-1,1),repeat=3):
        mask,d,D=orient(G,p,s);tr=sum(D[i,i] for i in range(3));vd=np.array([b.dot(v,x) for x in d]);K=vv*(1+tr)-np.sum(vd**2,axis=0);e=tr+.5*(tr*tr-np.sum(D**2,axis=(0,1)));H+=(K+e)*mask/8
        for i in range(3):
            P[i]-=2*((1+tr)*vd[i]-sum(D[i,j]*vd[j] for j in range(3)))*mask/8
            for j in range(3):S[i,j]+=(2*((1+tr-vv)*D[i,j]-sum(D[i,k]*D[j,k] for k in range(3))+vd[i]*vd[j])+(K-e)*(i==j))*mask/8
    return H,P,S
def validate():
    t=time.process_time();G=ref.Grid(.25,1);rng=np.random.default_rng(262);p=G.fix(G.vacuum()+.17*rng.normal(size=(4,G.n,G.n,G.n)));v=ref.tangent(p,rng.normal(scale=.05,size=p.shape))*G.inner;f=.07;e,g=evaluate(G,p);er,gr,ar=numpy_check(G,p,v,f);a=accel(G,p,v,f);q=ref.tangent(p,rng.normal(size=p.shape))*G.inner;q/=np.linalg.norm(q);eps=1e-5;fd=(evaluate(G,G.fix(p+eps*q))[0].sum()-evaluate(G,G.fix(p-eps*q))[0].sum())/(2*eps);de=abs(fd-np.sum(g*q))/max(abs(fd),1e-12)
    pr=p[:,:,::-1,:].copy();vr=v[:,:,::-1,:].copy();a0=accel(G,p,v,0);ar0=accel(G,pr,vr,0);H,P,_=fields(G,p,v);Hr,Pr,_=fields(G,pr,vr);fac=G.h**3/NORM;ch=G.vacuum();ch[1]=.2*(-1.)**np.indices(ch.shape[1:]).sum(0);ch=G.fix(ch);vac=G.vacuum();z=np.zeros_like(vac)
    # Exact energy time derivative with external term, estimated symmetrically using p_t=v and v_t=a.
    def total(p,v):return fields(G,p,v)[0].sum()*fac-f*np.sum(G.X[0]*(1-p[0])**POWER)*fac
    ep=1e-6;rate=(total(G.fix(p+ep*v),ref.tangent(G.fix(p+ep*v),v+ep*a))-total(G.fix(p-ep*v),ref.tangent(G.fix(p-ep*v),v-ep*a)))/(2*ep)
    r=dict(E_error=float(abs(e-er).max()),gradient_error=float(abs(g-gr).max()),acceleration_error=float(abs(a-ar).max()),directional_gradient_rel=float(de),mirror_accel_error=float(abs(ar0-a0[:,:,::-1,:]).max()),mirror_H_error=float(abs(H.sum()-Hr.sum())*fac),mirror_P_error=float(np.max(abs(P.sum((1,2,3))*np.array([1,-1,1])-Pr.sum((1,2,3))))*fac),checkerboard_E=float(evaluate(G,ch)[0].sum()),vacuum_E=float(evaluate(G,vac)[0].sum()),vacuum_accel=float(abs(accel(G,vac,z,f)).max()),forced_energy_rate=float(rate),python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform(),cpu_s=time.process_time()-t);save('validation' if POWER==1 else 'validation_coreforce',r);print(r,flush=True);assert max(r[k] for k in ('E_error','gradient_error','acceleration_error','mirror_accel_error','mirror_H_error','mirror_P_error','vacuum_E','vacuum_accel'))<1e-10 and de<1e-5 and abs(rate)<1e-7 and r['checkerboard_E']>0
def restpath(h,L):return RAW/f'rest_h{h:g}_L{L:g}.npy'
def prepare(h,L,maxiter):
    t=time.process_time();G=ref.Grid(h,L);a2.evaluate=evaluate;p,r=a2.optimize(G,G.hedgehog(ref.PROF),0,maxiter);np.save(restpath(h,L),p);r.update(h=h,L=L,actual_L=float(G.ax[-1]),degree=a2.geometric(p),sha256=hashlib.sha256(restpath(h,L).read_bytes()).hexdigest(),cpu_s=time.process_time()-t);save(f'static_h{h:g}_L{L:g}',r);print('rest',h,L,r,flush=True)
def step(G,p,v,dt,f):
    k1v=accel(G,p,v,f);k1p=v;k2p=v+dt*k1v/2;k2v=accel(G,G.fix(p+dt*k1p/2),k2p,f);k3p=v+dt*k2v/2;k3v=accel(G,G.fix(p+dt*k2p/2),k3p,f);k4p=v+dt*k3v;k4v=accel(G,G.fix(p+dt*k3p),k4p,f);q=G.fix(p+dt*(k1p+2*k2p+2*k3p+k4p)/6);return q,ref.tangent(q,v+dt*(k1v+2*k2v+2*k3v+k4v)/6)*G.inner
def measure(G,p,v,t,f):
    H,P,S=fields(G,p,v);fac=G.h**3/NORM;Q=float(np.sum((1-p[0])**POWER)*fac);V=float(-f*np.sum(G.X[0]*(1-p[0])**POWER)*fac);return dict(t=t,H=float(H.sum()*fac),P=(P.sum((1,2,3))*fac).tolist(),centroid=[float(np.sum(H*x)/H.sum()) for x in G.X],regions=[b.region(G,H,P,S,R) for R in (3,5)],Q=Q,V=V,F=f*Q)
def fit(t,y):return float(np.polyfit(t,y,1)[0])
def run(h,L,u,axis,T,sign,dtfac):
    start=time.process_time();G=ref.Grid(h,L);p=np.load(restpath(h,L));E0=float(evaluate(G,p)[0].sum());fac=G.h**3/NORM;Q0=float(np.sum((1-p[0])**POWER)*fac);f=sign*.025*E0/Q0;n=np.array(axis,float);n/=np.linalg.norm(n);p,v=b.init_boost(G,p,u,n) if u else (p,np.zeros_like(p));d0=a2.geometric(p);dt=dtfac*h;steps=round(T/dt);every=max(1,round(.1/dt));rows=[];tag=f'h{h:g}_L{L:g}_u{u:g}_axis{"".join(map(str,axis))}_T{T:g}_force{sign:+d}_dt{dtfac:g}'+('_coreforce' if POWER==2 else '')
    for i in range(steps+1):
        if i%every==0 or i==steps:
            rows.append(measure(G,p,v,i*dt,f))
            if i%max(every,round(1/dt))==0:print(tag,'t',i*dt,'cpu',time.process_time()-start,flush=True)
        if i<steps:p,v=step(G,p,v,dt,f)
    t=np.array([r['t'] for r in rows]);C=np.array([r['centroid'] for r in rows])@n;P=np.array([r['P'] for r in rows]);H=np.array([r['H'] for r in rows]);V=np.array([r['V'] for r in rows]);F=np.array([r['F'] for r in rows]);s=b.summary(rows,E0,n,u);s['total_energy_drift']=float(max(abs((H+V)/(H[0]+V[0])-1)));s['displacement']=float(C[-1]-C[0]);s['speed_over_P_H']=float(fit(t,C)/np.mean(P@n/H)) if u else None
    if u:
        w=min(1.,t[-1]/2);a=t<=w+1e-12;z=t>=t[-1]-w-1e-12;s['early_speed']=fit(t[a],C[a]);s['late_speed']=fit(t[z],C[z]);s['relative_speed_change']=s['late_speed']/s['early_speed']-1
    if sign:
        coef=np.polyfit(t,C,2);a=2*coef[0];imp=np.r_[0.,np.cumsum((F[:-1]+F[1:])*np.diff(t)/2)];s.update(acceleration=float(a),mass_over_E0_force=float(np.mean(F)/(a*E0)),P_slope_over_mean_force=float(fit(t,P@n)/np.mean(F)),impulse_rel_error=float(max(abs((P-P[0])@n-imp))/abs(imp[-1])),end_impulse=float(imp[-1]),fit_residual_max=float(max(abs(C-np.polyval(coef,t)))),free_fit_initial_velocity=float(coef[1]))
    dst=RAW/(tag+'_end.npy');np.save(dst,p);r=dict(h=h,L=L,u=u,axis=axis,T=T,dt=dt,f=f,force_sign=sign,force_power=POWER,E0=E0,Q0=Q0,rows=rows,summary=s,degree_initial=d0,degree_final=a2.geometric(p),end_sha256=hashlib.sha256(dst.read_bytes()).hexdigest(),cpu_s=time.process_time()-start);save(tag,r);print('DONE',tag,s,r['cpu_s'],flush=True)
def main():
    global POWER
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['validate','static','run']);ap.add_argument('--h',type=float,default=.2);ap.add_argument('--L',type=float,default=7);ap.add_argument('--maxiter',type=int,default=300);ap.add_argument('--u',type=float,default=.3);ap.add_argument('--axis',type=int,nargs=3,default=[1,0,0]);ap.add_argument('--T',type=float,default=.8);ap.add_argument('--sign',type=int,default=0,choices=[-1,0,1]);ap.add_argument('--power',type=int,choices=[1,2],default=1);ap.add_argument('--dtfac',type=float,default=.05);ap.add_argument('--cpu-limit',type=int,default=1800);a=ap.parse_args();POWER=a.power;resource.setrlimit(resource.RLIMIT_CPU,(a.cpu_limit,a.cpu_limit+1));load()
    if a.phase=='validate':validate()
    elif a.phase=='static':prepare(a.h,a.L,a.maxiter)
    else:run(a.h,a.L,a.u,a.axis,a.T,a.sign,a.dtfac)
if __name__=='__main__':main()
