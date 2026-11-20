"""EX-GC22: steady prescribed sink and full Gc flow, independent of EX-GC20 initialization.

Our construction. Run from repo root: python src/exgc22_drag.py fixed; python
src/exgc22_drag.py full --workers 2. All production settings in spec/EX-GC22_plan.md.
"""
import argparse
import json
import os
import platform
import ctypes
import subprocess
import tempfile
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import scipy
from scipy.fft import dstn, idstn, rfftn, irfftn
from scipy.ndimage import gaussian_filter
import ex008_longrange as e
import exgc09_lattice as g9
import exgc20_wave as wave

OUT = Path(e.ROOT) / 'results/EX-GC22'
KERNEL=None


def load_kernel():
    global KERNEL
    tmp=tempfile.TemporaryDirectory(prefix='exgc22-kernel-')
    path=Path(tmp.name)/'chem.so'
    subprocess.run(['c++','-O3','-shared','-fPIC',str(Path(__file__).with_name('exgc22_kernel.cpp')),'-o',str(path)],check=True)
    lib=ctypes.CDLL(str(path)); KERNEL=lib.chem
    ptr=np.ctypeslib.ndpointer(dtype=np.float64,flags='C_CONTIGUOUS')
    KERNEL.argtypes=[ctypes.c_int]+[ctypes.c_double]*7+[ptr,ptr,ptr]
    KERNEL.restype=None
    return tmp


def dump(name, obj):
    OUT.mkdir(exist_ok=True, parents=True)
    (OUT / name).write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n')


def fixed(K, v, tj, c=2., side=32., h=.5):
    n = round(side/h); z = (np.arange(n)-n//2)*h
    r2 = z[:, None, None]**2+z[None, :, None]**2+z[None, None, :]**2
    s = np.exp(-r2/(2*.8**2)); s *= 5.32*K/.2/(s.sum()*h**3)
    q = 2*np.pi*np.fft.fftfreq(n, h); qz=2*np.pi*np.fft.rfftfreq(n, h)
    q2=q[:,None,None]**2+q[None,:,None]**2+qz[None,None,:]**2
    w=v*qz[None,None,:]; sh=rfftn(s); sh[0,0,0]=0
    den=c*c*q2-w*w-1j*w/tj; den[0,0,0]=1
    gh=-(1/tj-1j*w)*sh/den
    jh=1j*c*c*qz[None,None,:]*sh/den
    g=1+irfftn(gh,s=s.shape); j=irfftn(jh,s=s.shape)
    a=2*(s*j/g).sum()/s.sum()/20
    al=2*(s*j).sum()/s.sum()/20
    # Compute residual directly from Fourier equations, and the first-order v coefficient.
    # J = i c² k s / den; div J = -c² |k|² s / den.
    res1=np.max(np.abs(-1j*w*gh-c*c*q2*sh/den+sh))
    res2=np.max(np.abs((1/tj-1j*w)*jh+1j*c*c*qz[None,None,:]*gh))
    d=c*c*tj
    linear_jh=-qz[None,None,:]**2*sh/(d*np.where(q2>0,q2,1)**2)
    linear_j=irfftn(linear_jh,s=s.shape)
    lin_gamma=-2*(s*linear_j).sum()/s.sum()/20
    return dict(K=K,v=v,tauJ=tj,c=c,D=d,side=side,h=h,M=float(s.sum()*h**3),
                depletion=float(1-g[n//2,n//2,n//2]),Gmin=float(g.min()),a=float(a),
                gamma=float(-a/v) if v else None, gamma_G1=float(-al/v) if v else None,
                gamma_linear=float(lin_gamma),external_a=float(2*.02*(s/g).sum()/s.sum()/20),
                residual=float(max(res1,res2)/max(1.,abs(sh).max())))


def fixed_main():
    t=time.monotonic(); cpu=time.process_time()
    rows=[fixed(k,v,t) for k in (.2,.05,.01,.001) for t in (.25,2,20) for v in (.01,.025,.05)]
    controls=[fixed(.2,v,t) for t in (.25,2,20) for v in (0,-.05)]
    controls += [fixed(.2,.05,t,side=64) for t in (.25,20)]
    controls += [fixed(.2,.05,t,h=.25) for t in (.25,20)]
    dump('fixed.json',dict(rows=rows,controls=controls,wall_s=time.monotonic()-t,cpu_s=time.process_time()-cpu,
                          platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__))
    for r in rows:
        if r['v']==.05: print('fixed',r['K'],r['D'],r['depletion'],r['gamma'],r['external_a'],flush=True)


class Box:
    def __init__(self,h,L,D,K):
        self.h=h; self.D=D; self.K=K; self.ax,self.R=g9.grid(h,L); self.n=len(self.ax)
        self.shape=self.R.shape
        z=self.ax; self.Z=np.broadcast_to(z,self.shape)
        q=4*np.sin(np.pi*np.arange(1,self.n+1)/(2*(self.n+1)))**2/h**2
        self.inv=1/(q[:,None,None]+q[None,:,None]+q[None,None,:])
        self.f0=(e.A,e.B/e.A,0.)
        self.dt=.0015*(h/.5)**2; self.kf=round(.0075/self.dt)
        self.dtV=self.dt*self.kf
    def poisson(self,W,G):
        # Dirichlet ghost box G=1, exact stationary discrete G for fixed W.
        for _ in range(100):
            new=1-idstn(dstn(self.K*W*G,type=1,norm='ortho')*self.inv/self.D,type=1,norm='ortho')
            err=np.max(abs(new-G)); G=new
            if err<1e-10: return G
        raise RuntimeError('Poisson iteration failed')
    def bc(self,E):
        z=self.ax; Z=self.Z
        return [1+E*Z[0],1+E*Z[-1],1+E*Z[:,0],1+E*Z[:,-1],
                np.full(Z[:,:,0].shape,1+E*(z[0]-self.h)),np.full(Z[:,:,-1].shape,1+E*(z[-1]+self.h))]
    def vel(self,u,V,push):
        s=self.K*u[3]*u[2]; wgt=s/(s+1e-3*s.max()); sig=1.5/self.h
        sm=gaussian_filter(s,sig,mode='constant')+1e-12
        Vbar=np.array([gaussian_filter(s*vv,sig,mode='constant')/sm for vv in V])
        return V+self.dtV*(wgt*np.array(push)/20-(V-Vbar))
    def chem(self,u,V):
        if KERNEL is not None:
            out=np.empty((3,*self.shape)); KERNEL(self.n,self.h,self.K,e.A,e.B,e.DY,e.DW,e.C,u,V,out); return out
        return self.chem_numpy(u,V)
    def chem_numpy(self,u,V):
        X,Y,W,G=u; s=self.K*G*W; adv=[-wave.div([(u[j]-self.f0[j])*vv for vv in V],self.h) for j in range(3)]
        return np.array([g9.lap(X,self.h)+e.A-(1+e.B)*X+X*X*Y+s+adv[0],
                         e.DY*g9.lap(Y,self.h)+e.B*X-X*X*Y+adv[1],
                         e.DW*g9.lap(W,self.h)-W+e.C*Y*W*W-W**3+adv[2]])
    def measure(self,t,u,V):
        s=self.K*u[3]*u[2]; ss=s.sum(); W=u[2]; wc=(W*self.Z).sum()/W.sum()
        vc=np.array([(s*vv).sum()/ss for vv in V]); spread=np.sqrt((s*((V-vc[:,None,None,None])**2).sum(0)).sum()/ss)
        return dict(t=t,M=float(ss*self.h**3),Wmax=float(W.max()),Wcen_z=float(wc),Vcore_z=float(vc[2]),
                    Vspread=float(spread),depletion=float(1-u[3].min()),Wint=float(W.sum()*self.h**3))
    def settle(self,T=30):
        _,_,u=g9.init(self.h,float(self.ax[-1])); V=np.zeros((3,*self.shape)); u[3]=self.poisson(u[2],u[3])
        rows=[]; n=round(T/self.dt); every=round(1/self.dt)
        for i in range(n+1):
            if i%(self.kf*10)==0 or i==n:
                u[3]=self.poisson(u[2],u[3])
            if i%self.kf==0 or i==n:
                # Use identical face-centred force as production.
                push=[2*j/u[3] for j in wave.centre_flux([-self.D*g for g in wave.grad_faces(u[3],self.h,self.bc(0))])]
                V=self.vel(u,V,push)
            if i%every==0 or i==n: rows.append(self.measure(i*self.dt,u,V))
            if i<n: u[:3]+=self.dt*self.chem(u,V)
        return u,V,rows
    def evolve(self,u0,V0,kick=0.,flux=0.,T=20):
        u=u0.copy(); V=V0.copy(); V[2]+=kick; E=-flux/self.D; u[3]+=E*self.Z; bc=self.bc(E)
        J=[-self.D*g for g in wave.grad_faces(u[3],self.h,bc)] if self.D>1 else None
        tj=self.D/4
        rows=[]; n=round(T/self.dt); every=round(.5/self.dt)
        for i in range(n+1):
            if i%every==0 or i==n: rows.append(self.measure(i*self.dt,u,V))
            if i==n: break
            if i%self.kf==0:
                push=[2*j/u[3] for j in wave.centre_flux(J if J is not None else [-g for g in wave.grad_faces(u[3],self.h,bc)])]
                V=self.vel(u,V,push)
            X,Y,W,G=u; s=self.K*G*W
            u[:3]+=self.dt*self.chem(u,V)
            if J is None: u[3]+=self.dt*(wave.div_faces(wave.grad_faces(G,self.h,bc),self.h)-s)
            else:
                gf=wave.grad_faces(G,self.h,bc)
                u[3]+=self.dt*(-wave.div_faces(J,self.h)-s)
                J=[jj+self.dt*(-4*gg-jj/tj) for jj,gg in zip(J,gf)]
            if not np.isfinite(u).all() or W.max()<.01: raise RuntimeError('particle lost')
        return rows


def full_case(args):
    K,D,h,flux_test=args; t=time.monotonic(); cpu=time.process_time(); b=Box(h,6,D,K)
    cache=Path(e.ROOT)/'local_raw'/'EX-GC22'; cache.mkdir(parents=True,exist_ok=True)
    cachefile=cache/f'settle_K{K:g}_D{D:g}_h{h:g}.npz'
    settlefile=OUT/f'settle_K{K:g}_D{D:g}_h{h:g}.json'
    if cachefile.exists() and settlefile.exists():
        with np.load(cachefile) as dat: u,V=dat['u'],dat['V']
        settle=json.loads(settlefile.read_text())['settle']
    else:
        u,V,settle=b.settle(); np.savez_compressed(cachefile,u=u,V=V)
    print('settled',K,D,h,settle[-1],flush=True)
    dump(f'settle_K{K:g}_D{D:g}_h{h:g}.json',dict(K=K,D=D,h=h,settle=settle))
    rows=b.evolve(u,V,kick=.05)
    out=dict(K=K,D=D,h=h,dt=b.dt,settle=settle,kick=rows)
    if flux_test:
        out['flux']=b.evolve(u,V,flux=.02)
        out['null']=b.evolve(u,V,T=10)
    for lo,hi,name in [(2,10,'early'),(10,20,'late')]:
        a=[r for r in rows if lo<=r['t']<=hi and r['Vcore_z']>0]
        out['gamma_'+name]=float(-np.polyfit([r['t'] for r in a],np.log([r['Vcore_z'] for r in a]),1)[0]) if len(a)>2 else None
    if flux_test:
        a=[r for r in out['flux'] if r['t']<=5]
        out['flux_accel_early']=float(np.polyfit([r['t'] for r in a],[r['Vcore_z'] for r in a],1)[0])
    out.update(wall_s=time.monotonic()-t,cpu_s=time.process_time()-cpu,platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__)
    dump(f'full_K{K:g}_D{D:g}_h{h:g}.json',out)
    print('full',K,D,h,settle[-1],out['gamma_early'],out['gamma_late'],out.get('flux_accel_early'),flush=True)
    return out


def full_main(workers):
    cases=[(K,D,.5,K in (.2,.001)) for D in (1,80) for K in (.2,.05,.01,.001)]
    with ProcessPoolExecutor(max_workers=workers) as pool:
        list(pool.map(full_case,cases))


if __name__=='__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('mode',choices=['fixed','full','refine','smoke']); parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    if args.mode=='fixed': fixed_main()
    else:
        tmp=load_kernel(); b=Box(.5,2,1,.2); rng=np.random.default_rng(22)
        u=rng.uniform(.5,2,(4,*b.shape)); V=rng.normal(0,.05,(3,*b.shape))
        assert np.max(np.abs(b.chem(u,V)-b.chem_numpy(u,V)))<1e-10
        print('C++ chemistry stencil matches NumPy smoke test',flush=True)
        if args.mode=='smoke':
            for D in (1,80):
                b=Box(.5,2,D,.01); u,V,_=b.settle(T=.03); b.evolve(u,V,kick=.05,T=.03)
            print('both integration branches pass short smoke test')
        elif args.mode=='full': full_main(args.workers)
        else: full_case((.2,80,.25,False))
