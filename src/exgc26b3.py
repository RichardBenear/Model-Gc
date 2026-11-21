"""Smoothly ramped core-localized force; inherited B2 spatial discretization."""
import argparse, hashlib, json, resource, time
from pathlib import Path
import exgc26b2 as b
import numpy as np

OUT=b.ROOT/'results/EX-GC26B3'
RAW=b.ROOT/'local_raw/EX-GC26B3'

def ramp(t):
    return .5*(1-np.cos(np.pi*np.clip(t,0,1)))

def step(G,p,v,t,dt,f):
    a=b.accel
    k1p=v;k1v=a(G,p,v,f*ramp(t))
    k2p=v+dt*k1v/2;k2v=a(G,G.fix(p+dt*k1p/2),k2p,f*ramp(t+dt/2))
    k3p=v+dt*k2v/2;k3v=a(G,G.fix(p+dt*k2p/2),k3p,f*ramp(t+dt/2))
    k4p=v+dt*k3v;k4v=a(G,G.fix(p+dt*k3p),k4p,f*ramp(t+dt))
    q=G.fix(p+dt*(k1p+2*k2p+2*k3p+k4p)/6)
    return q,b.ref.tangent(q,v+dt*(k1v+2*k2v+2*k3v+k4v)/6)*G.inner

def run(h,sign):
    start=time.process_time(); b.POWER=2;b.load();G=b.ref.Grid(h,7);path=b.restpath(h,7)
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    assert digest==json.loads((b.OUT/f'static_h{h:g}_L7.json').read_text())['sha256']
    p=np.load(path);v=np.zeros_like(p);E0=float(b.evaluate(G,p)[0].sum());fac=h**3/b.NORM
    Q0=float(np.sum((1-p[0])**2)*fac);f=sign*.025*E0/Q0;dt=.05*h;steps=round(5/dt);every=max(1,round(.1/dt));rows=[];degree=b.a2.geometric(p)
    tag=f'h{h:g}_force{sign}';OUT.mkdir(parents=True,exist_ok=True);RAW.mkdir(parents=True,exist_ok=True)
    for i in range(steps+1):
        t=i*dt
        if i%every==0 or i==steps:
            r=b.measure(G,p,v,t,f*ramp(t));r['D']=float(np.sum(G.X[0]*(1-p[0])**2)*fac);rows.append(r)
            if i%every==0 and (i//every)%10==0:print(tag,t,time.process_time()-start,flush=True)
        if i<steps:p,v=step(G,p,v,t,dt,f)
    end=RAW/(tag+'_end.npy');np.save(end,p)
    result=dict(h=h,L=7,force=sign,dt=dt,T=steps*dt,E0=E0,Q0=Q0,f=f,rows=rows,degree_initial=degree,degree_final=b.a2.geometric(p),rest_sha256=digest,end_sha256=hashlib.sha256(end.read_bytes()).hexdigest(),cpu_s=time.process_time()-start)
    (OUT/(tag+'.json')).write_text(json.dumps(result,indent=1,allow_nan=False)+'\n');print('DONE',tag,result['cpu_s'],flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--h',type=float,choices=[.2,.15],required=True);ap.add_argument('--force',type=int,choices=[0,1],required=True);a=ap.parse_args();cap=1050 if a.h==.2 else 2450;resource.setrlimit(resource.RLIMIT_CPU,(cap,cap+1));run(a.h,a.force)
