"""Static-source diagnostic for GC35M, holding the knots fixed."""
from itertools import product
import argparse
import json
from pathlib import Path
import time
import exgc35h as H
import exgc35m as M
import numpy as np

P=M.P


def forcing(p,h):
    x,y,s,u,_,_=H.components(p,h)
    n=len(p)
    k=2*np.pi*np.array(list(product(range(-5,6),repeat=3)))/(n*h)
    k=k[(np.linalg.norm(k,axis=1)<=2.05)&np.any(k!=0,axis=1)]
    ids=tuple((np.rint(k*n*h/(2*np.pi)).astype(int)%n).T)
    result=np.zeros((len(k),4,4),complex)
    result[:,0,0]=P.sf.fftn(x+y)[ids]/n**3
    for j,(i,b) in enumerate(H.PAIRS):
        phase=np.exp(-1j*k@(h*np.asarray(P.N.OFFS[j+4])))
        c=P.sf.fftn((s+u)[...,j])[ids]*phase/n**3
        result[:,i+1,b+1]=result[:,b+1,i+1]=c
    return 2*np.sin(k*h/2)/h,result.reshape(-1,16)@P.N.PM*P.KAP


def evolve(K,f):
    S=P.System(K)
    adj=S.RB.conj().swapaxes(1,2)/P.WT[None,:,None]
    h,v,_=S.init(f)
    scale=np.linalg.norm(f[:,:4])
    vals=[]
    dt=.025

    def rhs(x):
        a,b=x
        r,acc=S.cres(a,b,f)
        return b+M.ALPHA*r[:,:1]*M.TRACE,acc-M.BETA*P.mv(adj,r)

    x=(h,v)
    for i in range(round(20/dt)+1):
        a=rhs(x)
        physical,_=S.cres(x[0],a[0],f)
        vals.append(float(np.linalg.norm(physical)/scale))
        b=rhs(M.A.plus(x,a,dt/2))
        c=rhs(M.A.plus(x,b,dt/2))
        d=rhs(M.A.plus(x,c,dt))
        x=tuple(u+dt*(ka+2*kb+2*kc+kd)/6 for u,ka,kb,kc,kd in zip(x,a,b,c,d))
    return dict(peak=max(vals),final=vals[-1],at5=vals[200],at10=vals[400])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("degree",type=int)
    args=parser.parse_args()
    start=time.process_time()
    k=H.D.Knot(.2,args.degree,6)
    p,relax=k.relax(450)
    K,f=forcing(p.numpy(),.2)
    g=f.copy()
    g[:,4:]=0
    result=dict(degree=args.degree,force=relax["force"],full=evolve(K,f),energy=evolve(K,g),
                cpu_seconds=time.process_time()-start)
    root=Path(__file__).resolve().parents[1]/"results"/"EX-GC35M"
    root.mkdir(parents=True,exist_ok=True)
    (root/f"static_b{args.degree}.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)


if __name__=="__main__":
    main()
