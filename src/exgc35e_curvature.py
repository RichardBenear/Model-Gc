"""Quadratic leading-metric contribution to the synchronous tidal tensor."""
import argparse
from pathlib import Path
import time
import exgc35e as E
import numpy as np

P,C=E.P,E.C


def correction(experiment,dt,end):
    x=experiment
    S=x.system
    readouts={round(t/dt):t for t in E.READOUTS if t<=end}
    result={}
    def f(t):
        return x.flat(t)[3]
    with x.stencil():
        source=f(0)
        m,v,_=S.init(source)
        y=m,v
        def rhs(state,src):
            return state[1],S.acc(*state,src[:,P.IH])
        for i in range(round(end/dt)+1):
            if i in readouts:
                vmat=C.A.tensor(x.metric(y[1]))[...,1:,1:]
                quadratic=vmat@vmat/4
                field=np.stack([quadratic[...,a,b] for a,b in P.N.PAIRS],axis=-1)
                result[str(readouts[i])]=P.sf.fftn(field,axes=(0,1,2))[x.ids]/x.n**3
            if i==round(end/dt):
                break
            mid,fin=f((i+.5)*dt),f((i+1)*dt)
            a=rhs(y,source)
            b=rhs(C.A.plus(y,a,dt/2),mid)
            c=rhs(C.A.plus(y,b,dt/2),mid)
            d=rhs(C.A.plus(y,c,dt),fin)
            y=tuple(z+dt*(ka+2*kb+2*kc+kd)/6 for z,ka,kb,kc,kd in zip(y,a,b,c,d))
            source=fin
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="EX-GC35E")
    args=parser.parse_args()
    start=time.process_time()
    path=Path(__file__).resolve().parents[1]/"results"/args.output
    path.mkdir(parents=True,exist_ok=True)
    for tag,n,lattice,dt,end in (("coarse",24,True,.05,20),("fine",48,True,.05,20),
                                 ("continuum",24,False,.05,20),("half",24,True,.025,4)):
        arr=correction(E.Projected(n,lattice),dt,end)
        np.savez_compressed(path/f"{tag}_quadratic.npz",**arr)
        print(tag,time.process_time()-start,flush=True)
    (path/"quadratic_cpu.txt").write_text(str(time.process_time()-start))
