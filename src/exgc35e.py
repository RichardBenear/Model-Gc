"""GC35E: affine constraint projection with matched free and continuum controls."""
import argparse
import json
from pathlib import Path
import time
import exgc35c as C
import numpy as np

A, P = C.A, C.P
READOUTS=(4,8,12,16,20)


class Projected(C.Experiment):
    def __init__(self,n,lattice):
        super().__init__(n,lattice,True)
        # The supported nonlinear products are evaluated on 18 quadrature
        # points per axis with the physical lattice derivative and translations.
        self.n=18
        self.grid=C.R.StandingPacket(self.n,"generic",lattice)
        self.grid.h=self.h
        for x in self.grid.KK:
            self.grid.a[abs(x*P.LBOX/(2*np.pi))>2.01]=0
        self.ids=tuple((self.modes%self.n).T)
        self.symbols=[2*np.sin(v*self.h/2)/self.h if lattice else v for v in self.grid.KK]
        self.grid.w=np.sqrt(1+sum(v*v for v in self.symbols))
        self.translation=[np.exp(1j*v*self.h) for v in self.grid.KK]
        self.use_translation=lattice

    def evolve(self,dt,end):
        S=self.system
        k=np.linalg.norm(self.K,axis=1)
        weight=np.sqrt(np.tile(P.WT,2))
        matrix=np.concatenate([S.RC,S.RB*k[:,None,None]],axis=2)
        pinv=np.linalg.pinv(matrix/weight)/weight[:,None]
        unit=self.K/k[:,None]
        plane=np.eye(3)-unit[:,:,None]*unit[:,None,:]

        def tt(values):
            tensor=np.empty((len(values),3,3),dtype=values.dtype)
            for j,(a,b) in enumerate(P.N.PAIRS):
                tensor[:,a,b]=tensor[:,b,a]=values[:,j]
            out=plane@tensor@plane
            return out-np.trace(out,axis1=-2,axis2=-1)[:,None,None]*plane/2

        def norm(g):
            return float(np.sqrt(np.sum(P.WT*np.abs(g)**2)))

        def ttnorm(g):
            return float(np.linalg.norm(tt(g)))

        current=self.flat(0)
        m,mv,_=S.init(current[3])
        z=np.zeros_like(current[0])
        state=(m,mv,z,z,np.zeros_like(m),np.zeros_like(mv),np.zeros((*z.shape,4)),np.zeros_like(m),np.zeros_like(mv))
        initial=self.assemble(state,current)
        m2,m2v,_=S.init(initial[1])
        state=(*state[:4],m2,m2v,*state[6:])
        free=(m2.copy(),m2v.copy())
        readouts={}
        size=int(round(end/dt))
        window={int(round(t/dt)):t for t in READOUTS if t<=end}
        # Four samples give a second-order backward estimate at the final time.
        needed={j for i in window for j in (i-3,i-2,i-1,i,i+1) if 0<=j<=size}
        positions={}
        rows=[]
        totals=np.zeros(2)
        for i in range(size+1):
            t=i*dt
            first,f2,_=self.assemble(state,current)
            old,acc_free=S.cres(*free,f2)
            residual,_=S.cres(state[4],state[5],f2)
            correction=-P.mv(pinv,residual)
            dh,dhv=correction[:,:6],correction[:,6:]*k[:,None]
            state=(*state[:4],state[4]+dh,state[5]+dhv,*state[6:])
            after,acc=S.cres(state[4],state[5],f2)
            scale=np.linalg.norm(current[3][:,P.IC])
            totals+=norm(dh),norm(dhv)
            row=dict(t=t,before=float(A.EPSILON*np.linalg.norm(residual)/scale),
                     after=float(A.EPSILON*np.linalg.norm(after)/scale),
                     free=float(A.EPSILON*np.linalg.norm(old)/scale),
                     dh=norm(dh),dv=norm(dhv),tt_dh=ttnorm(dh),tt_dv=ttnorm(dhv))
            rows.append(row)
            if i in needed:
                positions[i]=(state[4].copy(),free[0].copy())
            if i in window:
                readouts[window[i]]=(-acc/2,-acc_free/2,-dhv/(2*dt),tt(-acc/2))
                print(row,flush=True)
            if i==size:
                break
            # The source hierarchy is triangular. Changing the second metric
            # changes only its own velocity and acceleration entries.
            def with_current_projection(s):
                return (*s[:4],state[5],acc,*s[6:])
            first=with_current_projection(first)
            middle,target=self.flat(t+dt/2),self.flat(t+dt)
            second,f2mid,_=self.assemble(A.plus(state,first,dt/2),middle)
            third,f2mid2,_=self.assemble(A.plus(state,second,dt/2),middle)
            fourth,f2end,_=self.assemble(A.plus(state,third,dt),target)
            state=tuple(x+dt/6*(a+2*b+2*c+d) for x,a,b,c,d in zip(state,first,second,third,fourth))
            def rhs(q,f):
                return q[1], S.acc(*q,f[:,P.IH])
            a=rhs(free,f2)
            b=rhs(A.plus(free,a,dt/2),f2mid)
            c=rhs(A.plus(free,b,dt/2),f2mid2)
            d=rhs(A.plus(free,c,dt),f2end)
            free=tuple(x+dt/6*(w+2*y+2*z+v) for x,w,y,z,v in zip(free,a,b,c,d))
            current=target
        arrays={}
        for i,t in window.items():
            def reconstructed(j):
                if i<size:
                    return -(positions[i+1][j]-2*positions[i][j]+positions[i-1][j])/(2*dt**2)
                return -(2*positions[i][j]-5*positions[i-1][j]+4*positions[i-2][j]-positions[i-3][j])/(2*dt**2)
            tidal,free_tidal,impulse,tt_tidal=readouts[t]
            acc=reconstructed(0)
            arrays.update({f"{t}_nominal":tidal,f"{t}_free":free_tidal,
                           f"{t}_impulse":impulse,f"{t}_reconstructed":acc,f"{t}_reconstructed_free":reconstructed(1),
                           f"{t}_tt":tt_tidal,f"{t}_reconstructed_tt":tt(acc)})
        return dict(history=rows,cumulative=totals.tolist()),arrays


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("tag")
    parser.add_argument("--n",type=int,default=24)
    parser.add_argument("--dt",type=float,default=.05)
    parser.add_argument("--end",type=int,default=20)
    parser.add_argument("--continuum",action="store_true")
    parser.add_argument("--output",default="EX-GC35E")
    args=parser.parse_args()
    start=time.process_time()
    run,arrays=Projected(args.n,not args.continuum).run(args.dt,args.end)
    run.update(settings=vars(args),cpu_seconds=time.process_time()-start)
    path=Path(__file__).resolve().parents[1]/"results"/args.output
    path.mkdir(parents=True,exist_ok=True)
    (path/f"{args.tag}.json").write_text(json.dumps(run,separators=(",",":"))+"\n")
    np.savez_compressed(path/f"{args.tag}.npz",**arrays)


if __name__=="__main__":
    main()
