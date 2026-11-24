"""Longer GC35A runs with band-limited quadrature and mode diagnostics."""
import argparse
from contextlib import contextmanager
import json
from pathlib import Path
import time

import exgc35a as A
import numpy as np

R, P = A.R, A.P


class Experiment(A.Experiment):
    def __init__(self, physical_n, lattice, compensate):
        super().__init__(24, lattice)
        self.h = P.LBOX / physical_n
        self.grid.h = self.h
        self.compensate = compensate
        self.modes = P.mask_modes()
        k, K = P.lattice_K(self.modes, self.h)
        self.K = K if lattice else k
        self.system = P.System(self.K)
        self.phase = np.exp(1j*k @ (self.h*self.offset).T) if lattice else np.ones((len(k),10))
        self.symbols = [2*np.sin(v*self.h/2)/self.h if lattice else v for v in self.grid.KK]
        self.grid.w = np.sqrt(1+sum(v*v for v in self.symbols))
        self.translation = [np.exp(1j*v*self.h) for v in self.grid.KK]
        self.use_translation = lattice and physical_n != self.n

    def shift(self, x, a, s=1):
        fft = P.sf.fftn(x)
        return P.sf.ifftn(fft * (self.translation[a] if s==1 else self.translation[a].conj())).real

    @contextmanager
    def stencil(self):
        if not self.use_translation:
            yield
            return
        # The original stencil uses these shifts on scalar grid values only.
        patches = [(P.O,"fw",self.shift),
                   (P.N,"grads",lambda f,h:[(self.shift(f,i)-f)/h for i in range(3)]),
                   (R,"forward",self.shift),
                   (R,"backward",lambda f,i:self.shift(f,i,-1)),
                   (R,"lap",lambda f,grid,lattice:P.sf.ifftn(-sum(k*k for k in self.symbols)*P.sf.fftn(f)).real)]
        original = [(module,key,getattr(module,key)) for module,key,_ in patches]
        try:
            for module,key,fn in patches:
                setattr(module,key,fn)
            yield
        finally:
            for module,key,fn in original:
                setattr(module,key,fn)

    def assemble(self, state, flat):
        dv, f, df, _ = super().assemble(state, flat)
        if not self.compensate:
            dv = (*dv[:6],np.zeros_like(dv[6]),*dv[7:])
        return dv, f, df

    def run(self, dt, end):
        with self.stencil():
            return self.evolve(dt,end)

    def evolve(self,dt,end):
        current=self.flat(0)
        m,mv,_=self.system.init(current[3])
        z=np.zeros_like(current[0])
        state=(m,mv,z,z,np.zeros_like(m),np.zeros_like(mv),np.zeros((*z.shape,4)),np.zeros_like(m),np.zeros_like(mv))
        initial=self.assemble(state,current)
        m2,m2v,_=self.system.init(initial[1])
        state=(*state[:4],m2,m2v,*state[6:])
        histories=[]
        count=round(end/dt)
        interval=round(.25/dt)
        for i in range(count+1):
            t=i*dt
            first,f2,df=self.assemble(state,current)
            if i%interval==0 or i==count:
                scale=np.linalg.norm(current[3][:,P.IC])
                r1,_=self.system.cres(state[0],state[1],current[3])
                residual,acc=self.system.cres(state[4],state[5],f2)
                power=np.linalg.norm(residual,axis=1)
                ids=np.argsort(power)[-8:][::-1]
                modes=[]
                for j in ids:
                    nhat=self.K[j]/np.linalg.norm(self.K[j])
                    v=residual[j,1:]
                    long=nhat@v
                    modes.append(dict(k=self.modes[j].tolist(),norm=float(power[j]),energy=float(abs(residual[j,0])),
                                      longitudinal=float(abs(long)),transverse=float(np.linalg.norm(v-long*nhat))))
                rms1=float(np.linalg.norm(state[0]))
                rms2=float(np.linalg.norm(state[4]))
                row=dict(t=t,linear=float(np.linalg.norm(r1)/scale),next_scaled=float(A.EPSILON*np.linalg.norm(power)/scale),
                         absolute=float(np.linalg.norm(power)),leading_source=float(scale),metric1=A.EPSILON*rms1,
                         metric2=A.EPSILON**2*rms2,metric_ratio=A.EPSILON*rms2/rms1,
                         tidal=float(np.linalg.norm(acc)/2),dominant=modes)
                histories.append(row)
                if abs(t-round(t))<1e-10:
                    print("t",t, "residual",row["next_scaled"],"ratio",row["metric_ratio"],"mode",modes[0]["k"],flush=True)
            if i==count:
                break
            mid,target=self.flat(t+dt/2),self.flat(t+dt)
            second=self.assemble(A.plus(state,first,dt/2),mid)[0]
            third=self.assemble(A.plus(state,second,dt/2),mid)[0]
            fourth=self.assemble(A.plus(state,third,dt),target)[0]
            state=tuple(x+dt/6*(a+2*b+2*c+d) for x,a,b,c,d in zip(state,first,second,third,fourth))
            current=target
        S=self.system
        zero=np.zeros_like(S.AiC)
        eye=np.broadcast_to(np.eye(6),zero.shape)
        generator=np.concatenate([np.concatenate([zero,eye],axis=-1),
                                  np.concatenate([S.AiC,S.AiB],axis=-1)],axis=-2)
        ev,vec=np.linalg.eig(generator)
        defect=generator@vec-vec*ev[:,None,:]
        constraint=np.concatenate([S.RC,S.RB],axis=-1)
        piv=np.unravel_index(ev.real.argmax(),ev.shape)
        mode_list={tuple(m["k"]) for row in histories for m in row["dominant"][:2]}
        stats=[]
        for j,mode in enumerate(self.modes):
            if tuple(mode) not in mode_list and j!=piv[0]:
                continue
            stats.append(dict(k=mode.tolist(),eigenvalues=[[float(x.real),float(x.imag)] for x in ev[j]],
                              eigen_residual=float(np.max(np.linalg.norm(defect[j],axis=0))),
                              constraint_norms=np.linalg.norm(constraint[j]@vec[j],axis=0).tolist(),
                              symbol_norm=float(np.linalg.norm(self.K[j]))))
        return dict(history=histories,spectrum=stats,max_vacuum_rate=float(ev.real.max()))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("--n",type=int,default=24)
    ap.add_argument("--dt",type=float,default=.05)
    ap.add_argument("--end",type=float,default=20)
    ap.add_argument("--continuum",action="store_true")
    ap.add_argument("--off",action="store_true")
    args=ap.parse_args()
    start=time.process_time()
    run=Experiment(args.n,not args.continuum,not args.off).run(args.dt,args.end)
    run.update(settings=vars(args),cpu_seconds=time.process_time()-start)
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35C"
    path.mkdir(parents=True,exist_ok=True)
    (path/f"{args.tag}.json").write_text(json.dumps(run,separators=(",",":"))+"\n")


if __name__=="__main__":
    main()
