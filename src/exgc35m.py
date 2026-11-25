"""GC35M: local synchronous Y1 damping on the compensated scalar packet."""
import argparse
import json
from pathlib import Path
import time
import exgc35e as E
import numpy as np

P,A=E.P,E.A
ALPHA=BETA=.12
TRACE=np.array([1.,1.,1.,0.,0.,0.])


class Damped(E.Projected):
    def __init__(self,n):
        super().__init__(n,True)
        self.adjoint=self.system.RB.conj().swapaxes(1,2)/P.WT[None,:,None]

    def correction(self,h,v,f):
        residual,_=self.system.cres(h,v,f)
        return ALPHA*residual[:,:1]*TRACE,-BETA*P.mv(self.adjoint,residual)

    def assemble(self,state,flat):
        d,src,df=super().assemble(state,flat)
        dh,dv=self.correction(state[4],state[5],src)
        return (*d[:4],d[4]+dh,d[5]+dv,*d[6:]),src,df

    def stability(self,dt):
        # Check the whole physical Brillouin cube, including the near-corner norm.
        kmax=2/self.h
        axis=np.linspace(-kmax,kmax,9)
        K=np.stack(np.meshgrid(axis,axis,axis,indexing="ij"),axis=-1).reshape(-1,3)
        S=P.System(K)
        rb=S.RB
        adj=rb.conj().swapaxes(1,2)/P.WT[None,:,None]
        eye=np.broadcast_to(np.eye(6),S.AiC.shape)
        upper=ALPHA*TRACE[None,:,None]*S.RC[:,:1,:]
        lower=S.AiC-BETA*adj@S.RC
        upper_right=eye+ALPHA*TRACE[None,:,None]*rb[:,:1,:]
        lower_right=-BETA*adj@rb+S.AiB
        block=np.concatenate([np.concatenate([upper,upper_right],axis=2),np.concatenate([lower,lower_right],axis=2)],axis=1)
        z=dt*np.linalg.eigvals(block)
        return dict(max_rk4_amp=float(np.max(abs(1+z+z*z/2+z**3/6+z**4/24))),
                    max_real=float(np.max(z.real)/dt))

    def evolve(self,dt,end):
        S=self.system
        current=self.flat(0)
        m,mv,_=S.init(current[3])
        zero=np.zeros_like(current[0])
        state=(m,mv,zero,zero,np.zeros_like(m),np.zeros_like(mv),np.zeros((*zero.shape,4)),np.zeros_like(m),np.zeros_like(mv))
        src=self.assemble(state,current)[1]
        m2,v2,_=S.init(src)
        state=(*state[:4],m2,v2,*state[6:])
        rows=[]
        size=round(end/dt)
        at={round(t/dt):t for t in E.READOUTS if t<=end}
        needed={j for i in at for j in (i-3,i-2,i-1,i,i+1) if 0<=j<=size}
        positions={}
        for i in range(size+1):
            first,f2,_=self.assemble(state,current)
            x,v=state[4:6]
            res,_=S.cres(x,first[4],f2)
            aux,_=S.cres(x,v,f2)
            dh,dv=self.correction(x,v,f2)
            factor=A.EPSILON/np.linalg.norm(current[3][:,P.IC])
            row=dict(t=i*dt,physical=float(factor*np.linalg.norm(res)),
                     auxiliary=float(factor*np.linalg.norm(aux)),
                     rhs_h=float(np.linalg.norm(dh)),rhs_v=float(np.linalg.norm(dv)))
            rows.append(row)
            if i in needed:
                positions[i]=x.copy()
            if i in at:
                print(row,flush=True)
            if i==size:
                break
            mid,target=self.flat((i+.5)*dt),self.flat((i+1)*dt)
            second=self.assemble(A.plus(state,first,dt/2),mid)[0]
            third=self.assemble(A.plus(state,second,dt/2),mid)[0]
            fourth=self.assemble(A.plus(state,third,dt),target)[0]
            state=tuple(x+dt/6*(a+2*b+2*c+d) for x,a,b,c,d in zip(state,first,second,third,fourth))
            current=target
        arrays={}
        for i,t in at.items():
            dd=(positions[i+1]-2*positions[i]+positions[i-1]) if i<size else (
                2*positions[i]-5*positions[i-1]+4*positions[i-2]-positions[i-3])
            arrays[str(t)]=-dd/(2*dt**2)
        return rows,arrays


if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("n",type=int)
    ap.add_argument("--reference",default="EX-GC35E_verified")
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[1]/"results"
    out=root/"EX-GC35M"
    out.mkdir(parents=True,exist_ok=True)
    start=time.process_time()
    x=Damped(args.n)
    check=x.stability(.05)
    assert check["max_rk4_amp"]<1.00001,check
    rows,readouts=x.run(.05,20)
    with np.load(root/args.reference/"continuum.npz") as ref:
        tides={k:float(P.tidal_err(v,ref[k+"_reconstructed"])) for k,v in readouts.items()}
    result=dict(n=args.n,coefficients=[ALPHA,BETA],stability=check,history=rows,tides=tides,cpu_seconds=time.process_time()-start)
    (out/f"n{args.n}.json").write_text(json.dumps(result,separators=(",",":"))+"\n")
    print(tides,flush=True)
