"""GC35D: compare relaxed Skyrme stress across grid spacing and box size."""
import argparse
import json
from itertools import product
from pathlib import Path
import time

import exgc35b_refinement as B
import numpy as np

P, torch = B.P, B.torch
TOLERANCE = 1e-5


class Knot(B.Knot):
    def __init__(self,h,degree,half_width):
        # Resample the same smooth starts on the requested physical box.
        ref=B.Knot(h,degree)
        self.h,self.B=h,degree
        half=round(half_width/h)
        self.n=2*half+1
        x=h*np.arange(-half,half+1)
        self.x=np.stack(np.meshgrid(x,x,x,indexing="ij"),axis=-1)
        self.mask=torch.tensor(np.max(abs(self.x),axis=-1)<half*h-1e-12)
        self.vac=ref.vac
        if self.n==ref.n:
            self.initial=ref.initial
            return
        r=np.linalg.norm(self.x,axis=-1)
        X,Y,Z=np.moveaxis(self.x,-1,0)
        if degree==1:
            sample=ref.inflate(torch.tensor(ref.initial)).numpy()
            axis=sample.shape[0]//2
            profile=np.arccos(np.clip(sample[axis:-1,axis,axis,0],-1,1))
            radii=np.arange(len(profile))*h
            f=np.where(r<=radii[-1],np.interp(r,radii,profile),2.1516747905508162/np.maximum(r,1e-12)**2)
            unit=np.divide(self.x,r[...,None],out=np.zeros_like(self.x),where=r[...,None]>0)
        else:
            u=(X+1j*Y)/np.maximum(r+Z,1e-14)
            u=np.where((X==0)&(Y==0)&(Z<0),1e12,u)**2
            a=abs(u)**2
            unit=np.stack([2*u.real,2*u.imag,1-a],axis=-1)/(1+a)[...,None]
            f=np.pi*np.exp(-r*r/(2*1.3**2))
        p=np.concatenate((np.cos(f)[...,None],np.sin(f)[...,None]*unit),axis=-1)
        self.initial=p[self.mask.numpy()].ravel()

    def relax(self,budget):
        start=time.process_time()
        interior=self.mask.numpy().ravel()
        history=[]
        best=self.initial.copy()
        last=None
        calls=0
        iters=0

        def objective(x):
            nonlocal last,calls
            q=x.reshape(-1,4)
            norm=np.linalg.norm(q,axis=1)
            u=q/norm[:,None]
            p=np.zeros((self.n**3,4))
            p[:,0]=1
            p[interior]=u
            energy,g=B.energy_gradient(p,self.n,self.h)
            g=g[interior]
            tangent=g-u*np.sum(u*g,axis=1)[:,None]
            force=np.max(abs(tangent))/self.h**3
            last=(energy,float(force))
            calls+=1
            lam=1e-5
            return energy+lam*np.sum((norm-1)**2),(tangent/norm[:,None]+2*lam*(norm-1)[:,None]*u).ravel()

        class StopRelaxation(Exception):
            pass

        def checkpoint(x):
            nonlocal best,iters
            iters+=1
            best=x.copy()
            e,f=last
            elapsed=time.process_time()-start
            if iters%100==0 or f<=TOLERANCE or elapsed>=budget:
                row=dict(iteration=iters,energy=e,force=f,cpu=elapsed)
                history.append(row)
                print(row,flush=True)
            if f<=TOLERANCE or elapsed>=budget or e<.5*self.B:
                raise StopRelaxation

        objective(best)
        history.append(dict(iteration=0,energy=last[0],force=last[1],cpu=0.))
        try:
            answer=B.minimize(objective,best,method="L-BFGS-B",jac=True,callback=checkpoint,
                              options=dict(maxiter=2**31-1,maxfun=2**31-1,ftol=0.,gtol=0.,maxcor=8,maxls=45))
            best=answer.x
            reason=str(answer.message)
        except StopRelaxation:
            reason="force tolerance" if last[1]<=TOLERANCE else "topology loss suspected" if last[0]<.5*self.B else "per-case compute limit"
        objective(best)
        return self.inflate(torch.tensor(best)),dict(history=history,iterations=iters,evaluations=calls,
                                                      force=last[1],stop=reason,cpu_seconds=time.process_time()-start)

    def measure(self,p):
        h,n=self.h,self.n
        d2,d4=self.density(p)
        density=(d2+d4).numpy()
        dx=[(torch.roll(p,-1,a)-torch.roll(p,1,a))/(2*h) for a in range(3)]
        degree=float(torch.linalg.det(torch.stack([p,*dx],dim=-2)).sum()*h**3/(2*np.pi**2))
        T=0.
        for stencil in self.stencils(p):
            D=stencil@stencil.transpose(-1,-2)
            tr=D.diagonal(dim1=-2,dim2=-1).sum(-1)
            quartic=(tr**2-(D*D).sum((-1,-2)))/2
            T=T+2*D+2*(tr[...,None,None]*D-D@D)-(tr+quartic)[...,None,None]*torch.eye(3)
        T=(T/(8*B.NORM)).numpy()
        total=T.sum(axis=(0,1,2))*h**3
        modes=np.array(list(product(range(-5,6),repeat=3)))
        k=2*np.pi*modes/(n*h)
        keep=(np.linalg.norm(k,axis=1)<=2.05)&np.any(modes!=0,axis=1)
        modes,k=modes[keep],k[keep]
        ids=tuple((modes%n).T)
        K=2*np.sin(k*h/2)/h
        up=np.zeros((len(K),4,4),complex)
        up[:,0,0]=P.sf.fftn(density)[ids]/n**3
        for j,(a,b) in enumerate(P.N.PAIRS):
            arr=T[...,a,b]
            if a!=b:
                arr=(arr+np.roll(arr,-1,a)+np.roll(arr,-1,b)+np.roll(np.roll(arr,-1,a),-1,b))/4
            off=np.asarray(P.N.OFFS[4+j])
            up[:,a+1,b+1]=up[:,b+1,a+1]=P.sf.fftn(arr)[ids]/n**3*np.exp(-1j*k@(h*off))
        stress=up[:,1:,1:]
        source=up.reshape(-1,16)@P.N.PM*P.KAP
        kernel=P.kernel_batch(K,0.)
        pinv=np.linalg.pinv(kernel,rcond=1e-12)
        def forced(src):
            return float(np.linalg.norm(P.mv(kernel,P.mv(pinv,src))-src)/np.linalg.norm(src))
        only=np.zeros_like(source)
        only[:,0]=source[:,0]
        div=np.einsum("ma,mab->mb",K,stress)
        e=float(density.sum()*h**3)
        return dict(energy=e,E2=float(d2.sum()*h**3),E4=float(d4.sum()*h**3),degree=degree,
                    integrated_stress=total.tolist(),virial=float(np.trace(total)/e),
                    full_residual=forced(source),energy_residual=forced(only),
                    stress_divergence=float(np.linalg.norm(div)/np.linalg.norm(up[:,0,0])),
                    source_modes=len(K))


def metric_variation():
    rng=np.random.default_rng(3504)
    F=rng.normal(size=(3,4))
    D=F@F.T
    S=rng.normal(size=(3,3))
    S=(S+S.T)/2
    tr=np.trace(D)
    q=(tr**2-np.trace(D@D))/2
    T=2*D+2*(tr*D-D@D)-(tr+q)*np.eye(3)
    def energy(eps):
        g=np.eye(3)+eps*S
        M=np.linalg.solve(g,D)
        return np.sqrt(np.linalg.det(g))*(np.trace(M)+(np.trace(M)**2-np.trace(M@M))/2)
    eps=1e-5
    predicted=-np.sum(T*S)/2
    fd=(energy(eps)-energy(-eps))/(2*eps)
    return float(abs(fd-predicted)/abs(predicted))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("tag")
    ap.add_argument("h",type=float)
    ap.add_argument("degree",type=int,choices=(1,2))
    ap.add_argument("--half",type=float,default=6)
    ap.add_argument("--cpu",type=float,default=650)
    args=ap.parse_args()
    knot=Knot(args.h,args.degree,args.half)
    p,relax=knot.relax(args.cpu)
    result=dict(settings=vars(args),n=knot.n,half_actual=(knot.n-1)*args.h/2,
                metric_variation_relative=metric_variation(),relaxation=relax,**knot.measure(p))
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35D"
    path.mkdir(parents=True,exist_ok=True)
    (path/f"{args.tag}.json").write_text(json.dumps(result,indent=2)+"\n")
    print({k:v for k,v in result.items() if k!="relaxation"},flush=True)


if __name__=="__main__":
    main()
