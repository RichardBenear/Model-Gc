"""EX-GC35B: static GC26-style knots and leading initial response on Y1."""
import os
for key in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS"):
    os.environ[key]="1"

import json
import platform
import time
from itertools import product
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

import exgc34p as P

torch.set_num_threads(1)
torch.set_default_dtype(torch.float64)
NORM=12*np.pi**2
ORIENT=tuple(product((-1,1),repeat=3))
MASS_EPS=1e-4


class Knot:
    def __init__(self,h,B):
        self.h=h
        self.n=round(12/h)+1
        x=np.linspace(-6,6,self.n)
        self.x=np.stack(np.meshgrid(x,x,x,indexing="ij"),axis=-1)
        r=np.linalg.norm(self.x,axis=-1)
        X,Y,Z=np.moveaxis(self.x,-1,0)
        if B==1:
            unit=np.divide(self.x,r[...,None],out=np.zeros_like(self.x),where=r[...,None]>0)
        else:
            u=(X+1j*Y)/np.maximum(r+Z,1e-14)
            u=np.where((X==0)&(Y==0)&(Z<0),1e12,u)**2
            a=abs(u)**2
            unit=np.stack([2*u.real,2*u.imag,1-a],axis=-1)/(1+a)[...,None]
        f=np.pi*np.exp(-r*r/(2*1.3**2))
        p=np.concatenate((np.cos(f)[...,None],np.sin(f)[...,None]*unit),axis=-1)
        self.mask=torch.tensor(np.max(abs(self.x),axis=-1)<6-1e-12)
        self.vac=torch.tensor([1.,0,0,0])
        self.initial=p[self.mask.numpy()].ravel()
        self.B=B

    def inflate(self,x):
        arr=self.vac.expand(self.n,self.n,self.n,4).clone()
        arr[self.mask]=torch.nn.functional.normalize(x.reshape(-1,4),dim=-1)
        return arr

    def stencils(self,p):
        for signs in ORIENT:
            yield torch.stack([(torch.roll(p,-s,a)-p)*(s/self.h)
                               for a,s in enumerate(signs)],dim=-2)

    def density(self,p,metric=None,clock=None):
        if metric is None:
            inv=torch.eye(3)
            weight=1.
        else:
            inv=torch.linalg.inv(metric)
            weight=clock*torch.sqrt(torch.linalg.det(metric))
        e2=e4=0.
        for stencil in self.stencils(p):
            D=stencil@stencil.transpose(-1,-2)
            mixed=inv@D
            tr=mixed.diagonal(dim1=-2,dim2=-1).sum(-1)
            e2=e2+tr
            e4=e4+.5*(tr**2-(mixed*mixed.transpose(-1,-2)).sum((-1,-2)))
        return weight*e2/(8*NORM),weight*e4/(8*NORM)

    def objective(self,x):
        u=torch.tensor(x,requires_grad=True)
        e2,e4=self.density(self.inflate(u))
        total=(e2+e4).sum()*self.h**3
        total.backward()
        return total.item(),u.grad.numpy()

    def rest(self):
        answer=minimize(self.objective,self.initial,method="L-BFGS-B",jac=True,
                        options=dict(maxiter=550,maxcor=8,ftol=2e-10,gtol=3e-7))
        p=self.inflate(torch.tensor(answer.x))
        e2,e4=self.density(p)
        dx=[(torch.roll(p,-1,a)-torch.roll(p,1,a))/(2*self.h) for a in range(3)]
        baryon=torch.linalg.det(torch.stack([p,*dx],dim=-2))
        density=(e2+e4).numpy()
        T=0.
        for stencil in self.stencils(p):
            D=stencil@stencil.transpose(-1,-2)
            tr=D.diagonal(dim1=-2,dim2=-1).sum(-1)
            quartic=(tr**2-(D*D).sum((-1,-2)))/2
            T=T+2*D+2*(tr[...,None,None]*D-D@D)-(tr+quartic)[...,None,None]*torch.eye(3)
        info=dict(B=self.B,h=self.h,n=self.n,E=float(density.sum()*self.h**3),E2=float(e2.sum()*self.h**3),
                  E4=float(e4.sum()*self.h**3),degree=float(baryon.sum()*self.h**3/(2*np.pi**2)),
                  iterations=int(answer.nit),success=bool(answer.success),optimizer_message=str(answer.message),
                  gradient_max=float(np.max(abs(answer.jac))))
        print("rest",info,flush=True)
        return p,density,(T/(8*NORM)).numpy(),info

    def gravity(self,density,stress):
        n,h=self.n,self.h
        modes=np.arange(-4,5)
        triads=np.array(list(product(modes,repeat=3)))
        triads=triads[np.any(triads!=0,axis=1)]
        ids=tuple((triads%n).T)
        k=2*np.pi*triads/(n*h)
        K=2*np.sin(k*h/2)/h
        up=np.zeros((len(k),4,4),complex)
        up[:,0,0]=P.sf.fftn(density)[ids]/n**3
        phases=np.exp(-1j*k @ (h*np.asarray(P.N.OFFS[4:])).T)
        for j,(a,b) in enumerate(P.N.PAIRS):
            if a==b:
                arr=stress[...,a,b]
            else:
                arr=(stress[...,a,b]+np.roll(stress[...,a,b],-1,a)+np.roll(stress[...,a,b],-1,b)
                    +np.roll(np.roll(stress[...,a,b],-1,a),-1,b))/4
            up[:,a+1,b+1]=up[:,b+1,a+1]=P.sf.fftn(arr)[ids]/n**3*phases[:,j]
        source=up.reshape(-1,16) @ P.N.PM*P.KAP
        M=P.kernel_batch(K,0.)
        inv=np.linalg.pinv(M,rcond=1e-12)
        full=P.mv(inv,source)
        full_res=float(np.linalg.norm(P.mv(M,full)-source)/np.linalg.norm(source))
        only=source.copy()
        only[:,4:]=0
        energy=P.mv(inv,only)
        continuum=P.mv(np.linalg.pinv(P.kernel_batch(k,0),rcond=1e-12),only)
        kept=np.linalg.norm(k,axis=1)<=1.5
        tidal=-energy[:,0,None,None]*K[:,:,None]*K[:,None,:]
        tidal_ref=-continuum[:,0,None,None]*k[:,:,None]*k[:,None,:]
        error=float(np.linalg.norm((tidal-tidal_ref)[kept])/np.linalg.norm(tidal_ref[kept]))
        coeff=full if full_res<=1e-3 else energy
        spectrum=np.zeros((n,n,n,10),complex)
        spectrum[ids]=coeff*n**3
        metric=P.sf.ifftn(spectrum,axes=(0,1,2)).real
        reconstructed=-2*np.sum(K*K,axis=1)*energy[:,0]/P.KAP
        lowest=np.linalg.norm(triads,axis=1)==1
        form=up[:,0,0]/(density.mean())
        monopole=float(np.mean((reconstructed[lowest]/form[lowest]).real)*(n*h)**3)
        return metric,dict(full_residual=full_res,energy_residual=float(np.linalg.norm(P.mv(M,energy)-only)/np.linalg.norm(only)),
                           tidal_error=error,monopole=monopole,feedback_source="full" if full_res<=1e-3 else "energy_only")

    def feedback(self,p,geometry):
        d=torch.tensor(geometry[...,0])
        mat=torch.zeros((*d.shape,3,3))
        for j,(a,b) in enumerate(P.N.PAIRS):
            mat[...,a,b]=mat[...,b,a]=torch.tensor(geometry[...,4+j])
        eye=torch.eye(3)
        base=p.detach().clone().requires_grad_()

        def accel(strength):
            g=eye+strength*mat
            clock=1+strength*d
            inv=torch.linalg.inv(g)
            e2,e4=self.density(base,g,clock)
            grad=torch.autograd.grad((e2+e4).sum(),base)[0]*NORM
            inertia=0.
            for stencil in self.stencils(base):
                mixed=inv@stencil
                tr=(stencil*mixed).sum((-1,-2))
                inertia=inertia+(1+tr)[...,None,None]*torch.eye(4)-stencil.transpose(-1,-2)@mixed
            inertia=inertia * (torch.sqrt(torch.linalg.det(g))/(4*clock))[...,None,None]
            v=torch.linalg.solve(inertia,-grad[...,None])[...,0]
            normal=torch.linalg.solve(inertia,base[...,None])[...,0]
            force=v-normal*((v*base).sum(-1)/(normal*base).sum(-1))[...,None]
            return force[self.mask].detach().numpy()

        estimates=[(accel(eps)-accel(-eps))/(2*eps) for eps in (MASS_EPS,MASS_EPS/2)]
        return dict(norm=float(np.sqrt(self.h**3)*np.linalg.norm(estimates[1])),
                    finite_relative=float(np.linalg.norm(estimates[0]-estimates[1])/np.linalg.norm(estimates[1])),
                    metric_strength=float(MASS_EPS*np.max(abs(geometry))))


def main():
    start=time.process_time()
    rows=[]
    for h in (.375,.3):
        for B in (1,2):
            knot=Knot(h,B)
            p,rho,T,row=knot.rest()
            metric,checks=knot.gravity(rho,T)
            row.update(checks,feedback=knot.feedback(p,metric))
            print("gravity",row,flush=True)
            rows.append(row)
    fine=rows[2:]
    ratios=dict(energy=fine[1]["E"]/fine[0]["E"],mass=fine[1]["monopole"]/fine[0]["monopole"])
    gates=dict(B0=all(abs(abs(row["degree"])-row["B"])<=.35 for row in fine) and ratios["energy"]<2,
               B1_control=all(row["energy_residual"]<=1e-10 and row["tidal_error"]<=.05 for row in rows)
                          and abs(ratios["mass"]/ratios["energy"]-1)<=1e-10,
               B1_full=all(row["full_residual"]<=1e-3 for row in fine),
               B2=all(row["feedback"]["norm"]>0 and row["feedback"]["finite_relative"]<=.01 for row in rows))
    res=dict(rows=rows,gates=gates,fine_ratios=ratios,cpu_seconds=time.process_time()-start,
             python=platform.python_version(),torch=torch.__version__,numpy=np.__version__)
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35B"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(res,indent=2)+"\n")
    print(gates,ratios,flush=True)


if __name__=="__main__":
    main()
