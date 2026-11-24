"""EX-GC35A: finite-time two-way scalar response through second gravitational order."""
import os
for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[name] = "1"

import json
import platform
import time
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid

import exgc34r as R

P = R.P
P.LBOX, P.SIG, P.NCUT = 12.0, 2.0, 8
EPSILON = 1e-3
ETA = P.N.ETA
PAIR = P.N.PAIRS


def plus(a, b, factor):
    return tuple(x+factor*y for x,y in zip(a,b))


def tensor(data):
    out = np.zeros((*data.shape[:-1],4,4))
    for j,(a,b) in enumerate(PAIR):
        out[...,a+1,b+1] = out[...,b+1,a+1] = data[...,j]
    return out


def connection(dh, inverse):
    lower = (dh.swapaxes(-3,-2)+np.moveaxis(dh,-3,-1)-dh)/2
    return inverse @ lower.reshape(*lower.shape[:-3],4,16)


def apply_inverse(inv, gamma):
    return (inv @ gamma.reshape(*gamma.shape[:-3],4,16)).reshape(gamma.shape)


class Experiment:
    def __init__(self,n,lattice):
        self.n,self.lattice = n,lattice
        self.grid = R.StandingPacket(n,"generic",lattice)
        h = self.grid.h
        for k in self.grid.KK:
            self.grid.a[abs(k*P.LBOX/(2*np.pi))>2.01] = 0
        modes = P.mask_modes()
        k,self.K = P.lattice_K(modes,h)
        if not lattice:
            self.K=k
        self.ids=tuple((modes%n).T)
        self.system=P.System(self.K)
        self.offset=np.asarray(P.N.OFFS)
        self.phase=np.exp(1j*k @ (h*self.offset).T) if lattice else np.ones((len(k),10))
        self.symbols=[2*np.sin(v*h/2)/h if lattice else v for v in self.grid.KK]
        self.cube = np.maximum.reduce(np.meshgrid(*([abs(np.arange(n)*h-P.LBOX/2)]*3),indexing="ij"))<3

    def differentiate(self,x,i):
        fac = self.symbols[i][(...,)+(None,)*(x.ndim-3)]
        return P.sf.ifftn(1j*fac*P.sf.fftn(x,axes=(0,1,2)),axes=(0,1,2)).real

    def metric(self,coeff,placed=False):
        spectra=np.zeros((self.n,)*3+(6,),complex)
        spectra[self.ids]=coeff*(self.phase[:,4:] if placed else 1)*self.n**3
        return P.sf.ifftn(spectra,axes=(0,1,2)).real

    def spatial_source_hat(self,tensor_up):
        coord=tensor_up.reshape(*tensor_up.shape[:3],16) @ P.N.PM
        return P.sf.fftn(coord,axes=(0,1,2))[self.ids]/self.n**3*P.KAP

    def lattice_source(self,phi,vel):
        n,b,d,o=P.O.msources2(phi,vel,self.grid.h)
        raw=np.stack([n,*b,*d,*o],axis=-1)
        fft=P.sf.fftn(raw,axes=(0,1,2))/self.n**3/self.grid.h**3
        source=2*P.KAP*fft[self.ids]*self.phase.conj()
        # Convert dephased coefficients to a common site for nonlinear products.
        collocated=np.zeros_like(fft)
        collocated[self.ids]=fft[self.ids]*self.phase.conj()
        collocated[0,0,0]=fft[0,0,0]
        raw=P.sf.ifftn(collocated*self.n**3,axes=(0,1,2)).real
        up=tensor(raw[...,4:])
        for i,(a,b) in enumerate(PAIR):
            if a==b:
                up[...,a+1,a+1]*=2
        up[...,0,0]=-raw[...,0]
        up[...,0,1:]=up[...,1:,0]=raw[...,1:4]
        return up,source

    def scalar_source(self,phi,vel):
        if self.lattice:
            return self.lattice_source(phi,vel)
        grad=np.stack([vel,*[self.differentiate(phi,i) for i in range(3)]],axis=-1)
        upper=grad@ETA
        up=upper[...,None]*upper[...,None,:]-.5*(np.sum(grad*upper,axis=-1)+phi**2)[...,None,None]*ETA
        return up,self.spatial_source_hat(up)

    def geometry(self,m,mv,ma):
        H, Hd, Hdd=[tensor(self.metric(v)) for v in (m,mv,ma)]
        dh=np.stack([Hd,*[self.differentiate(H,a) for a in range(3)]],axis=-3)
        dhdot=np.stack([Hdd,*[self.differentiate(Hd,a) for a in range(3)]],axis=-3)
        inv=-ETA@H@ETA
        invdot=-ETA@Hd@ETA
        g1=connection(dh,ETA).reshape(H.shape[:3]+(4,4,4))
        g1dot=connection(dhdot,ETA).reshape(g1.shape)
        g2=apply_inverse(inv@ETA,g1)
        g2dot=apply_inverse(invdot@ETA,g1)+apply_inverse(inv@ETA,g1dot)

        def ricci(g,gd):
            out=gd[...,0,:,:].copy()
            tr=np.einsum("...ama->...m",g)
            out[..., :,0]-=np.einsum("...ama->...m",gd)
            for a in range(1,4):
                out+=self.differentiate(g[...,a,:,:],a-1)
                out[..., :,a]-=self.differentiate(tr,a-1)
            return out

        r1=ricci(g1,g1dot)
        r2=ricci(g2,g2dot)+np.einsum("...rrs,...smn->...mn",g1,g1)-np.einsum("...rns,...smr->...mn",g1,g1)
        tr1=np.einsum("ab,...ab->...",ETA,r1)
        tr2=np.einsum("ab,...ab->...",ETA,r2)+np.sum(inv*r1,axis=(-1,-2))
        volume=np.trace(H,axis1=-2,axis2=-1)/2
        Einstein1=ETA@r1@ETA-.5*tr1[...,None,None]*ETA
        density2=(ETA@r2@ETA+inv@r1@ETA+ETA@r1@inv-.5*(tr2[...,None,None]*ETA+tr1[...,None,None]*inv)
                  +volume[...,None,None]*Einstein1)
        return inv,volume,g1,density2

    def flat(self,t):
        phi,vel=self.grid.field(t)
        up,f=self.scalar_source(phi,vel)
        return phi,vel,up,f

    def assemble(self,state,flat):
        m,mv,p,pv,m2,m2v,comp,delta,delta_v=state
        phi,vel,T0,f0=flat
        ma=self.system.acc(m,mv,f0[:,P.IH])
        inv,vol,gamma,G2=self.geometry(m,mv,ma)
        grad=np.stack([vel,*[self.differentiate(phi,i) for i in range(3)]],axis=-1)
        du=np.einsum("...ij,...j->...i",inv,grad)
        up=grad@ETA
        trace=np.sum(grad*du,axis=-1)
        kinetic=np.sum(grad*up,axis=-1)+phi**2
        metric_T=(du[...,None]*up[...,None,:]+up[...,None]*du[...,None,:]
                  -.5*(trace[...,None,None]*ETA+kinetic[...,None,None]*inv)+vol[...,None,None]*T0)
        tplus,dfplus=self.scalar_source(phi+p,vel+pv)
        tminus,dfminus=self.scalar_source(phi-p,vel-pv)
        dtensor=(tplus-tminus)/2
        df=(dfplus-dfminus)/2
        C0=-T0.mean(axis=(0,1,2))
        compdot=-np.einsum("...vab,ab->...v",gamma,C0)
        for a in range(3):
            compdot[...,0]-=self.differentiate(comp[...,a+1],a)
        C1=np.zeros_like(G2)
        C1[...,0,:]=C1[..., :,0]=comp
        effective=dtensor+metric_T+C1-(G2+G2.swapaxes(-1,-2))/2
        f2=df+self.spatial_source_hat(metric_T+C1-(G2+G2.swapaxes(-1,-2))/2)
        H=self.metric(m,True).transpose(3,0,1,2)
        Hdot=self.metric(mv,True).transpose(3,0,1,2)
        field_lap=R.lap(phi,self.grid,self.lattice)
        driving=R.spatial(phi,H,self.grid,self.lattice)-sum(H[:3])/2*field_lap-sum(Hdot[:3])/2*vel
        dp=R.lap(p,self.grid,self.lattice)-p+driving
        derivative=(mv,ma,pv,dp,m2v,self.system.acc(m2,m2v,f2[:,P.IH]),compdot,
                    delta_v,self.system.acc(delta,delta_v,df[:,P.IH]))
        return derivative,f2,df,effective

    def run(self,dt):
        current=self.flat(0)
        m,mv,_=self.system.init(current[3])
        zeros=np.zeros_like(current[0])
        state=(m,mv,zeros,zeros,np.zeros_like(m),np.zeros_like(mv),
               np.zeros((*zeros.shape,4)),np.zeros_like(m),np.zeros_like(mv))
        initial=self.assemble(state,current)
        second,second_v,_=self.system.init(initial[1])
        state=(*state[:4],second,second_v,*state[6:])
        summary,tidal,energies,fluxes=[],{},[],[]
        peak=0.0
        for i in range(round(4/dt)+1):
            t=i*dt
            first,f2,df,effective=self.assemble(state,current)
            scale=np.linalg.norm(current[3][:,P.IC])
            r1,_=self.system.cres(state[0],state[1],current[3])
            r2,acc=self.system.cres(state[4],state[5],f2)
            bad,tangent=self.system.cres(state[7],state[8],df)
            cube=self.cube
            mass=np.sum((effective[...,0,0]-effective[...,0,0].mean())[cube])*self.grid.h**3
            outflow=sum(np.sum(self.differentiate(effective[...,0,a+1],a)[cube])
                        for a in range(3))*self.grid.h**3
            energies.append(float(mass))
            fluxes.append(float(outflow))
            peak=max(peak,EPSILON*float(np.max(abs(self.metric(state[0])))))
            summary.append(dict(t=t,linear=float(np.linalg.norm(r1)/scale),
                                next_scaled=float(EPSILON*np.linalg.norm(r2)/scale),
                                direct_scaled=float(EPSILON*np.linalg.norm(bad)/scale),
                                combined=float(np.linalg.norm(r1+EPSILON*r2)/scale)))
            if t in (2,4):
                tidal[t]=-.5*acc
                summary[-1].update(tidal_norm=float(np.linalg.norm(acc)/2),
                                   reciprocal_tidal_norm=float(np.linalg.norm(tangent)/2),
                                   mean_effective=np.mean(effective,axis=(0,1,2)).tolist())
            if i==round(4/dt):
                break
            mid,end=self.flat(t+dt/2),self.flat(t+dt)
            second=self.assemble(plus(state,first,dt/2),mid)[0]
            third=self.assemble(plus(state,second,dt/2),mid)[0]
            fourth=self.assemble(plus(state,third,dt),end)[0]
            state=tuple(x+dt/6*(a+2*b+2*c+d) for x,a,b,c,d in zip(state,first,second,third,fourth))
            current=end
        flux=dict(initial=energies[0],end=energies[-1],integral=float(trapezoid(fluxes,dx=dt)),
                  times=[i*dt for i in range(len(energies))],energy=energies,outflow=fluxes)
        result=dict(n=self.n,dt=dt,lattice=self.lattice,history=summary,flux=flux,metric_amplitude=peak)
        print(self.n,dt,self.lattice,summary[-1],flush=True)
        return result,tidal


def main():
    start=time.process_time()
    (coarse,a),(fine,b),(half,_),(reference,c)=[Experiment(n,lattice).run(dt)
                                              for n,dt,lattice in ((24,.05,True),(48,.05,True),(24,.025,True),(24,.05,False))]
    maximum=lambda run,key:max(r[key] for r in run["history"])
    spatial=[dict(t=t,coarse=P.tidal_err(a[t],c[t]),fine=P.tidal_err(b[t],c[t])) for t in (2,4)]
    cf,ff,hf=(x["history"][-1]["next_scaled"] for x in (coarse,fine,half))
    gates=dict(A0=maximum(reference,"next_scaled")/EPSILON<=1e-3 and
                  max(maximum(x,"linear") for x in (coarse,fine,half))<=1e-5,
               A1=maximum(fine,"next_scaled")<=1e-3 and cf/ff>=2 and abs(cf/hf-1)<=.1,
               A2=all(x["fine"]<=.15 and x["coarse"]>=2*x["fine"] for x in spatial)
                  and all(r["reciprocal_tidal_norm"]>0 for r in fine["history"] if "reciprocal_tidal_norm" in r))
    result=dict(gates=gates,coarse=coarse,fine=fine,half=half,continuum=reference,spatial=spatial,
                final_residual_ratio=cf/ff,timestep_fraction=abs(cf/hf-1),
                cpu_seconds=time.process_time()-start,python=platform.python_version(),numpy=np.__version__)
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35A"
    path.mkdir(parents=True,exist_ok=True)
    (path/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(gates, spatial, flush=True)


if __name__=="__main__":
    main()
