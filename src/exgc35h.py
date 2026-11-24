"""GC35H: a staggered Skyrme stress candidate and a chain-rule counterexample."""
import argparse
from itertools import product
import json
from pathlib import Path
import time
import exgc35d as D
import numpy as np

P, B = D.P, D.B
PAIRS = P.N.PAIRS
NORM = B.NORM


def fw(a,axis,step=1):
    return np.roll(a,-step,axis)


def dot(a,b):
    return np.sum(a*b,axis=-1)


def center_face(x,a,b):
    return (x+fw(x,a)+fw(x,b)+fw(fw(x,a),b))/4


def gradients(q,h):
    gp=[(fw(q,a)-q)/h for a in range(3)]
    gm=[fw(gp[a],a,-1) for a in range(3)]
    return gp,gm


def components(q,h):
    gp,gm=gradients(q,h)
    squares=[(dot(u,u)+dot(v,v))/2 for u,v in zip(gp,gm)]
    rho2=sum(squares)
    sigma=[dot(gp[a],gm[a])-sum(squares[b] for b in range(3) if b!=a) for a in range(3)]
    area={}
    shears={}
    for a,b in PAIRS[3:]:
        area[a,b]=sum(dot(u,u)*dot(v,v)-dot(u,v)**2
                      for u in (gp[a],gm[a]) for v in (gp[b],gm[b]))/4
        u=(gp[a]+fw(gp[a],b))/2
        v=(gp[b]+fw(gp[b],a))/2
        c=3-a-b
        ws=(center_face(gp[c],a,b),center_face(gm[c],a,b))
        shears[a,b]=[np.stack((u,v,w),axis=-2) for w in ws]
        sigma.append(2*dot(u,v))
    rho4=sum(area.values())
    quartic=[sum((1 if c in key else -1)*qarea for key,qarea in area.items()) for c in range(3)]
    for ff in shears.values():
        terms=[]
        for f in ff:
            u,v,w=f[...,0,:],f[...,1,:],f[...,2,:]
            terms.append(dot(w,w)*dot(u,v)-dot(u,w)*dot(v,w))
        quartic.append(sum(terms))
    return rho2/NORM,rho4/NORM,np.stack(sigma,axis=-1)/NORM,np.stack(quartic,axis=-1)/NORM,area,shears


def divergence(T,h):
    rows=[(fw(T[...,a],a)-T[...,a])/h for a in range(3)]
    for j,(a,b) in enumerate(PAIRS[3:],start=3):
        x=T[...,j]
        rows[a]+=(x-fw(x,b,-1))/h
        rows[b]+=(x-fw(x,a,-1))/h
    return np.stack(rows,axis=-1)


def audit():
    rng=np.random.default_rng(3508)
    h=.2
    q=rng.normal(size=(6,6,6,4))
    q/=np.linalg.norm(q,axis=-1,keepdims=True)
    _,_,sigma,quartic,area,shears=components(q,h)
    u=rng.normal(size=(*q.shape[:3],6))
    vol=h**3
    weights=np.array([.5]*3+[1.]*3)
    def energy(eps):
        d=eps*u[...,:3]
        o=eps*u[...,3:]
        lap=np.zeros(q.shape[:3])
        xyz=np.zeros((3,*q.shape[:3]))
        matter=sum((-2)*P.O.mlagrangian2(v,lap,lap,xyz,np.moveaxis(d,-1,0),np.moveaxis(o,-1,0),h)
                   for v in np.moveaxis(q,-1,0))
        mass=P.O.MASS**2*np.sum((1+d.sum(-1)/2)*dot(q,q))*vol
        g=1+d
        e4=sum(np.sum(np.sqrt(np.prod(g,axis=-1))*Q/(g[...,a]*g[...,b])) for (a,b),Q in area.items())
        for j,ff in enumerate(shears.values()):
            g=np.broadcast_to(np.eye(3),(*q.shape[:3],3,3)).copy()
            g[...,0,1]=g[...,1,0]=o[...,j]
            ginv=np.linalg.inv(g)
            weight=np.sqrt(np.linalg.det(g))
            for f in ff:
                gram=f@f.swapaxes(-1,-2)
                m=ginv@gram
                tr=np.trace(m,axis1=-1,axis2=-2)
                qq=(tr**2-np.trace(m@m,axis1=-1,axis2=-2))/2
                e4+=.5*np.sum(weight*qq)
        return np.array([(matter-mass)/NORM,e4*vol/NORM])
    delta=1e-5
    fd=(energy(delta)-energy(-delta))/(2*delta)
    expect=-vol*np.array([np.sum(weights*u*T) for T in (sigma,quartic)])
    gp,_=gradients(q,h)
    lap=P.N.lap(q,h)
    predicted=np.stack([dot(lap+fw(lap,a),gp[a])/NORM for a in range(3)],axis=-1)
    measured=divergence(sigma,h)
    err=np.linalg.norm(measured-predicted)/np.linalg.norm(predicted)
    p=np.zeros((7,3,3,4))
    p[...,0]=np.array([1,-1,1,1,1,1,1])[:,None,None]
    _,_,sig4,sta4,_,_=components(p,h)
    lap=P.N.lap(p,h)
    lam=dot(lap,p)
    force=lap-lam[...,None]*p
    reaction=np.stack([(fw(lam,a)-lam)*(1-dot(p,fw(p,a)))/(h*NORM) for a in range(3)],axis=-1)
    return dict(finite_relative=(abs(fd-expect)/np.maximum(abs(expect),1e-30)).tolist(),
                sigma_unconstrained_ward=float(err),rank_one_force=float(np.max(abs(force))),
                rank_one_quartic=float(np.max(abs(sta4))),rank_one_ward=float(np.max(abs(divergence(sig4,h)))),
                rank_one_identity_difference=float(np.max(abs(divergence(sig4,h)-reaction))))


def forced(density,stress,h):
    n=density.shape[0]
    modes=np.array(list(product(range(-5,6),repeat=3)))
    k=2*np.pi*modes/(n*h)
    keep=(np.linalg.norm(k,axis=1)<=2.05)&np.any(modes!=0,axis=1)
    modes,k=modes[keep],k[keep]
    ids=tuple((modes%n).T)
    K=2*np.sin(k*h/2)/h
    up=np.zeros((len(k),4,4),complex)
    up[:,0,0]=P.sf.fftn(density)[ids]/n**3
    for j,(a,b) in enumerate(PAIRS):
        off=np.asarray(P.N.OFFS[4+j])
        up[:,a+1,b+1]=up[:,b+1,a+1]=P.sf.fftn(stress[...,j])[ids]/n**3*np.exp(-1j*k@(h*off))
    src=up.reshape(-1,16)@P.N.PM*P.KAP
    kernel=P.kernel_batch(K,0.)
    returned=P.mv(kernel,P.mv(np.linalg.pinv(kernel,rcond=1e-12),src))
    return float(np.linalg.norm(returned-src)/np.linalg.norm(src))


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("degree",type=int,choices=(1,2))
    a=parser.parse_args()
    start=time.process_time()
    check=audit()
    assert max(check["finite_relative"])<1e-6 and check["sigma_unconstrained_ward"]<1e-12,check
    h=.2
    K=D.Knot(h,a.degree,6)
    p,relaxation=K.relax(650)
    legacy=K.measure(p)
    q=p.numpy()
    rho2,rho4,sigma,quartic,_,_=components(q,h)
    rho=rho2+rho4
    mass=float(rho.sum()*h**3)
    flat_error=abs(mass/legacy["energy"]-1)
    c2,c4=divergence(sigma,h),divergence(quartic,h)
    result=dict(degree=a.degree,legacy=legacy,relaxation=relaxation,audit=check,
                energy=mass,flat_energy_relative=flat_error,native_residual=forced(rho,sigma+quartic,h),
                sigma_only_stress_residual=forced(rho,sigma,h),quartic_only_stress_residual=forced(rho,quartic,h),
                density_only_residual=forced(rho,np.zeros_like(sigma),h),
                divergence_sigma=float(np.linalg.norm(c2)),divergence_quartic=float(np.linalg.norm(c4)),
                divergence_total=float(np.linalg.norm(c2+c4)),cpu_seconds=time.process_time()-start)
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35H"
    path.mkdir(parents=True,exist_ok=True)
    (path/f"b{a.degree}.json").write_text(json.dumps(result,indent=2)+"\n")
    print({k:v for k,v in result.items() if k not in ("legacy","relaxation")},flush=True)


if __name__=="__main__":
    main()
