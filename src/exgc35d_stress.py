"""Independent quadratic/quartic stress and plaquette-placement checks."""
import json
from pathlib import Path
import exgc35d as D
import numpy as np
import torch


def audit():
    rng=np.random.default_rng(35041)
    D2err,D4err,placement=0.,0.,0.
    for _ in range(20):
        f=rng.standard_normal((3,4))
        gram=f@f.T
        tr=np.trace(gram)
        quartic=(tr**2-np.trace(gram@gram))/2
        t2=2*gram-tr*np.eye(3)
        t4=2*(tr*gram-gram@gram)-quartic*np.eye(3)
        metric=torch.eye(3,requires_grad=True)
        inv=torch.linalg.inv(metric)
        volume=torch.sqrt(torch.linalg.det(metric))
        contracted=inv@torch.tensor(gram)
        q2=torch.trace(contracted)
        q4=(q2.square()-torch.trace(contracted@contracted))/2
        g2=torch.autograd.grad(volume*q2,metric,retain_graph=True)[0].numpy()
        g4=torch.autograd.grad(volume*q4,metric)[0].numpy()
        D2err=max(D2err,float(np.linalg.norm(2*g2+t2)/np.linalg.norm(t2)))
        D4err=max(D4err,float(np.linalg.norm(2*g4+t4)/np.linalg.norm(t4)))
        n=17
        mode=rng.integers(-4,5,size=3)
        site=np.exp(2j*np.pi*np.einsum("a,aijk->ijk",mode,np.indices((n,n,n)))/n)
        coeff=2*np.pi*mode/n
        for a,b in ((0,1),(1,2),(0,2)):
            placed=(site+np.roll(site,-1,a)+np.roll(site,-1,b)+np.roll(np.roll(site,-1,a),-1,b))/4
            measured=placed*np.exp(-.5j*(coeff[a]+coeff[b]))
            expected=site*np.cos(coeff[a]/2)*np.cos(coeff[b]/2)
            placement=max(placement,float(np.max(abs(measured-expected))))
    return dict(sigma_stress_relative=D2err,quartic_stress_relative=D4err,plaquette_phase_max=placement)


if __name__=="__main__":
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35D"
    path.mkdir(parents=True,exist_ok=True)
    result=audit()
    (path/"stress_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
