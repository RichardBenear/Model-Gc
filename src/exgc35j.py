"""GC35J: discrete chain rule and staggered quartic test."""
import json
from pathlib import Path
import time
import exgc35h as H
import numpy as np
import torch


def audit(q,v,h,mu,soft):
    s=H.dot(q,q)
    if soft:
        potential=mu*(s-1)**2/4
        force=mu*(s-1)[...,None]*q
    else:
        potential=mu*s/2
        force=mu*q
    lap=H.P.N.lap(q,h)
    accel=lap-force
    diffs=[]
    predicted=[]
    scale=[]
    contraction=[]
    for a in range(3):
        nxt=lambda x:H.fw(x,a)
        chord=(nxt(q)-q)/h
        dv=(nxt(v)-v)/h
        pdot=lambda acc:-.5*(H.dot(acc+nxt(acc),chord)+H.dot(v+nxt(v),dv))
        dvsite=(nxt(potential)-potential)/h
        diffs.append(pdot(accel)-pdot(lap)-dvsite)
        predicted.append(.5*H.dot(force+nxt(force),chord)-dvsite)
        scale.append(np.linalg.norm(dvsite))
        contraction.append(np.linalg.norm(.5*H.dot(force+nxt(force),chord)))
    epsilon=1e-5
    rp=H.dot(q+epsilon*v,q+epsilon*v)
    rm=H.dot(q-epsilon*v,q-epsilon*v)
    diff=mu*((rp-1)**2-(rm-1)**2)/(8*epsilon) if soft else mu*(rp-rm)/(4*epsilon)
    exact=H.dot(force,v)
    return dict(mu=mu,soft=soft,energy_chain_fd=float(np.linalg.norm(diff-exact)/max(np.linalg.norm(exact),1e-100)),
                momentum_max=float(np.max(np.abs(diffs))),momentum_rel=float(np.linalg.norm(diffs)/max(np.linalg.norm(scale),np.linalg.norm(contraction),1e-100)),
                predicted_difference=float(np.max(abs(np.asarray(diffs)-predicted))))


def quartic(q,h):
    t=torch.tensor(q,requires_grad=True)
    energy=0.
    for orient in H.B.ORIENT:
        ds=torch.stack([(torch.roll(t,-s,a)-t)*(s/h) for a,s in enumerate(orient)],dim=-2)
        gram=ds@ds.transpose(-1,-2)
        tr=gram.diagonal(dim1=-2,dim2=-1).sum(-1)
        energy+=(tr.square()-(gram*gram).sum((-1,-2))).sum()/(16*H.NORM)
    g=torch.autograd.grad(energy,t)[0].detach().numpy()
    measured=H.divergence(H.components(q,h)[3],h)
    gp,_=H.gradients(q,h)
    expected=np.stack([-.5*H.dot(g+H.fw(g,a),gp[a]) for a in range(3)],axis=-1)
    return dict(absolute=float(np.linalg.norm(measured-expected)),
                relative=float(np.linalg.norm(measured-expected)/max(np.linalg.norm(expected),1e-100)),
                measured=float(np.linalg.norm(measured)))


def main():
    start=time.process_time()
    rng=np.random.default_rng(3510)
    q=rng.normal(size=(7,7,7,4))*.5
    v=rng.normal(size=q.shape)*.4
    p=np.zeros_like(q)
    p[...,0]=np.array([1,-1,1,1,1,1,1])[:,None,None]
    rows=[]
    for h in (.4,.2):
        rows.append(dict(h=h,quadratic=audit(q,v,h,4,False),
                         soft=[audit(q,v,h,mu,True) for mu in (1.,4.,16.)],
                         rank_one=audit(p,0*p,h,4,True),quartic=quartic(q,h),rank_one_quartic=quartic(p,h)))
    result=dict(rows=rows,gate=all(r["momentum_rel"]<=1e-12 for x in rows for r in x["soft"]),
                cpu_seconds=time.process_time()-start)
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35J"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)


if __name__=="__main__":
    main()
