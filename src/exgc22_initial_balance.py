"""Post-hoc t=0 instantaneous weighted-velocity balance from settled checkpoints.

No dynamics run, and no fit. At an initially symmetric resting state the recoil
and cohesion integrals should be zero even after a uniform boost. The measured
flow velocity is a consumption-weighted average whose derivative also contains
a reweighting term. This term need not vanish if the internal velocity is radial.
"""
import json
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
import ex008_longrange as e
import exgc09_lattice as g9
import exgc20_wave as wave

root=Path(e.ROOT); out=[]
for p in sorted((root/'local_raw/EX-GC22').glob('settle_*.npz')):
    rec=json.loads((root/'results/EX-GC22'/p.with_suffix('.json').name).read_text()); K,D,h=(rec[k] for k in ('K','D','h'))
    with np.load(p) as f: u=f['u']; V=f['V']
    V[2]+=.05; X,Y,W,G=u; s=K*G*W; ss=s.sum(); vc=(s*V[2]).sum()/ss
    bc=[np.ones(G.shape[1:]) for _ in range(6)]
    J=[-D*g for g in wave.grad_faces(G,h,bc)]
    gd=-wave.div_faces(J,h)-s
    wd=e.DW*g9.lap(W,h)-W+e.C*Y*W*W-W**3-sum(np.gradient(W*V[i],h,axis=i) for i in range(3))
    sd=K*(G*wd+W*gd)
    recoil= (s*(s/(s+1e-3*s.max())) *2*wave.centre_flux(J)[2]/G).sum()/ss/20
    sm=gaussian_filter(s,1.5/h,mode='constant')+1e-12
    cohesion=-(s*(V[2]-gaussian_filter(s*V[2],1.5/h,mode='constant')/sm)).sum()/ss
    reweight=(sd*(V[2]-vc)).sum()/ss
    # Diagnostic only: uniform convective transport of the pre-existing local V
    # would contribute -v0 ∂z V. This is not a proposed production rule or proof
    # of momentum conservation; check whether it cancels the initial reweighting.
    transport=-(s*vc*np.gradient(V[2],h,axis=2)).sum()/ss
    out.append(dict(K=K,D=D,h=h,initial_V=float(vc),recoil=float(recoil),cohesion=float(cohesion),reweight=float(reweight),
                    initial_gamma=float(-(recoil+cohesion+reweight)/vc),uniform_transport_term=float(transport),
                    initial_gamma_with_uniform_transport=float(-(recoil+cohesion+reweight+transport)/vc),
                    G_equilibrium_max_residual=float(abs(gd).max())))
    print(out[-1])
(root/'results/EX-GC22/initial_balance.json').write_text(json.dumps(out,indent=2)+'\n')
