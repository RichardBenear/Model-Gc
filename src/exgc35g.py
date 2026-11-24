"""GC35G: compare the h=.125 run to the original GC35E comparison data."""
from pathlib import Path
import argparse
import json
import time
import exgc35e as E
from exgc35e_analysis import projector
from exgc35e_curvature import correction
import numpy as np


def relative(a,b,w=1):
    return float(np.linalg.norm((a-b)*w)/np.linalg.norm(b*w))


def compare(new,reference,old):
    run=json.loads((new/"finer.json").read_text())
    prior=json.loads((reference/"tidal.json").read_text())
    prev=json.loads((reference/"fine.json").read_text())
    with np.load(new/"finer.npz") as p, np.load(reference/"continuum.npz") as r, np.load(reference/"continuum_quadratic.npz") as q:
        start=time.process_time()
        geom=correction(E.Projected(96,True),.05,20)
        tt,bound=projector(96,True)
        ctt,_=projector(24,False)
        rows=[]
        w=np.sqrt(E.P.WT)
        for t in E.READOUTS:
            s=str(t)
            a,b=p[f"{t}_reconstructed"],r[f"{t}_reconstructed"]
            tidal=relative(a,b,w)
            tensor=relative(tt(a),ctt(b))
            baseline=next(x for x in prior["tidal"] if x["t"]==t and x["tag"]=="fine")
            rows.append(dict(t=t,tidal=tidal,TT=tensor,
                             coordinate=relative(a+geom[s],b+q[s],w),
                             nominal=relative(p[f"{t}_nominal"],r[f"{t}_nominal"],w),
                             free=relative(p[f"{t}_reconstructed_free"],r[f"{t}_reconstructed_free"],w),
                             impulse=float(np.linalg.norm(w*p[f"{t}_impulse"])/np.linalg.norm(w*b)),
                             projected_to_free_tt=relative(tt(a),tt(p[f"{t}_reconstructed_free"])),
                             h25_ratio=baseline["reconstructed"]/tidal,h25_tt_ratio=baseline["reconstructed_tt"]/tensor))
    peak=lambda data,key:max(x[key] for x in data["history"])
    checkpoints=[next(x for x in run["history"] if x["t"]==t) for t in E.READOUTS]
    before,after=peak(run,"before"),peak(run,"after")
    late=abs(rows[-1]["tidal"]/rows[-2]["tidal"]-1)
    original=json.loads((old/"fine.json").read_text())
    reproduction=max(abs(next(x for x in prev["history"] if x["t"]==t)["free"]/
                        next(x for x in original["history"] if x["t"]==t)["next_scaled"]-1) for t in (4,20))
    return dict(rows=rows,checkpoints=checkpoints,before=before,after=after,free=peak(run,"free"),
                cumulative=run["cumulative"],before_ratio=peak(prev,"before")/before,tt_projector=bound,
                production_cpu=run["cpu_seconds"],analysis_cpu=time.process_time()-start,
                original_free_replay_relative=reproduction,late_change=late,
                gates=dict(pre=before<=1e-3,post=after<=1e-10,
                           tidal=all(x["tidal"]<=.02 for x in rows)))


if __name__=="__main__":
    root=Path(__file__).resolve().parents[1]/"results"
    p=argparse.ArgumentParser()
    p.add_argument("--reference",default="EX-GC35E")
    a=p.parse_args()
    result=compare(root/"EX-GC35G",root/a.reference,root/"EX-GC35C")
    (root/"EX-GC35G"/"comparison.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
