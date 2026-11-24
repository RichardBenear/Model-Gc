"""GC35E: compare tidal readouts and audit the full complex TT projector."""
import argparse
import json
from pathlib import Path
import exgc35c as C
import exgc35e as E
import numpy as np


def projector(n,lattice):
    X=C.Experiment(n,lattice,True)
    S=X.system
    k=np.linalg.norm(X.K,axis=1)
    unit=X.K/k[:,None]
    plane=np.eye(3)-unit[:,:,None]*unit[:,None,:]
    def tt(v):
        tensor=np.empty((len(v),3,3),dtype=v.dtype)
        for j,(a,b) in enumerate(E.P.N.PAIRS):
            tensor[:,a,b]=tensor[:,b,a]=v[:,j]
        out=plane@tensor@plane
        return out-np.trace(out,axis1=-2,axis2=-1)[:,None,None]*plane/2
    weight=np.sqrt(np.tile(E.P.WT,2))
    matrix=np.concatenate([S.RC,S.RB*k[:,None,None]],axis=2)
    pinv=np.linalg.pinv(matrix/weight)/weight[:,None]
    tt_ops=np.stack([tt(pinv[:,start:start+6,i]) for start in (0,6) for i in range(4)],axis=1)
    error=np.linalg.norm(tt_ops.reshape(len(k),-1),axis=1)/np.linalg.norm(pinv.reshape(len(k),-1),axis=1)
    return tt,float(np.max(error))


def summarize(folder):
    refs={}
    arrays={}
    quadratics={}
    tags=("coarse","fine","continuum","half")
    for tag in tags:
        with (folder/f"{tag}.json").open() as file:
            refs[tag]=json.load(file)
        with np.load(folder/f"{tag}.npz") as data:
            arrays[tag]={k:data[k] for k in data if not k.endswith("_tt")}
        with np.load(folder/f"{tag}_quadratic.npz") as data:
            quadratics[tag]={k:data[k] for k in data}
    transforms={}
    audits={}
    for tag,meta in refs.items():
        params=meta["settings"]
        transforms[tag],audits[tag]=projector(params["n"],not params["continuum"])

    def rel(a,b,w=1):
        return float(np.linalg.norm((a-b)*w)/np.linalg.norm(b*w))

    rows=[]
    for t in E.READOUTS:
        continuum=arrays["continuum"]
        ctt=transforms["continuum"]
        ti=continuum[f"{t}_reconstructed"]
        for tag in ("coarse","fine","half"):
            if t>refs[tag]["settings"]["end"]:
                continue
            x=arrays[tag]
            d=x[f"{t}_reconstructed"]
            nominal=x[f"{t}_nominal"]
            weight=np.sqrt(E.P.WT)
            tt=transforms[tag]
            q=quadratics[tag][str(t)]
            qc=quadratics["continuum"][str(t)]
            rows.append(dict(t=t,tag=tag,reconstructed=rel(d,ti,weight),
                             coordinate_tidal=rel(d+q,ti+qc,weight),
                             coordinate_tidal_tt=rel(tt(d+q),ctt(ti+qc)),
                             nominal=rel(nominal,continuum[f"{t}_nominal"],weight),
                             reconstructed_tt=rel(tt(d),ctt(ti)),
                             nominal_tt=rel(tt(nominal),ctt(continuum[f"{t}_nominal"])),
                             free=rel(x[f"{t}_reconstructed_free"],continuum[f"{t}_reconstructed_free"],weight),
                             projected_to_free_tt=rel(tt(d),tt(x[f"{t}_reconstructed_free"])),
                             reconstructed_to_nominal=rel(d,nominal,weight),
                             impulse_ratio=float(np.linalg.norm(x[f"{t}_impulse"]*weight)/np.linalg.norm(ti*weight)),
                             continuum_projection_effect=rel(ti,continuum[f"{t}_reconstructed_free"],weight)))
    def peak(tag,key):
        return max(row[key] for row in refs[tag]["history"])
    coarse=refs["coarse"]["history"]
    half=refs["half"]["history"]
    t4_coarse=next(row for row in coarse if abs(row["t"]-4)<1e-10)
    t4_half=half[-1]
    fine_rows=[row for row in rows if row["tag"]=="fine"]
    gates=dict(post=max(peak(k,"after") for k in tags)<=1e-10,
               pre=peak("fine","before")<=1e-3 and peak("coarse","before")/peak("fine","before")>=2,
               tidal=all(row["reconstructed"]<=.02 for row in fine_rows))
    dt4=dict(coarse_before=t4_coarse["before"],half_before=t4_half["before"],
             reconstructed_relative=rel(arrays["half"]["4_reconstructed"],arrays["coarse"]["4_reconstructed"],np.sqrt(E.P.WT)),
             nominal_relative=rel(arrays["half"]["4_nominal"],arrays["coarse"]["4_nominal"],np.sqrt(E.P.WT)))
    return dict(tidal=rows,tt_projector_max_relative=audits,gates=gates,dt4=dt4,
                quadratic_cpu_seconds=float((folder/"quadratic_cpu.txt").read_text()))


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--output",default="EX-GC35E")
    args=parser.parse_args()
    path=Path(__file__).resolve().parents[1]/"results"/args.output
    result=summarize(path)
    (path/"tidal.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
