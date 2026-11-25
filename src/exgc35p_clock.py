"""GC35P: non-overlapping directed carrier paths in the GC35O clock loop."""
from bisect import bisect_right
import hashlib
import json
from pathlib import Path
import random
import time
import exgc35o as O
import numpy as np

PATHS=((22,23,24,25,26),(22,23,24,25,26),(25,24,23,22,21))


def run(step,seed,feedback=True):
    size=O.SIZE
    x=np.arange(size)*O.H
    q=.36*np.exp(-((x-3)/1.4)**2)*np.cos(2*(x-3))
    p=-(np.roll(q,-1)-np.roll(q,1))/(2*O.H)
    rho=np.array([O.density(q[i],p[i],q[i-1],q[(i+1)%size]) for i in range(size)])
    background=float((rho.sum()+len(PATHS)*O.TOKEN_WEIGHT)/size)
    fields=[[(0.,(float(u),float(v),0.,0.))] for u,v in zip(q,p)]
    token_times=[[0.] for _ in PATHS]
    indices=[{site:i for i,site in enumerate(path)} for path in PATHS]
    crossing=[[j for j,idx in enumerate(indices) if i in idx] for i in range(size)]
    f_due=[step]*size
    c_due=[O.HOP]*len(PATHS)
    rng=random.Random(seed)
    records={}
    witness=None
    clamps=0

    def lapse(n):
        nonlocal clamps
        if not feedback:
            return 1.
        val=min(1.25,max(.75,1+n))
        clamps+=val!=1+n
        return val

    def no_unseen_carrier(i,t):
        for j in crossing[i]:
            distance=indices[j][i]-(len(token_times[j])-1)
            if distance>=0 and c_due[j]+max(0,distance-1)*O.HOP/1.25<=t:
                return False
        return True

    def field_event(i):
        j,k=(i-1)%size,(i+1)%size
        t,(u,v,n,w)=fields[i][-1]
        when=f_due[i]
        d=when-t
        il,(ul,_,nl,_)=O.earlier(fields[j],t)
        ir,(ur,_,nr,_)=O.earlier(fields[k],t)
        count=0
        deps=[]
        for c in crossing[i]:
            version=bisect_right(token_times[c],t)-1
            target=indices[c][i]
            count+=version==target
            if version>=target:
                deps.append((c,min(version,target+1)))
        N=lapse(n)
        area=3-2*N
        speed=N/np.sqrt(area)
        left,right=lapse(nl),lapse(nr)
        l=.5*(speed+left/np.sqrt(3-2*left))
        r=.5*(speed+right/np.sqrt(3-2*right))
        vv=v+d*(r*(ur-u)-l*(u-ul))/O.H**2
        ww=w+d*((nr-2*n+nl)/O.H**2-O.COUPLING*(
            O.density(u,v,ul,ur)/area+O.TOKEN_WEIGHT*count-background))
        nn=n+d*ww
        fields[i].append((when,(u+d*speed*vv,vv,nn,ww)))
        f_due[i]=when+step/lapse(nn)
        records[(0,i,len(fields[i])-1)]=(j,il,k,ir,tuple(deps))

    def carrier_event(j):
        when=c_due[j]
        stage=len(token_times[j])
        dest=PATHS[j][stage]
        ref,state=O.earlier(fields[dest],when)
        token_times[j].append(when)
        c_due[j]=when+O.HOP/lapse(state[2]) if stage+1<len(PATHS[j]) else float("inf")
        records[(1,j,stage)]=(dest,ref)

    while min(min(f_due),min(c_due))<=O.END:
        ready=[]
        for i,when in enumerate(f_due):
            t=fields[i][-1][0]
            if when<=O.END and min(f_due[i-1],f_due[(i+1)%size])>t and no_unseen_carrier(i,t):
                ready.append((when,0,i))
        for j,when in enumerate(c_due):
            if when<=O.END and f_due[PATHS[j][len(token_times[j])]]>when:
                ready.append((when,1,j))
        if not ready:
            raise RuntimeError("deadline deadlock")
        when,kind,index=min(ready) if seed is None else rng.choice(ready)
        if kind==0:
            field_event(index)
        else:
            carrier_event(index)
        if witness is None and when-min(h[-1][0] for h in fields)>4*step:
            witness=[when,min(h[-1][0] for h in fields)]
    data=np.array([O.earlier(a,O.END)[1] for a in fields])
    digest=hashlib.sha256()
    for arr in fields:
        digest.update(np.asarray([[t,*x] for t,x in arr]).tobytes())
    for arr in token_times:
        digest.update(np.asarray(arr).tobytes())
    digest.update(str(sorted(records.items())).encode())
    events=[len(a)-1 for a in fields]
    assert np.isfinite(data).all()
    return data,dict(digest=digest.hexdigest(),events=sum(events),site_event_range=[min(events),max(events)],
                     hops=token_times,terminal_counts=[sum(path[len(t)-1]==q for path,t in zip(PATHS,token_times))
                                                      for q in (21,26)],
                     witness=witness,clamps=int(clamps),lapse_range=(1+np.array([data[:,2].min(),data[:,2].max()])).tolist())


if __name__=="__main__":
    start=time.process_time()
    runs=[]
    for step in (.04,.02):
        ref,a=run(step,None)
        flat,b=run(step,None,False)
        trials=[]
        for seed in range(8):
            out,row=run(step,3520+seed)
            trials.append(dict(equal=(out.tobytes()==ref.tobytes() and row["digest"]==a["digest"]),**row))
        runs.append(dict(step=step,reference=a,flat=b,schedules=trials,feedback=O.compare(ref,flat)))
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35P"
    p.mkdir(parents=True,exist_ok=True)
    result=dict(runs=runs,cpu_seconds=time.process_time()-start)
    (p/"reciprocal.json").write_text(json.dumps(result,indent=2)+"\n")
    print([(r["step"],sum(x["equal"] for x in r["schedules"]),r["reference"]) for r in runs],flush=True)
