"""GC35S emergent-clock observables on exclusion-ring depth cuts."""
from bisect import bisect_right
import hashlib
import json
from pathlib import Path
import random
import time
import numpy as np
import exgc35l as L


def sim(bits,cut,seed,cap=None):
    n=len(bits)
    rng=random.Random(seed)
    state=[i if b else -1 for i,b in enumerate(bits)]
    roots=[None]*n
    serial=[0]*n
    events={}
    while cap is None or len(events)<cap:
        ready=[]
        for i in range(n):
            j=(i+1)%n
            if state[i]<0 or state[j]>=0:
                continue
            parents=tuple(sorted({x for x in (roots[i],roots[j]) if x is not None}))
            d=1+max((events[p][0] for p in parents),default=0)
            if d<=cut:
                ready.append((i,j,d,parents))
        if not ready:
            break
        i,j,d,p=rng.choice(ready)
        serial[i]+=1
        key=(i,serial[i])
        events[key]=(d,p,state[i])
        roots[i]=roots[j]=key
        state[i],state[j]=-1,state[i]
    payload=sorted(events.items())
    return payload,state,hashlib.sha256(repr((payload,state)).encode()).hexdigest()


def layers(payload,n):
    forward=[[] for _ in range(n)]
    back=[[] for _ in range(n)]
    maximum=max(d for _,(d,_,_) in payload)
    reverse={}
    ordered=sorted(payload,key=lambda row:(row[1][0],row[0]))
    for key,(d,parents,token) in reversed(ordered):
        r=reverse.get(key,1)
        reverse[key]=r
        for p in parents:
            reverse[p]=max(reverse.get(p,1),r+1)
    for key,(d,_,_) in payload:
        for i in (key[0],(key[0]+1)%n):
            forward[i].append(d)
            back[i].append(maximum+1-reverse[key])
    return forward,back,ordered,maximum


def slope(data,lo,hi):
    x=np.arange(lo,hi+1,dtype=float)
    z=x-x.mean()
    return np.array([np.dot(z,np.searchsorted(sorted(y),x,side="right"))/np.dot(z,z) for y in data])


def fixture(n,density,style):
    if style=="balanced":
        return [int((i+1)*density)-int(i*density) for i in range(n)]
    result=[1]*round(n*density)+[0]*(n-round(n*density))
    random.Random(351600+n+round(n*density)).shuffle(result)
    return result


def run(bits,cut):
    a,state,digest=sim(bits,cut,0)
    trials=[sim(bits,cut,k) for k in range(1,8)]
    f,b,ordered,depth=layers(a,len(bits))
    lo=max(1,depth//4)
    hi=max(lo+1,3*depth//4)
    r=slope(f,lo,hi)
    rr=slope(b,lo,hi)
    early=slope(f,lo,(hi+lo)//2)
    late=slope(f,(hi+lo)//2,hi)
    return dict(digest=digest,matched=1+sum(x==a and y==state and h==digest for x,y,h in trials),schedules=8,events=len(a),
                depth=depth,window=[lo,hi],forward=float(r.mean()),reverse=float(rr.mean()),
                mean_abs_site_foliation=float(np.mean(abs(r-rr))),
                half_window_spread=float(abs(early.mean()-late.mean()))),f,ordered,r


def step_case(n,uniform):
    bits=fixture(n//2,.25,"balanced")+fixture(n//2,.75,"balanced")
    cut=24 if n==64 else 12
    row,ticks,ordered,r=run(bits,cut)
    lo,hi=row["window"]
    w=np.array(bits)
    hist=[]
    i=0
    for d in range(1,hi+1):
        while i<len(ordered) and ordered[i][1][0]==d:
            site=ordered[i][0][0]
            w[site]=0
            w[(site+1)%n]=1
            i+=1
        if d>=lo:
            hist.append(w.copy())
    occ=np.mean(hist,axis=0)
    travels={}
    for key,(d,p,t) in ordered:
        travels.setdefault(t,[]).append((key[0],d))
    reports=[]
    fb=layers(ordered,n)[1]
    rev=slope(fb,lo,hi)
    calibration=[u for u in uniform if u["n"]==n and u["style"]=="balanced"]
    for start,stop,rho in ((n//8,3*n//8,.25),(5*n//8,7*n//8,.75)):
        ref,_,_,_=run(fixture(n,rho,"balanced"),cut)
        rate=float(np.mean(r[start:stop]))
        probes=[]
        for seq in travels.values():
            for j in range(len(seq)-4):
                segment=seq[j:j+5]
                if not all(start<=i<stop for i,_ in segment) or not (lo<=segment[0][1]<segment[-1][1]<=hi):
                    continue
                elapsed=segment[-1][1]-segment[0][1]
                read=0
                for (_,a),(site,b) in zip(segment,segment[1:]):
                    ts=sorted(ticks[site])
                    read+=bisect_right(ts,b)-bisect_right(ts,a)
                probes.append([elapsed,read])
        bg_speed=float(np.mean([4/(e*ref["forward"]) for e,t in probes])) if probes else None
        own_speed=float(np.mean([4/t for e,t in probes])) if probes and all(t for e,t in probes) else None
        observed=float(occ[start:stop].mean())
        predicted=float(np.interp(observed,[u["rho"] for u in calibration],[u["forward"] for u in calibration]))
        reports.append(dict(calibrated_rate=predicted,calibration_relative=rate/predicted-1,
                            reverse_rate=float(rev[start:stop].mean()),sites=[start,stop],initial_density=rho,occupancy=float(occ[start:stop].mean()),
                            rate=rate,uniform_rate=ref["forward"],ratio=rate/ref["forward"],
                            crossings=probes,own_speed=own_speed,background_speed=bg_speed))
    return dict(n=n,**row,plateaus=reports,occupancy=occ.tolist())


if __name__=="__main__":
    start=time.process_time()
    uniform=[]
    for n in (32,64):
        for rho in (.125,.25,.5,.75,.875):
            for style in ("balanced","shuffled"):
                row,*_=run(fixture(n,rho,style),96)
                uniform.append(dict(n=n,rho=rho,style=style,**row))
    steps=[step_case(n,uniform) for n in (32,64)]
    witness=L.case((2,2))[1]["token"]
    def terminal_depth(branch):
        depths=[]
        for site,parents in witness[branch]["events"]:
            depths.append(1+max((depths[p] for p in parents if p>=0),default=0))
        return depths
    control=[terminal_depth(x) for x in ("a","b")]
    bits=fixture(32,.5,"balanced")
    count_hashes=[sim(bits,96,k,97)[2] for k in range(8)]
    result=dict(uniform=uniform,steps=steps,paired_control_event_depths=control,fixed_count_distinct_hashes=len(set(count_hashes)),
                cpu_seconds=time.process_time()-start)
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35S"
    p.mkdir(parents=True,exist_ok=True)
    (p/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print([(x["n"],x["rho"],round(x["forward"],5),x["matched"]) for x in uniform],flush=True)
