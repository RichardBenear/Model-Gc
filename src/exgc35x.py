"""GC35X background-gated probe in dilute exclusion cuts."""
from pathlib import Path
import json
import time
import numpy as np
import exgc35s as S


def measure(payload,n):
    f,b,ordered,D=S.layers(payload,n)
    lo,hi=D//4,3*D//4
    rows={}
    for name,ticks in (("forward",f),("reverse",b)):
        hits=0
        counts=[]
        for start in range(n):
            pos=start
            k=0
            for d in range(lo+1,hi+1):
                if d in ticks[pos]:
                    k+=1
                    pos=(pos+1)%n
            counts.append(k)
            hits+=k
        rows[name]=dict(depth_speed=hits/(n*(hi-lo)),proper_speed=1.0 if hits else None,
                        rates=S.slope(ticks,lo,hi).tolist(),from_sites=counts)
    own={}
    for (i,j),(d,parents,tag) in ordered:
        if lo<d<=hi:
            own[tag]=own.get(tag,0)+1
    rows["intrinsic_forward_depth_speed"]=sum(own.values())/(len(own)*(hi-lo))
    return rows


if __name__=="__main__":
    start=time.process_time()
    cases=[]
    for n in (32,64):
        for rho in (.125,.25):
            for style in ("balanced","shuffled"):
                bits=S.fixture(n,rho,style)
                rules=[measure(S.sim(bits,96,k)[0],n) for k in range(8)]
                cases.append(dict(n=n,rho=rho,style=style,exact=sum(x==rules[0] for x in rules),**rules[0]))
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35X"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(dict(cases=cases,cpu=time.process_time()-start),indent=2)+"\n")
    print([(x["n"],x["rho"],x["style"],x["exact"],x["forward"]["depth_speed"],x["reverse"]["depth_speed"],x["intrinsic_forward_depth_speed"]) for x in cases],flush=True)
