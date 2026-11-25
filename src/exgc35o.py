"""GC35O: reciprocal scalar clocks and a contested-carrier fixture."""
from bisect import bisect_right
import hashlib
import json
from pathlib import Path
import random
import time
import exgc35l as L
import numpy as np

SIZE=48
H=.25
END=6.
COUPLING=.18
HOP=.34
TOKEN_WEIGHT=.16


def earlier(history,t):
    index=bisect_right(history,t,key=lambda row:row[0])-1
    stamp,state=history[index]
    if index==0 or stamp==t:
        return index,state
    previous,old=history[index-1]
    frac=(t-stamp)/(stamp-previous)
    return index,tuple(x+frac*(x-y) for x,y in zip(state,old))


def density(q,p,left,right):
    return .5*p*p+.25*((q-left)**2+(right-q)**2)/H**2


def static_check():
    x=np.arange(SIZE)*H
    k=2*np.pi/(SIZE*H)
    source=np.cos(k*x)
    K=2*np.sin(k*H/2)/H
    n=-COUPLING*source/K**2
    lap=(np.roll(n,-1)+np.roll(n,1)-2*n)/H**2
    return float(np.linalg.norm(lap-COUPLING*source)/np.linalg.norm(COUPLING*source))


def run(step,mode,seed=None,carrier=None):
    x=np.arange(SIZE)*H
    q=.36*np.exp(-((x-3)/1.4)**2)*np.cos(2*(x-3))
    p=-(np.roll(q,-1)-np.roll(q,1))/(2*H)
    rho=np.array([density(q[i],p[i],q[i-1],q[(i+1)%SIZE]) for i in range(SIZE)])
    mark=SIZE//2-2
    bags=((),(-2,-2),(-3,-3,-3),(-4,-4),())
    born=[[(0.,len(b))] for b in bags]
    if carrier:
        rho[mark:mark+5]+=TOKEN_WEIGHT*np.array(list(map(len,bags)))
    bg=float(rho.mean())
    phi=np.zeros(SIZE)
    if mode=="preset":
        k=2*np.sin(np.pi*np.fft.fftfreq(SIZE))/H
        k[0]=1.
        f=-COUPLING*np.fft.fft(rho-bg)/k**2
        f[0]=0
        phi=np.fft.ifft(f).real
    histories=[[(0.,(float(a),float(b),float(c),0.))] for a,b,c in zip(q,p,phi)]
    rng=random.Random(seed)
    terms=[]
    hops=[]
    clamps=0
    witness=None
    high_lapse=0.
    def clock(n):
        nonlocal clamps,high_lapse
        high_lapse=max(high_lapse,abs(n))
        if mode=="off":
            return 1.
        value=min(1.25,max(.75,1+n))
        if value!=1+n:
            clamps+=1
        return value
    nexts=[step/clock(s[0][1][2]) for s in histories]
    deadlines=[float("inf")]+[HOP/clock(histories[mark+i][0][1][2]) if carrier else float("inf")
                                     for i in range(1,4)]+[float("inf")]
    def hop(j,t):
        nonlocal bags
        pair=bags[j][:2] if carrier=="oldest" else bags[j][-2:]
        bags=L.update(bags,j,pair,len(hops))
        hops.append((j,pair,t))
        deadlines[j]=float("inf")
        for b in (j-1,j,j+1):
            born[b].append((t,len(bags[b])))
            if 0<b<4 and len(bags[b])>=2 and deadlines[b]==float("inf"):
                n=earlier(histories[mark+b],t)[1][2]
                deadlines[b]=t+HOP/clock(n)

    while min(min(nexts),min(deadlines))<=END:
        i=min(range(SIZE),key=nexts.__getitem__)
        j=min(range(5),key=deadlines.__getitem__)
        if deadlines[j]<=min(END,nexts[i]):
            hop(j,deadlines[j])
            continue
        if seed is not None:
            ready=[j for j in range(SIZE) if nexts[j]<=END and
                   min(nexts[j-1],nexts[(j+1)%SIZE])>histories[j][-1][0]]
            if not ready:
                raise RuntimeError("no locally ready event")
            i=rng.choice(ready)
        j,k=(i-1)%SIZE,(i+1)%SIZE
        due=nexts[i]
        t,(u,v,n,w)=histories[i][-1]
        a,(ul,_,nl,_)=earlier(histories[j],t)
        b,(ur,_,nr,_)=earlier(histories[k],t)
        count=0
        if carrier and mark<=i<mark+5:
            records=born[i-mark]
            count=records[bisect_right(records,t,key=lambda row:row[0])-1][1]
        nexts[i]=float("inf")
        dt=due-t
        rate=clock(n)
        rl,rr=clock(nl),clock(nr)
        area=3-2*rate
        speed=rate/np.sqrt(area)
        ml=(speed+rl/np.sqrt(3-2*rl))/2
        mr=(speed+rr/np.sqrt(3-2*rr))/2
        vv=v+dt*(mr*(ur-u)-ml*(u-ul))/H**2
        uu=u+dt*speed*vv
        src=density(u,v,ul,ur)/area+TOKEN_WEIGHT*count
        nw=w+dt*((nr-2*n+nl)/H**2-COUPLING*(src-bg))
        nn=n+dt*nw
        if mode=="preset":
            nn,nw=n,0.
        if mode=="source_off":
            nn,nw=0.,0.
        histories[i].append((due,(uu,vv,nn,nw)))
        nexts[i]=due+step/clock(nn)
        terms.append((i,len(histories[i])-1,j,a,k,b))
        if seed is not None and witness is None:
            lag=min(s[-1][0] for s in histories)
            if due-lag>4*step:
                witness=[i,due,lag]
    result=np.array([earlier(s,END)[1] for s in histories])
    versions=np.asarray([len(s)-1 for s in histories])
    digest=hashlib.sha256()
    for s in histories:
        digest.update(np.asarray([[t,*y] for t,y in s]).tobytes())
    digest.update(np.asarray(sorted(terms),dtype=np.int64).tobytes())
    assert np.isfinite(result).all() and sum(map(len,bags))==7
    stats=dict(step=step,mode=mode,carrier=carrier,seed=seed,digest=digest.hexdigest(),
               events=int(versions.sum()),event_range=[int(versions.min()),int(versions.max())],
               minmax_lapse=(1+np.array([result[:,2].min(),result[:,2].max()])).tolist(),
               max_lapse_perturbation=high_lapse,clamps=clamps,witness=witness,
               terminal_counts=list(map(len,bags)) if carrier else None,hops=hops)
    return result,stats,bags


def compare(a,b):
    return dict(phi_L2=float(np.linalg.norm(a[:,0]-b[:,0])),
                lapse_L2=float(np.linalg.norm(a[:,2]-b[:,2])))


def main():
    start=time.process_time()
    rhs={}
    for h in (.04,.02):
        ref,recs,_=run(h,"dynamic")
        trials=[]
        for seed in range(8):
            value,rec,_=run(h,"dynamic",3515+seed)
            trials.append(dict(record=rec,exact=value.tobytes()==ref.tobytes() and rec["digest"]==recs["digest"]))
        ctr={}
        for m in ("off","preset","source_off"):
            result,rec,_=run(h,m)
            ctr[m]=dict(record=rec,change=compare(ref,result))
        rhs[str(h)]=dict(reference=recs,schedules=trials,controls=ctr)
        if h==.04:
            coarse=ref
    x,a,ba=run(.02,"dynamic",carrier="oldest")
    y,b,bb=run(.02,"dynamic",carrier="newest")
    z,c,_=run(.02,"off",carrier="oldest")
    ga=L.graph(3,[(i,p) for i,p,_ in a["hops"]],ba,False)
    gb=L.graph(3,[(i,p) for i,p,_ in b["hops"]],bb,False)
    result=dict(static_relative=static_check(),packet=rhs,step_change=compare(coarse,ref),
                carrier=dict(oldest=a,newest=b,flat_clocks=c,same_fields=x.tobytes()==y.tobytes(),
                             same_token_dag=L.equivalent(ga,gb),change_flat=compare(x,z),token_source_effect=compare(x,ref)),
                cpu_seconds=time.process_time()-start)
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35O"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(dict(static=result["static_relative"],matches=[sum(x["exact"] for x in r["schedules"]) for r in rhs.values()],
               carrier_same_DAG=result["carrier"]["same_token_dag"],carrier_same_fields=result["carrier"]["same_fields"],
               cpu=result["cpu_seconds"]),flush=True)


if __name__=="__main__":
    main()
