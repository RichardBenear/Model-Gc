"""GC35W recoil with a persistent particle and a versioned position relation."""
from itertools import combinations, combinations_with_replacement, product
import json
from pathlib import Path
import time
import networkx as nx
import exgc35l as L

SL=5


def moves(state,positions,threshold):
    pos,pid,waiting,free,marks=state
    for j in waiting:
        new=(pos,pid,tuple(i for i in waiting if i!=j),free+((j,-1),),marks)
        yield ("arrive",j),new,()
    for j,parent in free:
        delta=pos-positions[j]
        if abs(delta)==1:
            yield ("absorb",j,delta),(pos,pid,waiting,tuple(x for x in free if x[0]!=j),marks),tuple(x for x in (pid,parent) if x>=0)
    net=sum(s for _,s,_ in marks)
    direction=(net>0)-(net<0)
    if abs(net)>=threshold and 0<=pos+direction<SL:
        for chosen in combinations([m for m in marks if m[1]==direction],threshold):
            yield ("recoil",direction,tuple(j for j,_,_ in chosen)),(
                pos+direction,pid,waiting,free,tuple(m for m in marks if m not in chosen)),tuple(x for x in (pid,*(m[2] for m in chosen)) if x>=0)


def successor(op,new,e):
    pos,pid,wait,free,marks=new
    if op[0]=="arrive":
        free=free[:-1]+((op[1],e),)
    if op[0]=="absorb":
        marks=marks+((op[1],op[2],e),)
    if op[0]=="recoil":
        pid=e
    return (pos,pid,wait,free,marks)


def graph(state,ev,anchored):
    g=nx.MultiDiGraph()
    for i,(op,deps) in enumerate(ev):
        e=("e",i)
        g.add_node(e,label=op if anchored else op[0])
        for p in set(deps):
            g.add_edge(("e",p),e,kind="dependency")
    if anchored:
        pos,pid,wait,free,marks=state
        for who,p in [(("position",pos),pid),*((("free",j),p) for j,p in free),*((("deposit",j,s),p) for j,s,p in marks)]:
            g.add_node(who,label=who)
            if p>=0:
                g.add_edge(("e",p),who,kind="output")
    return g


def screen(positions,delayed,threshold):
    start=(2,-1,tuple(i for i,v in enumerate(delayed) if v),tuple((i,-1) for i,v in enumerate(delayed) if not v),())
    ref=None
    failures={}
    overlaps=paths=0
    def walk(state,ev):
        nonlocal ref,overlaps,paths
        options=list(moves(state,positions,threshold))
        if options:
            overlaps+=sum(a[0][0]=="recoil" and b[0][0]=="absorb" and abs(a[1][0]-positions[b[0][1]])!=1
                          for a in options for b in options)
            for op,st,p in options:
                walk(successor(op,st,len(ev)),ev+((op,p),))
            return
        paths+=1
        pos,pid,wait,free,marks=state
        ledger=(pos,tuple(sorted(j for j,p in free)),tuple(sorted((j,s) for j,s,p in marks)))
        absorbed=sum(e[0][0]=="absorb" for e in ev)
        consumed=threshold*sum(e[0][0]=="recoil" for e in ev)
        assert len(free)+len(wait)+absorbed==len(positions) and consumed+len(marks)==absorbed
        aggregate=(pos,tuple(sorted(positions[j] for j,p in free)),tuple(sorted(s for j,s,p in marks)))
        x=dict(ledger=ledger,aggregate=aggregate,events=ev)
        if ref is None:
            ref=(x,graph(state,ev,False),graph(state,ev,True))
            return
        if ledger!=ref[0]["ledger"]:
            failures.setdefault("ledger",[ref[0],x])
        if aggregate!=ref[0]["aggregate"]:
            failures.setdefault("aggregate",[ref[0],x])
        for mode,key in (("plain",1),("anchored",2)):
            if mode not in failures and not L.equivalent(ref[key],graph(state,ev,key==2)):
                failures[mode]=[ref[0],x]
    walk(start,())
    return dict(paths=paths,overlaps=overlaps,failures=failures)


if __name__=="__main__":
    t=time.process_time()
    output={}
    for h in (1,2):
        counts=dict(configs=0,paths=0,overlap_configs=0)
        witness={}
        overlap_witness={}
        for n in range(1,5):
            for pos in combinations_with_replacement(range(SL),n):
                for delays in product((0,1),repeat=n):
                    s=screen(pos,delays,h)
                    counts["configs"]+=1
                    counts["paths"]+=s["paths"]
                    counts["overlap_configs"]+=bool(s["overlaps"])
                    for k,v in s["failures"].items():
                        counts[k]=counts.get(k,0)+1
                        example=dict(positions=pos,delayed=delays,histories=v)
                        witness.setdefault(k,example)
                        if s["overlaps"]:
                            overlap_witness.setdefault(k,example)
        output[h]=dict(counts=counts,witnesses=witness,overlap_witnesses=overlap_witness)
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35W"
    p.mkdir(parents=True,exist_ok=True)
    (p/"results.json").write_text(json.dumps(dict(data=output,cpu=time.process_time()-t),indent=2)+"\n")
    print({h:v["counts"] for h,v in output.items()},flush=True)
