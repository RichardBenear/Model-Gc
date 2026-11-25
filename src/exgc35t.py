"""GC35T persistent-vertex absorber and shared-counter control."""
from collections import deque
from functools import lru_cache
from itertools import product
import json
from pathlib import Path
import time
import networkx as nx
import exgc35l as L


def screen(lengths,serial):
    n=len(lengths)
    end=tuple(x+1 for x in lengths)
    initial=((0,)*n,())
    frontier=deque([initial])
    seen={initial}
    terminals=set()
    transitions=0
    commuting=0
    overlap=0
    for_case={}
    def outgoing(state):
        adv,order=state
        for i,k in enumerate(adv):
            if k<end[i]:
                nxt=adv[:i]+(k+1,)+adv[i+1:]
                fire=k==lengths[i]
                z=order+(i,) if serial and fire else order
                yield i,k,fire,(nxt,z)
    while frontier:
        s=frontier.popleft()
        children=list(outgoing(s))
        if not children:
            terminals.add(s[1])
        transitions+=len(children)
        for a,(i,k,f,nxt) in enumerate(children):
            for j,l,g,other in children[a+1:]:
                if serial and f and g:
                    overlap+=1
                    for_case.setdefault("minimal",dict(state=s,competing=[i,j]))
                else:
                    left=next(z for m,_,_,z in outgoing(nxt) if m==j)
                    right=next(z for m,_,_,z in outgoing(other) if m==i)
                    assert left==right
                    commuting+=1
            if nxt not in seen:
                frontier.append(nxt)
                seen.add(nxt)
    @lru_cache(None)
    def paths(adv):
        if adv==end:
            return 1
        return sum(paths(adv[:i]+(v+1,)+adv[i+1:]) for i,v in enumerate(adv) if v<end[i])
    def graph(order,anchored):
        g=nx.MultiDiGraph()
        for i,l in enumerate(lengths):
            for k in range(l+1):
                node=("event",i,k)
                g.add_node(node,label="absorb" if k==l else "hop")
                if k:
                    g.add_edge(("event",i,k-1),node,kind="carrier")
                elif anchored:
                    src=("initial",i)
                    g.add_node(src,label=src)
                    g.add_edge(src,node,kind="carrier")
            if anchored:
                dest=("mass",i)
                g.add_node(dest,label=dest)
                g.add_edge(("event",i,l),dest,kind="mass")
        if serial:
            for i,j in zip(order,order[1:]):
                g.add_edge(("event",i,lengths[i]),("event",j,lengths[j]),kind="counter")
            if anchored:
                for name,i in (("initial-counter",order[0]),("final-counter",order[-1])):
                    name=(name,)
                    g.add_node(name,label=name)
                    e=("event",i,lengths[i])
                    pair=(name,e) if name=="initial-counter" else (e,name)
                    g.add_edge(*pair,kind="counter")
        return g
    orders=sorted(terminals)
    bad={}
    for label in ("plain","bordered"):
        ref=graph(orders[0],label=="bordered")
        for order in orders[1:]:
            if not L.equivalent(ref,graph(order,label=="bordered")):
                bad[label]=[orders[0],order]
                break
    return dict(lengths=lengths,reachable=len(seen),transitions=transitions,represented_histories=paths((0,)*n),
                independent_pairs=commuting,destructive_pairs=overlap,witness=for_case,
                terminal_orders=len(orders),failures=bad,conserved_total=n)


if __name__=="__main__":
    start=time.process_time()
    results={}
    for serial in (False,True):
        rows=[screen(t,serial) for n in range(1,5) for t in product(range(1,4),repeat=n)]
        # Also check the minimal direct absorber with no prior route hops.
        direct=[screen((0,)*n,serial) for n in range(1,5)]
        results["counter" if serial else "persistent"]=dict(direct=direct,routed=rows)
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35T"
    p.mkdir(parents=True,exist_ok=True)
    (p/"results.json").write_text(json.dumps(dict(cases=results,cpu_seconds=time.process_time()-start),indent=2)+"\n")
    print({k:{m:(len(v),sum(bool(x["failures"]) for x in v)) for m,v in b.items()} for k,b in results.items()},flush=True)
