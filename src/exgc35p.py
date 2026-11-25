"""GC35P carrier catalogue; see the registered rule definitions."""
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import time
import networkx as nx
import exgc35l as L


def multisets(total, slots):
    if slots==1:
        yield (total,)
        return
    for first in range(total+1):
        for tail in multisets(total-first,slots-1):
            yield (first,*tail)


def moves(kind,bags):
    for i in range(1,len(bags)-1):
        q=bags[i]
        if kind in ("B","C") and len(q)>=2:
            options=list(dict.fromkeys(combinations(range(len(q)),2))) if kind=="B" else [tuple(range(len(q)))]
            for take in options:
                half=len(take)//2
                yield i,take,[(i-1,-1)]*half+[(i+1,1)]*(len(take)-half),"split"
        elif kind in ("A","D","E"):
            for d in (-1,1):
                ids=tuple(j for j,(_,v) in enumerate(q) if d==v)
                if kind=="D" and ids:
                    yield i,ids,[(i+d,d)]*len(ids),"queue"
                elif kind in ("A","E"):
                    for j in ids:
                        yield i,(j,),[(i+d,d)],"hop"
            if kind=="E":
                for j,k in combinations(range(len(q)),2):
                    if q[j][1]*q[k][1]==-1:
                        yield i,(j,k),[(i,0)]*2,"capture"


def advance(bags,m,e):
    i,take,out,_=m
    parents=tuple(bags[i][j] for j in take)
    q=[list(bag) for bag in bags]
    q[i]=[x for j,x in enumerate(q[i]) if j not in take]
    for j,d in out:
        q[j].append((e,d))
    return tuple(tuple(b) for b in q),parents


def initial(config,n,directional):
    if directional:
        return ((),)+tuple(tuple([(-i-1,-1)]*config[2*(i-1)]+[(-i-1,1)]*config[2*i-1]) for i in range(1,n+1))+((),)
    return ((),)+tuple(tuple([(-i-1,0)]*v) for i,v in enumerate(config,1))+((),)


def ledger(bags):
    return tuple(tuple(sum(v==d for _,v in b) for d in (-1,0,1)) for b in bags)


def terminal_graph(events,bags,labelled):
    g=nx.MultiDiGraph()
    def parents(values,node):
        for p,d in values:
            if p<0:
                g.add_node(p,label=("initial",-p-1 if labelled else None))
            g.add_edge(p,node,kind=d)
    for e,(i,kind,ps) in enumerate(events):
        g.add_node(e,label=(kind,i if labelled else None))
        parents(ps,e)
    for i,b in enumerate(bags):
        if b:
            e=len(events)+i
            g.add_node(e,label=("terminal",i if labelled else None))
            parents(b,e)
    return g


def graph_forks(kind,bags,needed):
    first=None
    failures={}
    paths=0
    def visit(q,ev):
        nonlocal first,paths
        for m in moves(kind,q):
            r,p=advance(q,m,len(ev))
            visit(r,ev+((m[0],m[3],p),))
            if needed<=failures.keys():
                return
        if any(moves(kind,q)):
            return
        paths+=1
        record=dict(events=ev,terminal=q,ledger=ledger(q))
        if first is None:
            first=record
            return
        for k in needed-failures.keys():
            lab=k=="labelled"
            a=terminal_graph(first["events"],first["terminal"],lab)
            b=terminal_graph(ev,q,lab)
            if not L.equivalent(a,b):
                failures[k]=dict(a=first,b=record)
    visit(bags,())
    return failures,paths


def ledger_solver(kind):
    @lru_cache(None)
    def endings(q):
        bags=tuple(tuple((0,d) for d in b) for b in q)
        targets=set()
        for m in moves(kind,bags):
            r,_=advance(bags,m,0)
            target=tuple(tuple(sorted(d for _,d in b)) for b in r)
            targets.add(target)
        if not targets:
            return frozenset([ledger(bags)])
        return frozenset(x for r in targets for x in endings(r))
    return endings


def walk(kind,q,reverse):
    ev=[]
    while True:
        options=list(moves(kind,q))
        if not options:
            return q,ev
        m=options[-1] if reverse else options[0]
        q,p=advance(q,m,len(ev))
        ev.append((m[0],m[3],p))


def run(kind):
    direction=kind in ("A","D","E")
    solve=ledger_solver(kind)
    count=failed=paths=checks=0
    witnesses={}
    conserved=True
    for total in range(8):
        for n in range(1,5):
            for cfg in multisets(total,n*(2 if direction else 1)):
                q=initial(cfg,n,direction)
                count+=1
                # Ordinary chip firing has a site-count ledger. Route tags on
                # split outputs are for graph edges, not distinct ledger species.
                qcounts=tuple(tuple(sorted(d for _,d in b)) for b in q)
                if kind in ("A","D"):
                    outcomes=[((sum(cfg[::2]),0,0),)+((0,0,0),)*n+((0,0,sum(cfg[1::2])),)]
                else:
                    outcomes=solve(qcounts)
                if kind in ("B","C"):
                    outcomes={tuple(sum(row) for row in r) for r in outcomes}
                    conserved &= all(sum(r)==total for r in outcomes)
                else:
                    conserved &= all(sum(map(sum,r))==total for r in outcomes)
                if len(outcomes)>1:
                    failed+=1
                    witnesses.setdefault("ledger",dict(start=cfg,n=n,terminals=sorted(outcomes)[:2]))
                if kind=="A":
                    a,ea=walk(kind,q,False)
                    b,eb=walk(kind,q,True)
                    assert ledger(a)==ledger(b)==outcomes[0]
                    # The exact scorer is exercised for both tie orders on the
                    # small configurations. Larger cases use the disjoint-path proof.
                    if total<=3:
                        assert L.equivalent(terminal_graph(ea,a,True),terminal_graph(eb,b,True))
                        assert L.equivalent(terminal_graph(ea,a,False),terminal_graph(eb,b,False))
                        checks+=1
                else:
                    needed={"labelled","unlabelled"}-witnesses.keys()
                    if needed:
                        found,p=graph_forks(kind,q,needed)
                        paths+=p
                        for k,v in found.items():
                            witnesses[k]=dict(start=cfg,n=n,**v)
    return dict(configurations=count,ledger_failures=failed,graph_paths_until_witness=paths,
                exact_no_overlap_checks=checks,conserved=bool(conserved),witnesses=witnesses)


if __name__=="__main__":
    start=time.process_time()
    records={k:run(k) for k in "ABCDE"}
    result=dict(candidates=records,cpu_seconds=time.process_time()-start)
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35P"
    p.mkdir(parents=True,exist_ok=True)
    (p/"catalogue.json").write_text(json.dumps(result,indent=2)+"\n")
    print({k:{a:b for a,b in v.items() if a!="witnesses"} for k,v in records.items()},flush=True)
