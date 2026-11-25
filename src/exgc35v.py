"""GC35V bounded typed directed binary hypergraph critical pairs."""
from collections import Counter, defaultdict
from itertools import combinations, combinations_with_replacement as multisets, product
import json
from pathlib import Path
import time
import networkx as nx
import exgc35l as L

CAP=700


def choices(rule,edges,step=0):
    left,right=rule
    lookup=defaultdict(list)
    for i,(_,edge) in enumerate(edges):
        lookup[edge].append(i)
    vertices=set(v for _,(_,a,b) in edges for v in (a,b))
    ids=sorted(set(v for _,a,b in left for v in (a,b)))
    fresh=sorted(set(v for _,a,b in right for v in (a,b))-set(ids))
    found=set()
    for values in product(sorted(vertices),repeat=len(ids)):
        bind=dict(zip(ids,values))
        bind.update({v:max(max(vertices),100*(step+1))+j+1 for j,v in enumerate(fresh)})
        buckets=[lookup[t,bind[a],bind[b]] for t,a,b in left]
        for take in product(*buckets):
            if len(set(take))!=len(left):
                continue
            emit=tuple(sorted((t,bind[a],bind[b]) for t,a,b in right))
            found.add((tuple(sorted(take)),emit))
    return sorted(found)


def apply(state,match):
    edges,history=state
    take,out=match
    return tuple(x for i,x in enumerate(edges) if i not in take)+tuple((len(history),e) for e in out),history+(tuple(edges[j][0] for j in take),)


def graph(state,seed,types,anchor=True):
    edges,ev=state
    g=nx.MultiDiGraph()
    boundaries=set(v for _,a,b in seed for v in (a,b))
    for i,deps in enumerate(ev):
        g.add_node(("e",i),label="event")
        for p in deps:
            if p>=0:
                g.add_edge(("e",p),("e",i),kind="depend")
            elif anchor and seed[-p-1][0] in types:
                root=("i",-p-1)
                g.add_node(root,label=root)
                g.add_edge(root,("e",i),kind="depend")
    for j,(p,(t,a,b)) in enumerate(edges):
        e=("f",j)
        g.add_node(e,label=("edge",t))
        for k,v in (("tail",a),("head",b)):
            node=("v",v)
            g.add_node(node,label=node if v in boundaries else "fresh")
            g.add_edge(e,node,kind=k)
        if p>=0:
            g.add_edge(("e",p),e,kind="produce")
        elif anchor and seed[-p-1][0] in types:
            root=("i",-p-1)
            g.add_node(root,label=root)
            g.add_edge(root,e,kind="produce")
    return g


def join(rule,state,a,b,depth=3,anchor=True,root=None):
    seed=tuple(e for _,e in state[0]) if root is None else root
    types={x[0] for x in rule[1]}
    left=[apply(state,a)]
    right=[apply(state,b)]
    for level in range(1,depth+1):
        table=defaultdict(list)
        for x in left:
            shape=tuple(sorted(Counter(i[0] for _,i in x[0]).items()))
            table[shape].append(x)
        cache={}
        for y in right:
            shape=tuple(sorted(Counter(i[0] for _,i in y[0]).items()))
            if shape not in table:
                continue
            gy=graph(y,seed,types,anchor)
            if shape not in cache:
                cache[shape]=[graph(x,seed,types,anchor) for x in table[shape]]
            if any(L.equivalent(g,gy) for g in cache[shape]):
                return "join",level
        if level<depth:
            left=[apply(x,c) for x in left for c in choices(rule,x[0],len(x[1]))]
            right=[apply(x,c) for x in right for c in choices(rule,x[0],len(x[1]))]
            if max(len(left),len(right))>CAP:
                return "cap",level
    return "unresolved",depth


def starts(n):
    atlas=list(product(("C","F"),range(2),range(2)))
    for q in multisets(atlas,n):
        if any(t=="C" for t,_,_ in q):
            yield q


def make(seed):
    return (tuple((-i-1,e) for i,e in enumerate(seed)),())


def check(rule,maximum,minimal,reached=False):
    tested=0
    overlap=0
    first=None
    for n in range(1,maximum+1):
        for seed in starts(n):
            state=make(seed)
            groups=[(state,choices(rule,state[0]),3)]
            if reached:
                groups.extend((s,choices(rule,s[0]),2) for m in groups[0][1] for s in [apply(state,m)])
            for current,moves,budget in groups:
                for a,b in combinations(moves,2):
                    if minimal and (not set(a[0])&set(b[0]) or len(set(a[0])|set(b[0]))!=n):
                        continue
                    overlap+=bool(set(a[0])&set(b[0]))
                    tested+=1
                    status,h=join(rule,current,a,b,budget,root=seed)
                    if status!="join":
                        terminal_pair=all(not choices(rule,x[0],len(x[1])) for x in (apply(current,a),apply(current,b)))
                        first=dict(seed=seed,state=current,choices=[a,b],status=status,depth=h,terminal_pair=terminal_pair)
                        return dict(tested=tested,overlap=overlap,first=first)
    return dict(tested=tested,overlap=overlap,first=None)


if __name__=="__main__":
    t=time.process_time()
    initial=[
      (("C",a,b),) for a,b in ((0,0),(0,1))
    ]+[(("C",0,1),("F",a,b)) for a,b in product(range(2),repeat=2)]
    rows=[]
    for lhs in initial:
        vertices=sorted(set(v for t,a,b in lhs for v in (a,b)))
        pairs=list(product(vertices,repeat=2))
        for a,b in pairs:
            for size in range(3):
                for fs in multisets(pairs,size):
                    rhs=(("C",a,b),)+tuple(("F",i,j) for i,j in fs)
                    if Counter(lhs)==Counter(rhs):
                        continue
                    rule=lhs,rhs
                    first=check(rule,3,True)
                    finite=check(rule,4,False,True) if first["first"] is None and first["overlap"] else None
                    rows.append(dict(rule=rule,critical=first,finite=finite))
    controls=[]
    published=[
        ((("C",0,1),),(("C",0,1),("C",1,2))),
        ((("C",0,1),("C",2,1)),(("C",0,3),("C",1,3),("C",2,3)))
    ]
    for rule,seed in zip(published,[(("C",0,0),("C",0,0)),(("C",0,3),("C",1,3),("C",2,3))]):
        state=make(seed)
        pairs=[]
        for a,b in combinations(choices(rule,state[0]),2):
            pairs.append(dict(shared=bool(set(a[0])&set(b[0])),border=join(rule,state,a,b),
                              output_only=join(rule,state,a,b,anchor=False)))
        controls.append(dict(rule=rule,seed=seed,pairs=pairs))
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35V"
    out.mkdir(parents=True,exist_ok=True)
    data=dict(rows=rows,controls=controls,cpu=time.process_time()-t)
    (out/"results.json").write_text(json.dumps(data,indent=2)+"\n")
    print(len(rows),sum(x["finite"] is not None for x in rows),sum(x["finite"] is not None and x["finite"]["first"] is None for x in rows),controls,flush=True)
