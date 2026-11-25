import json
import time
from itertools import combinations, product
from pathlib import Path
import networkx as nx
from networkx.algorithms.isomorphism import categorical_node_match, categorical_multiedge_match


def choices(bag):
    return sorted(set(combinations(sorted(bag),2)))


def graph(n,events,bags,site_history):
    g=nx.MultiDiGraph()
    prev={}
    for i,(site,parents) in enumerate(events):
        g.add_node(i,label=f"event {site}")
        for p in parents:
            if p<0:
                g.add_node(p,label=f"initial {-p-1}")
            g.add_edge(p,i,kind="carrier")
        if site_history and site in prev:
            g.add_edge(prev[site],i,kind="site")
        prev[site]=i
    tail=len(events)
    for site,bag in enumerate(bags):
        if not bag:
            continue
        g.add_node(tail,label=f"terminal {site}")
        for p in bag:
            if p<0:
                g.add_node(p,label=f"initial {-p-1}")
            g.add_edge(p,tail,kind="carrier")
        tail+=1
    return g


def equivalent(a,b):
    if tuple(a.nodes(data=True))==tuple(b.nodes(data=True)) and tuple(sorted(a.edges(data="kind")))==tuple(sorted(b.edges(data="kind"))):
        return True
    return nx.is_isomorphic(a,b,node_match=categorical_node_match("label",""),edge_match=categorical_multiedge_match("kind",""))


def update(bags,i,pair,e):
    q=[list(b) for b in bags]
    for p in pair:
        q[i].remove(p)
    q[i-1].append(e)
    q[i+1].append(e)
    return tuple(tuple(b) for b in q)


def case(start):
    n=len(start)
    initial=((),)+tuple(tuple([-i-1]*v) for i,v in enumerate(start,1))+((),)
    known=None
    witness={}
    paths=0
    def walk(tok,fifo,te,fe):
        nonlocal known,paths
        eligible=[i for i in range(1,n+1) if len(tok[i])>=2]
        if eligible:
            for i in eligible:
                fp=fifo[i][:2]
                fn=update(fifo,i,fp,len(te))
                for pair in choices(tok[i]):
                    walk(update(tok,i,pair,len(te)),fn,te+((i,pair),),fe+((i,fp),))
            return
        paths+=1
        ledger=tuple(map(len,tok))
        counts=tuple(sum(site==i for site,_ in te) for i in range(1,n+1))
        payload=dict(events=te,site_events=fe,terminal_inputs=tok,site_terminal_inputs=fifo)
        if known is None:
            known=(ledger,counts,graph(n,te,tok,False),graph(n,fe,fifo,True),payload)
            return
        for index,name in enumerate(("ledger","firings","token","site")):
            if name in witness:
                continue
            if index==0: matches=ledger==known[0]
            elif index==1: matches=counts==known[1]
            else: matches=equivalent(known[index],graph(n,te if index==2 else fe,tok if index==2 else fifo,index==3))
            if not matches:
                witness[name]=dict(initial=start,a=known[-1],b=payload,ledger=ledger,firings=counts)
    walk(initial,initial,(),())
    return paths,witness


def shared_fork():
    # Follow fixed site firings until both candidates share an input.
    start=(2,3,2)
    b=((),)+tuple(tuple([-i-1]*v) for i,v in enumerate(start,1))+((),)
    prefix=[]
    for i in (1,3):
        p=b[i][:2]
        b=update(b,i,p,len(prefix));prefix.append((i,p))
    common=b[2][0]
    pair1=(common,0)
    pair2=(common,1)
    # The outer producers in the symmetric arrangement have distinct site labels.
    def complete(pair):
        x=update(b,2,pair,len(prefix));ev=prefix+[(2,pair)]
        while any(len(x[i])>=2 for i in range(1,4)):
            i=next(i for i in range(1,4) if len(x[i])>=2)
            p=x[i][:2];x=update(x,i,p,len(ev));ev.append((i,p))
        return ev,x
    a,xa=complete(pair1)
    b,xb=complete(pair2)
    return dict(initial=start,shared=common,pairs=[pair1,pair2],a=a,b=b,
                ledger_equal=tuple(map(len,xa))==tuple(map(len,xb)),
                token_equal=equivalent(graph(3,a,xa,False),graph(3,b,xb,False)))


if __name__=="__main__":
    start=time.process_time()
    counters={}
    first={}
    configurations=paths=0
    for weight in range(8):
        for n in range(1,5):
            for config in product(range(weight+1),repeat=n):
                if sum(config)!=weight:
                    continue
                p,w=case(config)
                configurations+=1;paths+=p
                for key,payload in w.items():
                    counters[key]=counters.get(key,0)+1
                    first.setdefault(key,payload)
    result=dict(configurations=configurations,paths=paths,mismatch_configurations=counters,
                first=first,shared_fork=shared_fork(),branch_control=[[1,0,0],[0,0,1]],
                cpu_seconds=time.process_time()-start)
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35L"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print({k:v for k,v in result.items() if k!="first"},flush=True)
