"""GC35R finite systematic rule search."""
import json
from itertools import product, combinations
from pathlib import Path
import time
import networkx as nx
import exgc35q as Q

LIMIT=25000


class Cycle(Exception):
    pass


class Budget(Exception):
    pass


def outcomes(start,rules,first=None):
    count=0
    def visit(word,roots,ev,active,forced=None):
        nonlocal count
        count+=1
        if count>LIMIT:
            raise Budget
        if word in active:
            raise Cycle
        m=Q.matches(word,rules) if forced is None else [forced]
        if not m:
            yield word,roots,ev
        for i,j in m:
            a,b=rules[j]
            stop=i+len(a)
            yield from visit(word[:i]+b+word[stop:],roots[:i]+(len(ev),)*len(b)+roots[stop:],
                             ev+(tuple(sorted(set(roots[i:stop]))),),active|{word})
    return visit(start,tuple(-i-1 for i in range(len(start))),(),set(),first)


def graph(row,mode):
    word,roots,ev=row
    g=nx.MultiDiGraph()
    for i,parents in enumerate(ev):
        g.add_node(i,label="event")
        for p in parents:
            if p>=0 or mode=="strict":
                if p<0:
                    g.add_node(p,label=("input",-p-1))
                g.add_edge(p,i,kind="producer")
    if mode!="plain":
        for i,p in enumerate(roots):
            out=len(ev)+i
            g.add_node(out,label=("output",i,word[i]))
            if p>=0 or mode=="strict":
                if p<0:
                    g.add_node(p,label=("input",-p-1))
                g.add_edge(p,out,kind="producer")
    return g


def score(word,rules,modes,first=None,second=None):
    failures={}
    refs={}
    baseline=None
    count=0
    for branch in ((first,second) if second is not None else (first,)):
        for row in outcomes(word,rules,branch):
            count+=1
            if baseline is None:
                baseline=row
                refs={m:graph(row,m) for m in modes}
            for m in modes:
                if m in failures:
                    continue
                if baseline[0]!=row[0] or not Q.prior.equivalent(refs[m],graph(row,m)):
                    failures[m]=(baseline,row)
    return failures,count


def catalogue(alphabet,limit):
    words=["".join(x) for n in range(1,limit+1) for x in product(alphabet,repeat=n)]
    rules=[(a,b) for a in words for b in words if len(b)<=len(a) and a!=b and a.count("1")==b.count("1")]
    counters={}
    examples={}
    survivors=[]
    def mark(k,s):
        counters[k]=counters.get(k,0)+1
        examples.setdefault(k,s)
    for size in (1,2):
        for rs in combinations(rules,size):
            overlaps=Q.overlaps(rs)
            if not any("1" in w for w,_,_ in overlaps):
                continue
            direct=any(all("1" in rs[x[1]][0] for x in (a,b)) for _,a,b in overlaps)
            mark("eligible",rs)
            if direct:
                mark("both_matches_consume_carrier",rs)
            bad={}
            try:
                for w,a,b in overlaps:
                    f,_=score(w,rs,("output","strict"),a,b)
                    for k,v in f.items():
                        bad.setdefault(k,dict(word=w,first=a,second=b,histories=v))
            except Cycle:
                mark("local_cycle",rs)
                continue
            except Budget:
                mark("local_cap",rs)
                continue
            if "output" in bad:
                mark("output_border_failure",dict(rules=rs,witness=bad["output"]))
                continue
            mark("output_border_survivor",rs)
            if "strict" not in bad:
                mark("strict_border_survivor",rs)
            if direct:
                mark("output_survivor_both_consume_carrier",rs)
            survivors.append((rs,direct,"strict" not in bad,bad))
    full=[]
    for rs,direct,strict,local in survivors:
        outcome="pass";witness=None;screened=hist=0
        try:
            for n in range(1,8):
                for chars in product(alphabet,repeat=n):
                    w="".join(chars)
                    fail,count=score(w,rs,("plain",))
                    screened+=1;hist+=count
                    if fail:
                        outcome="fail"
                        witness=dict(word=w,histories=fail["plain"])
                        break
                if witness:
                    break
        except Cycle:
            outcome="cycle"
        except Budget:
            outcome="cap"
        mark("full_"+outcome,rs)
        full.append(dict(rules=rs,both_consume_carrier=direct,strict_border=strict,local=local,
                         result=outcome,witness=witness,screened=screened,histories=hist))
    return dict(alphabet=alphabet,max_side=limit,rule_count=len(rules),counts=counters,examples=examples,full=full)


if __name__=="__main__":
    start=time.process_time()
    rows=[catalogue("01",3),catalogue("012",2)]
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35R"
    p.mkdir(parents=True,exist_ok=True)
    (p/"results.json").write_text(json.dumps(dict(rows=rows,cpu_seconds=time.process_time()-start),indent=2)+"\n")
    print([(r["alphabet"],r["rule_count"],r["counts"]) for r in rows],flush=True)
