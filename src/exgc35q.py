"""Finite words and exact event-dependency comparisons for GC35Q."""
from itertools import product
import json
from pathlib import Path
import time
import networkx as nx
import exgc35l as prior
import exgc35p as pools

RULES={
    "A":("01",(("10","01"),)),
    "B":("0LR",(("R0","0R"),("0L","L0"),("RL","LR"))),
    "C":("0LR",(("R0","0R"),("0L","L0"))),
    "D":("0GP",(("GP","P"),("PG","P"),("G0","0G"))),
    "E":("0GP",(("GP0","00P"),("PG","P"),("G0","0G"))),
}
SWAP=str.maketrans("RL","LR")


def matches(word,rules):
    return [(i,k) for i in range(len(word)) for k,(a,_) in enumerate(rules) if word.startswith(a,i)]


def outcomes(word,rules,forced=None):
    def visit(w,roots,events,first=None):
        steps=[first] if first is not None else matches(w,rules)
        if not steps:
            yield w,events
            return
        for i,k in steps:
            a,b=rules[k]
            end=i+len(a)
            ev=(i,k,tuple(sorted({v for v in roots[i:end] if v>=0})))
            outroot=(len(events),)*len(b)
            yield from visit(w[:i]+b+w[end:],roots[:i]+outroot+roots[end:],events+(ev,))
    return visit(word,(-1,)*len(word),(),forced)


def graph(events,word,rules,mode,mirror=False):
    g=nx.MultiDiGraph()
    for i,(site,k,parents) in enumerate(events):
        if mirror:
            site=len(word)-site-len(rules[k][0])
            k=(1,0,2)[k]
        g.add_node(i,label=None if mode=="plain" else (site,k))
        for p in parents:
            g.add_edge(p,i,kind="producer")
    return g


def equivalent(a,b,start,rules,mode):
    if mode=="terminal":
        return a[0]==b[0]
    if a[1]==b[1]:
        return True
    ga=graph(a[1],start,rules,mode)
    gb=graph(b[1],start,rules,mode)
    if prior.equivalent(ga,gb):
        return True
    return mode=="mirror" and start[::-1].translate(SWAP)==start and prior.equivalent(
        ga,graph(b[1],start,rules,mode,mirror=True))


def compare(start,rules,kind,first=None,second=None):
    modes=["terminal","labelled","plain"]+(["mirror"] if kind=="B" else [])
    a=outcomes(start,rules,first)
    baseline=next(a)
    fails={}
    count=[1,0]
    def check(ref,row):
        for mode in modes:
            if mode not in fails and not equivalent(ref,row,start,rules,mode):
                fails[mode]=dict(word=start,a=ref,b=row,first=first,second=second)
    for x in a:
        count[0]+=1
        check(baseline,x)
    if second is not None:
        for x in outcomes(start,rules,second):
            count[1]+=1
            check(baseline,x)
    return fails,count


def overlaps(rules):
    result=set()
    for i,(a,_) in enumerate(rules):
        for j,(b,_) in enumerate(rules):
            for d in range(len(a)):
                if d==0 and i>=j:
                    continue
                end=min(len(a)-d,len(b))
                if a[d:d+end]!=b[:end]:
                    continue
                word=a+b[end:]
                result.add((word,(0,i),(d,j)))
    return sorted(result,key=lambda x:(len(x[0]),x))


def screen(kind,alphabet,rules):
    stats={}
    witnesses={}
    words=total=0
    for n in range(1,7):
        for chars in product(alphabet,repeat=n):
            word="".join(chars)
            fails,counts=compare(word,rules,kind)
            words+=1
            total+=sum(counts)
            for mode,case in fails.items():
                stats[mode]=stats.get(mode,0)+1
                witnesses.setdefault(mode,case)
    pairs=[]
    for word,a,b in overlaps(rules):
        fails,counts=compare(word,rules,kind,a,b)
        pairs.append(dict(word=word,first=a,second=b,completions=counts,failures=fails))
    return dict(words=words,terminal_histories=total,failure_words=stats,first=witnesses,overlaps=pairs)


if __name__=="__main__":
    start=time.process_time()
    rows={k:screen(k,*v) for k,v in RULES.items()}
    rows["F"]={str(x):dict(zip(("terminal_histories","failures"),prior.case(x))) for x in ((3,),(2,2),(2,3,2))}
    for x in ((3,),(2,2),(2,3,2)):
        failures,paths=pools.graph_forks("B",pools.initial(x,len(x),False),{"labelled","unlabelled"})
        rows["F"][str(x)].update(pool_graph_failures=failures,pool_graph_paths=paths)
    p=Path(__file__).resolve().parents[1]/"results"/"EX-GC35Q"
    p.mkdir(parents=True,exist_ok=True)
    (p/"results.json").write_text(json.dumps(dict(candidates=rows,cpu_seconds=time.process_time()-start),indent=2)+"\n")
    print({k:(v["words"],v["failure_words"],[r["word"] for r in v["overlaps"]]) for k,v in rows.items() if k!="F"},flush=True)
