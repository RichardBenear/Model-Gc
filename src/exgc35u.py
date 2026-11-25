"""GC35U: producible initial ports and a transport countercheck."""
import json
from itertools import product
from pathlib import Path
import time
import exgc35r as R

original_score=R.score
original_graph=R.graph


def score(start,rules,*args):
    letters=set("".join(b for a,b in rules))
    def graph(row,mode):
        g=original_graph(row,mode)
        if mode=="strict":
            g.remove_nodes_from([x for x in g if x<0 and start[-x-1] not in letters])
        return g
    R.graph=graph
    try:
        return original_score(start,rules,*args)
    finally:
        R.graph=original_graph


def transport(hop):
    rules=(("01","1"),("12","1"),hop)
    for n in range(1,6):
        for w in product(".012",repeat=n):
            start="".join(w)
            histories=list(R.outcomes(start,rules))
            g=R.graph(histories[0],"plain")
            for row in histories[1:]:
                if histories[0][0]==row[0] and not R.Q.prior.equivalent(g,R.graph(row,"plain")):
                    return dict(rule=hop,witness=start,a=histories[0],b=row)
    return dict(rule=hop,witness=None)


if __name__=="__main__":
    t=time.process_time()
    R.score=score
    rows=[R.catalogue("01",3),R.catalogue("012",2)]
    data=dict(catalogue=rows,transport=[transport(("0.",".0")),transport((".2","2."))],
              cpu=time.process_time()-t)
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35U"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(data,indent=2)+"\n")
    print([(r["alphabet"],r["counts"]) for r in rows], data["transport"],flush=True)
