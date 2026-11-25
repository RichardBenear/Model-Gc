"""Post hoc diagnostic of EX-GC35L (Claude review, not a registered readout).

Re-scores the GC35L enumeration with two coarser graph comparisons:
  mirror     - labelled graphs, but also accept a match after reflecting site/sink labels
               (only a symmetry when the initial configuration is a palindrome)
  unlabelled - event/initial/terminal node kinds and edge kinds only; site labels dropped
               (closer to Wolfram-style causal-graph comparison, where events carry no place)
Same domain and rule as src/exgc35l.py, which it imports unchanged.
"""
import json, time
from itertools import product
from pathlib import Path
import networkx as nx
from networkx.algorithms.isomorphism import categorical_node_match, categorical_multiedge_match
import exgc35l as L


def relabel(g, f):
    h = g.copy()
    for v, d in h.nodes(data=True):
        d["label"] = f(d["label"])
    return h


def kind(lbl):
    return lbl.split()[0]


def mirror_label(n):
    def f(lbl):
        k, i = lbl.split()
        i = int(i)
        return f"{k} {n + 1 - i}"     # sites 0..n+1 incl. sinks; initial tokens carry their site
    return f


def iso(a, b):
    return nx.is_isomorphic(a, b, node_match=categorical_node_match("label", ""),
                            edge_match=categorical_multiedge_match("kind", ""))


def case(start):
    n = len(start)
    pal = tuple(start) == tuple(reversed(start))
    initial = ((),) + tuple(tuple([-i - 1] * v) for i, v in enumerate(start, 1)) + ((),)
    ref = {}
    fail = {}
    seen = set()

    def score(name, g):
        if name not in ref:
            ref[name] = (g, relabel(g, kind))
            return
        if fail.get(name + "/unlabelled"):
            return
        a, au = ref[name]
        if not iso(a, g):
            fail[name + "/labelled"] = True
            if not (pal and iso(relabel(a, mirror_label(n)), g)):
                fail[name + "/mirror"] = True
            if not iso(au, relabel(g, kind)):
                fail[name + "/unlabelled"] = True

    def walk(tok, fifo, te, fe):
        eligible = [i for i in range(1, n + 1) if len(tok[i]) >= 2]
        if eligible:
            for i in eligible:
                fp = fifo[i][:2]
                fn = L.update(fifo, i, fp, len(te))
                for pair in L.choices(tok[i]):
                    walk(L.update(tok, i, pair, len(te)), fn, te + ((i, pair),), fe + ((i, fp),))
            return
        score("token", L.graph(n, te, tok, False))
        key = fe
        if key not in seen:
            seen.add(key)
            score("site", L.graph(n, fe, fifo, True))

    walk(initial, initial, (), ())
    return fail


if __name__ == "__main__":
    t0 = time.process_time()
    counts, examples, total = {}, {}, 0
    for weight in range(8):
        for n in range(1, 5):
            for config in product(range(weight + 1), repeat=n):
                if sum(config) != weight:
                    continue
                total += 1
                for k in case(config):
                    counts[k] = counts.get(k, 0) + 1
                    examples.setdefault(k, list(config))
    res = dict(configurations=total, mismatch_configurations=counts, first_example=examples,
               cpu_seconds=time.process_time() - t0)
    out = Path(__file__).resolve().parents[1] / "review" / "EX-GC35L_posthoc_claude.json"
    out.write_text(json.dumps(res, indent=2) + "\n")
    print(res)
