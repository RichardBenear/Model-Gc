"""EX-GC36A separate recomputation path (Claude). Shared dependencies: networkx; CRR x6_generic.py (dynamics only, for
check 2); exgc35v.py (ChatGPT's 35V engine, for check 1). Nothing is imported from exgc36a.py.
1. Calibration: the published rule xy,zy->xw,yw,zw on the asymmetric 3-in star + leaf edge, through the 35V engine with
   vertex labels removed and no input anchors (as in review/EX-GC35V_control_check_claude.py), to depth 3.
2. Particle fork: from isolated phase B, enumerate all event sequences of length 1..4 per branch with the CRR simulator,
   and test whether ANY left and right final relation sets are isomorphic as plain hypergraphs (state confluence, a
   weaker condition than causal-graph equivalence)."""
import sys, pathlib, itertools, json
import networkx as nx
from networkx.algorithms.isomorphism import categorical_node_match, categorical_multiedge_match
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "external" / "crr"))
import exgc35v as V
import x6_generic as X6

out = {}
# 1
orig = V.graph
def unl(state, seed, types, anchor=True):
    g = orig(state, seed, types, anchor)
    for n, d in g.nodes(data=True):
        if n[0] == "v": d["label"] = "vtx"
    return g
V.graph = unl
rule = ((("C", 0, 1), ("C", 2, 1)), (("C", 0, 5), ("C", 1, 5), ("C", 2, 5)))
seed = (("C", 0, 3), ("C", 1, 3), ("C", 2, 3), ("C", 4, 0))
st = V.make(seed)
out["calibration_35v_engine"] = [V.join(rule, st, a, b, depth=3, anchor=False)
                                 for a, b in itertools.combinations(V.choices(rule, st[0]), 2)]
V.graph = orig

# 2
spec = X6.load_specs(str(ROOT / "external" / "crr" / "W21_CLOSED_PARTICLE_HITS.json"))[55744]
B = sorted(spec.phase_rels[1])
def sim_from(rels, n):
    s = X6.Sim(spec, 0)
    for r in list(s.rels): s._remove(r)
    s.gates_by_x = {}; s.gates_set = set()
    for r in rels: s._add(r)
    s.n = n; s._recompute_xs(set(v for r in rels for v in r)); return s
def hg(rels):
    g = nx.MultiDiGraph()
    for j, (a, b, c) in enumerate(sorted(rels)):
        g.add_node(("r", j), lab="r")
        for v in (a, b, c): g.add_node(("v", v), lab="v")
        g.add_edge(("r", j), ("v", a), kind="s"); g.add_edge(("r", j), ("v", b), kind="t"); g.add_edge(("r", j), ("v", c), kind="t")
    return g
def same(a, b):
    return len(a) == len(b) and nx.is_isomorphic(hg(a), hg(b), node_match=categorical_node_match("lab", None),
                                                 edge_match=categorical_multiedge_match("kind", None))
root = sim_from(B, 5)
g1, g2 = sorted(root.gates_set)
L = [root.clone()]; L[0].fire(g1)
R = [root.clone()]; R[0].fire(g2)
rows = []
for level in range(0, 5):
    hit = any(same(a.rels, b.rels) for a in L for b in R)
    rows.append(dict(level=level, left=len(L), right=len(R), any_isomorphic_state=hit))
    if hit or level == 4: break
    def ext(S):
        o = []
        for s in S:
            for g in sorted(s.gates_set):
                c = s.clone(); c.fire(g); o.append(c)
        return o
    L, R = ext(L), ext(R)
out["particle_state_confluence_L"] = rows
print(json.dumps(out, indent=1, default=list))
(ROOT / "results" / "EX-GC36A").mkdir(parents=True, exist_ok=True)
(ROOT / "results" / "EX-GC36A" / "recompute.json").write_text(json.dumps(out, indent=1, default=list) + "\n")
