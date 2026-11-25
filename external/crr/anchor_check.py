"""Pre-design check (Claude, 2026-11-25): run CRR x6_generic.py Rule 55744 to 5,000 events, compare with W37 row 0,
and count competing gates over the first 2,000 events. Run: python external/crr/anchor_check.py (about 3 minutes)."""
import sys, pathlib, collections, statistics
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from x6_generic import RuleSpec, Sim, rel
T = [[1,4,5],[2,1,3],[2,4,5],[3,1,4],[3,2,5],[3,4,5],[4,2,5],[5,2,3]]   # W21_CLOSED_PARTICLE_HITS template
s = Sim(RuleSpec(55744, T, []), 0)
widths = []; pairs = share = readwrite = 0
for i in range(5000):
    G = list(s.gates_set); widths.append(len(G))
    if i < 2000 and len(G) > 1:
        cons = {g: {rel(g[0],g[1],g[2]), rel(g[1],g[2],g[0]), rel(g[0],g[3],g[4])} for g in G}
        gate = {g: rel(g[4],g[1],g[2]) for g in G}
        for a in range(len(G)):
            for b in range(a+1, len(G)):
                pairs += 1
                if cons[G[a]] & cons[G[b]]: share += 1
                elif gate[G[a]] in cons[G[b]] or gate[G[b]] in cons[G[a]]: readwrite += 1
    if s.step() is None: break
inc = collections.Counter(v for r in s.rels for v in r)
print("nodes", s.n, "relations", len(s.rels), "max incidence", max(inc.values()),
      "max 2-section degree", max(len(v) for v in s.twosection().values()))
print("legal gates per state: mean", round(statistics.mean(widths), 1), "min", min(widths), "max", max(widths))
print("gate pairs (first 2000 events)", pairs, "share a consumed relation", share, "gate read vs consume", readwrite)
print("platform", sys.version.split()[0])
