"""EX-GC36B Step 1 separate path (Claude). Uses CRR x6_generic.py helpers (rel, decode_mask,
RuleSpec.outputs, load_specs) and re-implements only the orientation-free gate search. Nothing imported from exgc36b.
Counts rules whose every phase has >=1 orientation-free move and all such moves land in a phase up to isomorphism."""
import sys, pathlib, itertools, collections, json
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "external" / "crr"))
import x6_generic as X6

def gates_free(R):
    by = collections.defaultdict(set)
    for r in R: by[r[0]].add(r)
    out = set()
    for x, xr in by.items():
        for r1 in xr:
            for y, z in ((r1[1], r1[2]), (r1[2], r1[1])):
                if X6.rel(y, z, x) not in R or len({x, y, z}) < 3: continue
                for r3 in xr - {r1}:
                    for u, v in ((r3[1], r3[2]), (r3[2], r3[1])):
                        if len({x, y, z, u, v}) == 5 and X6.rel(v, y, z) in R: out.add((x, y, z, u, v))
    return out

def iso_phase(R, S, phases):
    Sset = set(S)
    for p in itertools.permutations(S):
        mp = {v: i for i, v in enumerate(p)}
        rr = {(mp[a], min(mp[b], mp[c]), max(mp[b], mp[c])) for a, b, c in R if {a, b, c} <= Sset}
        if rr in phases: return True
    return False

specs = X6.load_specs(str(ROOT / "external" / "crr" / "W21_CLOSED_PARTICLE_HITS.json"))
passing = []
for idx, sp in specs.items():
    ok = True
    for R0 in sp.phase_rels:
        gs = gates_free(R0)
        if not gs: ok = False; break
        for g in gs:
            x, y, z, u, v = g; R = set(R0) - {X6.rel(x, y, z), X6.rel(y, z, x), X6.rel(x, u, v)}
            R |= set(sp.outputs(y, z, u, v, 5))
            if not iso_phase(R, (y, z, u, v, 5), sp.phase_rels): ok = False; break
        if not ok: break
    if ok: passing.append(idx)
print(json.dumps(dict(rules=len(specs), C_closed=passing)))
(ROOT / "results" / "EX-GC36B" / "recompute.json").write_text(json.dumps(dict(rules=len(specs), C_closed=passing)) + "\n")
