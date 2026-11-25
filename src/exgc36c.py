"""EX-GC36C (Claude, 2026-11-25, EXPLORATORY): persistence of a CRR-rule particle in its own debris.
Run: python src/exgc36c.py  (writes results/EX-GC36C/results.json)"""
import json, sys, time, platform, pathlib, random, statistics
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "external" / "crr"))
import x6_generic as X6
from exgc36b import canon5, decode

T = 600; N = 200
DATA = {e["template_index"]: e for e in json.load(open(ROOT / "external" / "crr" / "W21_CLOSED_PARTICLE_HITS.json"))}


class SimC(X6.Sim):
    """x6 Sim with both orientations of {u,v} (version C). Everything else unchanged."""
    def _gates_for_x(self, x):
        out = set(); R = self.rels; xrels = self.by_src.get(x, ())
        if len(xrels) < 2: return out
        for r1 in xrels:
            _, p, q = r1
            for y, z in ((p, q), (q, p)):
                if len({x, y, z}) < 3 or X6.rel(y, z, x) not in R: continue
                for r3 in xrels:
                    if r3 == r1: continue
                    for u, v in ((r3[1], r3[2]), (r3[2], r3[1])):
                        if len({x, y, z, u, v}) < 5: continue
                        if X6.rel(v, y, z) in R: out.add((x, y, z, u, v))
        return out
    def gate_valid(self, g):
        x, y, z, u, v = g
        return (len(set(g)) == 5 and X6.rel(x, y, z) in self.rels and X6.rel(y, z, x) in self.rels
                and X6.rel(x, u, v) in self.rels and X6.rel(v, y, z) in self.rels)


def start(idx, version, phase, seed):
    e = DATA[idx]; spec = X6.RuleSpec(idx, e["template"], e["state_masks"])
    s = (SimC if version == "C" else X6.Sim)(spec, seed)
    for r in list(s.rels): s._remove(r)
    s.gates_by_x = {}; s.gates_set = set()
    for r in decode(e["state_masks"][phase]): s._add(r)
    s.n = 5; s._recompute_xs(set(range(5)))
    return s


def run(idx, version, phase, seed):
    e = DATA[idx]; classes = {canon5(decode(m), range(5)) for m in e["state_masks"]}
    s = start(idx, version, phase, seed)
    live = {tuple(range(5))}; renewals = 0; inside = 0; maxlive = 1; t = 0; halted = False
    while t < T:
        if not s.gates_set:
            halted = True; break
        g = s.rng.choice(tuple(s.gates_set))
        if any(set(g) <= set(S) for S in live): inside += 1
        w = s.fire(g); t += 1
        S2 = tuple(sorted((g[1], g[2], g[3], g[4], w)))
        live = {S for S in live if canon5(s.rels, S) in classes}
        if canon5(s.rels, S2) in classes:
            live.add(S2); renewals += 1
        maxlive = max(maxlive, len(live))
        if not live: break
    return dict(seed=seed, lifetime=t if not live else None, events=t, renewals=renewals, inside_fraction=inside / max(t, 1),
                max_live=maxlive, final_moves=len(s.gates_set), halted=halted, relations=len(s.rels))


def summarize(rows):
    lt = [r["lifetime"] if r["lifetime"] is not None else T + 1 for r in rows]
    alive500 = sum(1 for x in lt if x > 500)
    return dict(n=len(rows), alive_at_500=alive500, alive_at_T=sum(1 for x in lt if x > T),
                median_lifetime=statistics.median(lt), min_lifetime=min(lt), max_lifetime=max(lt),
                halted=sum(r["halted"] for r in rows), max_live=max(r["max_live"] for r in rows),
                mean_inside=statistics.mean(r["inside_fraction"] for r in rows))


def main():
    out = ROOT / "results" / "EX-GC36C"; out.mkdir(parents=True, exist_ok=True)
    res = dict(platform=dict(python=platform.python_version(), machine=platform.machine(), where="cloud container"), T=T, N=N, runs={})
    t0 = time.process_time()
    for idx, version in ((66427, "C"), (55744, "L"), (66427, "L")):
        for ph in range(len(DATA[idx]["state_masks"])):
            key = f"{idx}{version}_phase{ph}"
            rows = [run(idx, version, ph, sd) for sd in range(N)]
            res["runs"][key] = dict(summary=summarize(rows), rows=rows)
            print(key, res["runs"][key]["summary"], round(time.process_time() - t0, 1), flush=True)
            res["cpu_seconds"] = time.process_time() - t0
            (out / "results.json").write_text(json.dumps(res, indent=1) + "\n")


if __name__ == "__main__":
    main()
