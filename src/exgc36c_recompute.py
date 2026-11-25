"""EX-GC36C separate tracker (Claude). Replays schedules with the same simulator and seeds and, after each event,
scans every candidate 5-set for a phase support (36B canonical form). Candidates are all unions of the vertex sets of two
live relations that share exactly one vertex. Every phase (10-11 relations on 5 vertices) contains such a pair, so no support is
missed. (A first version scanned only the 12 most recent vertices; it missed supports on older vertices and was replaced.) Compares the first
event at which no support is found with the incremental tracker's lifetime (first 100 events only, for cost)."""
import sys, pathlib, json, itertools
ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import exgc36c as C
from exgc36b import canon5, decode

TMAX = 100
out = {}
for idx, version in ((66427, "C"), (55744, "L"), (66427, "L")):
    e = C.DATA[idx]; classes = {canon5(decode(m), range(5)) for m in e["state_masks"]}
    for ph in range(len(e["state_masks"])):
        rows = []
        for sd in range(5):
            s = C.start(idx, version, ph, sd); t = 0; dead = None
            while t < TMAX and s.gates_set:
                g = s.rng.choice(tuple(s.gates_set)); s.fire(g); t += 1
                R = list(s.rels); cand = set()
                for a, b in itertools.combinations(R, 2):
                    if len(set(a) & set(b)) == 1: cand.add(tuple(sorted(set(a) | set(b))))
                if not any(canon5(s.rels, S) in classes for S in cand):
                    dead = t; break
            ref = C.run(idx, version, ph, sd)
            rows.append(dict(seed=sd, scanner_first_empty=dead, tracker_lifetime=ref["lifetime"]))
        out[f"{idx}{version}_phase{ph}"] = rows
        print(f"{idx}{version}_phase{ph}", rows, flush=True)
(ROOT / "results" / "EX-GC36C").mkdir(parents=True, exist_ok=True)
(ROOT / "results" / "EX-GC36C" / "recompute.json").write_text(json.dumps(out, indent=1) + "\n")
