"""EX-GC59: recompute the light-loop predictor after fixing its singular self-term (the stored value was -inf: a loop link coincided with a
source link carrying the source's exponential tail). The wave delays themselves are the stored scored values; nothing is re-run.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc59 as E, exgc53 as T53
f = E.OUT / "EX-GC59_main.json"; d = json.loads(f.read_text()); n = d["n"]; row = d["runs"]["axis_16"]
src, kap, N, Nd = E.source_terms(n); st = T53.square_loop(E.pos(n, "axis_16"), 8); A = 16.0 ** 2
row["light"]["Omega_pred_loop"] = -E.beta_gr_loop(src["up"]["PS"], kap, st, n) / (2 * A); row["light"]["ratio"] = row["light"]["Omega_light"] / row["light"]["Omega_pred_loop"]
row["light"]["note"] = "predictor recomputed with the singular self-term dropped (exgc59_fixpred.py)"; f.write_text(json.dumps(d, indent=1, default=float)); print(row["light"])
