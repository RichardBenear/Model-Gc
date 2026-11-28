"""EX-GC63 recomputation -> results/EX-GC63/recompute.json. R_dyn by an independent least-squares design matrix (and a t >= 2
fit); 63B moments by explicit link loops over x x j. Shared: run records, packet builders, operators. EXPLORATORY. Our
construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC63"; L = lambda k: json.loads((OUT / f"EX-GC63_{k}.json").read_text())
def F1(V, tmin=0.0):
    p = L(f"A_{V}_pair")["runs"]["pair"]; o = L(f"A_{V}_one")["runs"]["one"]; t = np.array(p["t"])
    y = np.array([a[0] - b[0] for a, b in zip(p["Py"], o["Py"])]); k = t >= tmin - 1e-9
    A = np.stack([np.ones(k.sum()), t[k], t[k] ** 2], 1); return float(np.linalg.lstsq(A, y[k], rcond=None)[0][1])
res = dict(R=F1("0.2") / F1("0.0") / np.sqrt(0.96), R_tge2=F1("0.2", 2) / F1("0.0", 2) / np.sqrt(0.96))
import exgc63 as Q, exgc59 as E59, exgc53 as T53, exgc58 as X
n = 97; c = n // 2
for nm, chi in (("up", E59.UP), ("down", E59.DOWN)):
    p = E59.source(n, chi); s = Q.d2e1(p, E59.M); PS = T53.link_P(p) + 0.25 * X.curl_T(X.sdens(p)); mz = {}
    for lab, j in (("D2E1", s), ("PS", PS)):
        x = np.arange(n) - c; X3, Y3, _ = np.meshgrid(x, x, x, indexing="ij")
        mz[lab] = float(0.5 * (np.sum(X3 * j[1]) - np.sum(Y3 * j[0])))   # midpoint offsets are along the link's own axis, so they drop out of the z-moment
    res[nm] = dict(mz=mz, ratio_z=abs(mz["D2E1"] / mz["PS"]))
(OUT / "recompute.json").write_text(json.dumps(res, indent=1)); print(json.dumps(res, indent=1))
