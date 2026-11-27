"""EX-GC49 separate path and post hoc diagnostic.
R1: beta of each transported variable from the exact exterior form N = F^{-1}(-U), by a second finite difference of N^2
in U (no fit, no lattice). R2 (post hoc, not scored): the lattice fits share an offset that is linear in the source
strength Q (the same for all four f, so a linear-order lattice/boundary artefact). Separate it by fitting the x^2
coefficient a2(Q) = c1 Q + c2 Q^2 + c3 Q^3 across the four strengths; beta = 2 c2 / (a1/Q)^2.
EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json
from pathlib import Path
import numpy as np
import exgc49 as G

OUT = G.OUT


def r1():
    out = {}
    for v, V in G.VARS.items():
        h = 1e-4; n2 = lambda U: V["Fi"](-U) ** 2
        b = (n2(h) + n2(-h) - 2 * n2(0.0)) / (2 * h * h) / 2      # N^2 = 1 - 2U + 2 beta U^2
        out[v] = dict(beta_series=float(b), pred=V["beta"])
    return out


def r2():
    A = json.loads((OUT / "A.json").read_text())["result"]; out = {}
    for v, o in A.items():
        Qs = np.array([r["Q"] for r in o["rows"].values()])
        row = {}
        for key in ("fit_lattice", "fit_axis", "fit_diag", "fit_cont"):
            a1 = np.array([r[key]["coef"][0] for r in o["rows"].values()]); a2 = np.array([r[key]["coef"][1] for r in o["rows"].values()])
            M = np.stack([Qs, Qs ** 2, Qs ** 3], 1); c = np.linalg.lstsq(M, a2, rcond=None)[0]; k1 = np.mean(a1 / Qs)
            row[key] = dict(c_lin=float(c[0]), beta_sep=float(2 * c[1] / k1 ** 2), a1_over_Q=float(k1))
        out[v] = row
    return out


if __name__ == "__main__":
    res = dict(R1=r1(), R2_posthoc=r2())
    for k, v in res["R1"].items():
        print("R1", k, v)
    for k, v in res["R2_posthoc"].items():
        print("R2", k, {kk: round(vv["beta_sep"], 4) for kk, vv in v.items()}, "lin", {kk: round(vv["c_lin"], 5) for kk, vv in v.items()})
    (OUT / "recompute.json").write_text(json.dumps(dict(study="EX-GC49", part="recompute", label="EXPLORATORY; R2 post hoc, not scored", result=res), indent=1))
