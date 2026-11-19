"""EX-GC14 summary (reads records only): survival, knot degree, knot-core distance, far G vs point sink."""
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc09_summary as s9  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC14")


def main():
    gp = s9.green_profile(0.5, 6)
    summ = {}
    for f in sorted(glob.glob(os.path.join(OUT, "*_s*.json"))):
        d = json.load(open(f))
        r = d["rows"]
        late = [x for x in r if x["t"] >= 30]
        rec = dict(case=d["case"], lam=d["lam"], E=d["E"], seed=d["seed"], M0=r[0]["M"], M_final=r[-1]["M"],
                   alive=r[-1]["Wmax"] > 1, degrees=sorted({x["deg"] for x in r}),
                   M_late=float(np.mean([x["M"] for x in late])) if late else None,
                   knot_core_dist_late=float(np.mean([np.linalg.norm(np.subtract(x["Kcen"], x["Wcen"]))
                                                       for x in late])) if late and r[-1]["Wmax"] > 1 else None,
                   Wcen_z_final=r[-1]["Wcen"][2])
        if rec["alive"] and d["E"] == 0:
            rec["far_ratio"] = s9.far_ratio(d["Gprofile"], r[-1]["M"], gp, 3.5, 5)[0]
        summ[os.path.basename(f)[:-5]] = rec
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    for k, v in summ.items():
        print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a not in ("case",)})


if __name__ == "__main__":
    main()
