"""EX-GC42 separate recomputation of the variant-A coefficient from the stored axis profiles in far.json.
Independent of fit_A and the periodic Green's function: fits n(r) = -a/r + b + c r^2 (continuum 1/r, plus the periodic
box's constant and jellium terms) on r = 6..30 along +x. Shared dependencies: the stored n_axis (far_field) and dE0."""
import json
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1] / "results" / "EX-GC42"
far = json.loads((R / "far.json").read_text())["result"]
src = json.loads((R / "source.json").read_text())["result"]["cases"]
out = {}
for lam0, d in far.items():
    for key, row in d["cases"].items():
        hub, dl, fg = key.split("_")
        n = np.array(row["A"]["n_axis"]); r = np.arange(6, 31)
        X = np.column_stack([-1 / r, np.ones_like(r, float), r ** 2.0])
        a = np.linalg.lstsq(X, n[r], rcond=None)[0][0]
        ratio = a / (float(fg) / (4 * np.pi) * src[f"{lam0}_{hub}_{dl}"]["dE0"])
        out[f"{lam0}_{key}"] = dict(ratio_recomputed=float(ratio), ratio_fitA=row["A"]["ratio"],
                                    rel_diff=float(ratio / row["A"]["ratio"] - 1))
(R / "recompute.json").write_text(json.dumps(out, indent=1) + "\n")
print("max |rel diff|", max(abs(v["rel_diff"]) for v in out.values()))
for fg in ("0.01", "0.05"):
    for lam0 in ("0.1", "0.2"):
        v = [o["ratio_recomputed"] for k, o in out.items() if k.startswith(lam0) and k.endswith(fg)]
        print(lam0, fg, min(v), max(v), (max(v) - min(v)) / np.mean(v))
