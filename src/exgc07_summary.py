"""EX-GC07 summary figure and numbers: consumption M(t) of the clocked particle in the 3D (spectral-angle) runs."""
import gzip
import json
import os

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC07")
RUNS = [("k1_tV4", "κ = 1, τ_V = 4", (0.0512, 0.3083)), ("k1_tV20", "κ = 1, τ_V = 20", (0.300, 0.0)),
        ("k0.9_tV7", "κ = 0.9, τ_V = 7", (0.0228, 0.2921)), ("k0.9_tV6", "κ = 0.9, τ_V = 6", (None, None)),
        ("k0.9_tV50", "κ = 0.9, τ_V = 50", (0.1018, 0.0547)), ("k0.9_tV5.5_big", "κ = 0.9, τ_V = 5.5, kick ×50", (None, None))]
summary = {}
fig, ax = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True)
for i, (tag, lab, lin) in enumerate(RUNS):
    d = json.load(gzip.open(os.path.join(OUT, f"ts_{tag}.json.gz"), "rt"))
    t = np.array([r["t"] for r in d["rows"]]); M = np.array([r["M"] for r in d["rows"]])
    l1 = np.array([r["l1"] for r in d["rows"]])
    dead = M < 0.5 * M[0]
    rec = dict(label=lab, M0=M[0], collapse_t=float(t[np.argmax(dead)]) if dead.any() else None, T=float(t[-1]),
               linear=lin, wall_s=d.get("wall"), l1_max=float(l1.max()), centroid_max=float(max(
                   np.hypot(np.hypot(r["cx"], r["cy"]), r["cz"]) for r in d["rows"])))
    a = np.abs(M - M[0]); pk = [j for j in range(1, len(a) - 1) if a[j] >= a[j - 1] and a[j] >= a[j + 1] and 1e-6 < a[j] < 0.05 * M[0]]
    if len(pk) > 4:
        rec["growth_fit"] = float(np.polyfit(t[pk], np.log(a[pk]), 1)[0])
        rec["omega_fit"] = float(np.pi / np.median(np.diff(t[pk])))
    summary[tag] = rec
    ax[0].plot(t, M / M[0], lw=0.8, label=lab)
    ax[1].semilogy(t, np.maximum(a / M[0], 1e-8), lw=0.6, label=lab)
ax[0].set_ylabel("consumption M(t) / M(0)")
ax[0].legend(fontsize=7, ncol=2)
ax[1].set_ylabel("|M − M(0)| / M(0)"); ax[1].set_xlabel("time")
fig.suptitle("EX-GC07 (EXPLORATORY): breathing grows until the clocked core collapses; no circling seen", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC07_breathing.png"), dpi=150)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
for k, v in summary.items():
    print(k, {a: (round(b, 4) if isinstance(b, float) else b) for a, b in v.items() if a != "label"})
