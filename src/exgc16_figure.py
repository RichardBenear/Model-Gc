"""EX-GC16 figure: fall rate by particle type (weak, strong consumption) and the two-unit pull versus distance."""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc16_recoil as r  # noqa: E402
import exgc16_summary as s  # noqa: E402

rows = json.load(open(os.path.join(r.OUT, "summary.json")))
kinds = ["single", "double_site", "spread7", "pair_z", "pair_x"]
labels = ["single", "2x rate\none site", "spread\n7 sites", "two units\nstacked", "two units\nside by side"]
fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.2))
for j, (m0, col) in enumerate(((50, "#1f77b4"), (500, "#d62728"))):
    ys, es, ps = [], [], []
    for k in kinds:
        rw = next(x for x in rows if x["kind"] == k and x["m0"] == m0 and x["E"] == 0.02)
        ys.append(rw["fall_z"] / rw["ideal_fall_z"])
        es.append(rw["fall_z_sem"] / abs(rw["ideal_fall_z"]))
        ps.append(rw["pred_fall"][2] / rw["ideal_fall_z"])
    xx = np.arange(len(kinds)) + (j - 0.5) * 0.25
    a1.errorbar(xx, ys, es, fmt="o", color=col, label=f"measured, m0 = {m0}")
    a1.plot(xx, ps, "_", ms=14, mew=2, color=col, alpha=0.6, label=f"mean-field prediction, m0 = {m0}")
a1.axhline(1, color="k", lw=0.8, ls="--")
a1.set_xticks(range(len(kinds)), labels, fontsize=8)
a1.set_ylabel("fall rate / ideal (m0 E / 3)")
a1.set_title("Free fall in a G gradient, by particle type")
a1.legend(fontsize=7)

pr = sorted([x for x in rows if x["kind"] == "pair_x" and x["E"] == 0 and "pair_pull" in x], key=lambda x: x["d"])
d = np.array([x["d"] for x in pr])
pull = np.array([x["pair_pull"] for x in pr])
err = []
for x in pr:
    files = [os.path.join(r.OUT, f"pair_x_m500_E0_d{x['d']}_s{sd}.json") for sd in range(x["seeds"])]
    v = [(lambda q: (q["force_exp"][0][0] - q["force_exp"][1][0]) / 2)(json.load(open(f))) for f in files]
    err.append(np.std(v, ddof=1) / np.sqrt(len(v)))
dd = np.linspace(2, 9, 50)
pd = []
L = pr[0]["L"]
for dv in range(2, 10):
    u, _ = r.make_units("pair_x", L, 500, dv)
    F = s.pred_force(u, s.mean_field(u, L, 1000, 0))
    pd.append((F[0, 0] - F[1, 0]) / 2)
a2.errorbar(d, pull, err, fmt="o", color="#d62728", label="measured (token lattice)")
a2.plot(range(2, 10), pd, "-", color="#555", label="mean-field lattice prediction")
a2.plot(dd, 500 ** 2 / (12 * np.pi * 1000) / dd ** 2, ":", color="k", label=r"continuum $m_0^2/(12\pi\Omega d^2)$")
a2.set_xscale("log")
a2.set_yscale("log")
a2.set_xlabel("separation d (lattice units)")
a2.set_ylabel("pull of each unit toward the other")
a2.set_title("Two consumers, no gradient (m0 = 500)")
a2.legend(fontsize=7)
fig.suptitle("EX-GC16 (EXPLORATORY): recoil rule on a G token lattice — our construction", fontsize=10)
fig.tight_layout()
out = os.path.join(r.ROOT, "paper", "figures", "exgc16_recoil.png")
fig.savefig(out, dpi=130)
print(out)
