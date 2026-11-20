"""EX-GC18 figure: core (W) and knot positions, and M, for the recoil runs."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "EX-GC18")
runs = [("plain_K0.2_E0.02_rec1_s0", "rigid shift, K 0.2", "#1f77b4", "-"),
        ("plain_K0.4_E0.02_rec1_s0", "rigid shift, K 0.4", "#1f77b4", "--"),
        ("knot_K0.2_E0.02_rec1_s0", "knot kicked, K 0.2", "#d62728", "-"),
        ("knot_K0.4_E0.02_rec1_s0", "knot kicked, K 0.4", "#d62728", "--"),
        ("knot_K0.2_E0.02_rec1_s0_field", "W pushed (field), K 0.2", "#2ca02c", "-"),
        ("knot_K0.4_E0.02_rec1_s0_field", "W pushed (field), K 0.4", "#2ca02c", "--"),
        ("knot_K0.2_E0_rec1_s0", "knot, no gradient", "#7f7f7f", ":")]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for tag, lab, c, ls in runs:
    r = json.load(open(os.path.join(D, tag + ".json")))
    t = [x["t"] for x in r["rows"]]
    ax[0].plot(t, [x["Wcen"][2] for x in r["rows"]], ls, color=c, label=lab + " (core)")
    if r["rows"][0]["Kcen"] is not None and "field" not in tag and "E0_" not in tag:
        ax[0].plot(t, [x["Kcen"][2] for x in r["rows"]], ls, color=c, alpha=.4, lw=3, label=lab + " (knot)")
    ax[1].plot(t, [x["M"] / r["rows"][0]["M"] for x in r["rows"]], ls, color=c, label=lab)
ax[0].plot([0, 29], [0, -0.048 * 29], "k:", lw=.8, label="recoil speed 0.048")
ax[0].set_xlabel("t")
ax[0].set_ylabel("z position")
ax[0].set_title("Where the particle goes (gradient E = 0.02, down = lower G)")
ax[0].legend(fontsize=6)
ax[1].set_xlabel("t")
ax[1].set_ylabel("M / M(0)")
ax[1].set_title("Does the particle survive?")
ax[1].legend(fontsize=7)
fig.suptitle("EX-GC18 (EXPLORATORY): recoil rule on the knot particle — our construction", fontsize=10)
fig.tight_layout()
out = os.path.join(ROOT, "paper", "figures", "exgc18_knotrecoil.png")
fig.savefig(out, dpi=130)
print(out)
