"""EX-GC25 figure from results/EX-GC25/EX-GC25_runs.json."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rec = json.load(open(os.path.join(ROOT, "results", "EX-GC25", "EX-GC25_runs.json")))
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
for mode, mk, lab in (("energy", "o", "coupled to energy"), ("charge", "s", "coupled to a non-energy charge")):
    r = [x for x in rec["fall"]["rows"] if x["mode"] == mode]
    ax[0].plot([x["M"] for x in r], [x["a_over_g"] for x in r], mk, ms=8, label=lab + " (measured)")
    ax[0].plot([x["M"] for x in r], [-x["a_pred"] / rec["fall"]["g"] for x in r], "k+", ms=12)
    for x in r:
        ax[0].annotate(x["kink"], (x["M"], x["a_over_g"]), textcoords="offset points", xytext=(6, 4))
ax[0].plot([], [], "k+", label="prediction")
ax[0].set_xscale("log"); ax[0].set_xlabel("kink rest energy M"); ax[0].set_ylabel("fall acceleration / g")
ax[0].set_title("Fall in a uniform potential"); ax[0].legend(fontsize=8)
for k, mk in (("A", "o-"), ("B", "s-")):
    r = [x for x in rec["coast"]["rows"] if x["kink"] == k]
    ax[1].plot([x["w_over_h"] for x in r], [max(abs(x["loss"]), 1e-6) for x in r], mk, label=f"kink {k}")
ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("kink width / lattice spacing")
ax[1].set_ylabel("|speed lost| over t = 200 (u = 0.5)"); ax[1].set_title("Coasting with no rest frame"); ax[1].legend()
fig.suptitle("EX-GC25 (EXPLORATORY): relativistic phi^4 kinks in 1+1 D, our construction", fontsize=10)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "exgc25_framefree.png"), dpi=150)
