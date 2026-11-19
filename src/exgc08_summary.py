"""EX-GC08 figure: tangle energy during annealing of twisted shells under the two no-cut rules."""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC08")
runs = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "B2_*.json")))]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
col = {(1, "lift"): "#c0504d", (2, "lift"): "#3b6ea5", (1, "flag"): "#9bbb59", (2, "flag"): "#8064a2",
       (0, "lift"): "0.5"}
summary = []
for r in runs:
    k = (r["turns"], r["rule"])
    panel = 0 if r["rule"] == "lift" else 1
    sw = [h["sweep"] for h in r["hist"]]
    E = [h["E"] for h in r["hist"]]
    lab = f"{r['turns']} turn{'s' if r['turns'] != 1 else ''}, S = {r['S']}, seed {r['seed']}, {r['sweeps']} sweeps"
    ax[panel].plot(sw, E, color=col.get(k, "k"), lw=1.2 if r["S"] == "90+120" else 0.8,
                   ls="-" if r["S"] == "90+120" else ":", label=lab)
    summary.append(dict(S=r["S"], rule=r["rule"], turns=r["turns"], seed=r["seed"], sweeps=r["sweeps"],
                        E_init=r["E_init"], E_final=r["E_final"], z2_init=r["z2_init"], z2_final=r["z2_final"],
                        z2_constant=len({h["z2"] for h in r["hist"]}) == 1))
ax[0].set_title("Rule with the three-node check (no half-turn ambiguity)", fontsize=9)
ax[1].set_title("Pairwise rule only (controls)", fontsize=9)
for a in ax:
    a.set_xlabel("annealing sweep"); a.legend(fontsize=6.5)
ax[0].set_ylabel("tangle energy (mismatched edges)")
fig.suptitle("EX-GC08 (EXPLORATORY): one turn stays tangled, two turns untangle (discrete belt trick)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC08_belt_trick.png"), dpi=150)
json.dump(summary, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
for s in summary:
    print(s)
