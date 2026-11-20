"""Figures for the Model Gc white paper draft v0.2, computed from committed records only.

wp_loop_classes.png: EX-GC13 loop classes (results/EX-GC13/loop_*_s*.json).
wp_coasting_times.png: characteristic speed-loss times from EX-GC23 (note table), EX-GC24A (diagnostics.json) and
EX-GC25 (EX-GC25_runs.json). Time = interval / fractional speed loss (or -interval / ln retained).
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "paper", "figures")

# --- loop classes
loops = [("rot1", "one knot turned 360°", 1), ("rot2turns", "one knot turned 720°", 0),
         ("translate", "one knot carried round a circle", 0), ("rotB2", "degree-2 knot turned 360°", 0),
         ("exchange", "two knots swapped", 1), ("exchange_rotA", "swap + one knot turned 360°", 0),
         ("exchange_twice", "two knots swapped twice", 0), ("rotpair", "pair turned rigidly 360°", 0)]
M = np.zeros((len(loops), 3))
for i, (k, _, _) in enumerate(loops):
    for s in range(3):
        M[i, s] = json.load(open(os.path.join(ROOT, "results", "EX-GC13", f"loop_{k}_s{s}.json")))["loop_class"]
fig, ax = plt.subplots(figsize=(7.5, 4.2))
ax.imshow(M, cmap="RdBu_r", vmin=-0.6, vmax=1.6, aspect="auto")
for i, (_, lab, exp) in enumerate(loops):
    for s in range(3):
        ax.text(s, i, f"{int(M[i, s])}", ha="center", va="center", fontsize=12, color="k")
    ax.text(3.0, i, f"expected {exp}", va="center", fontsize=9)
ax.set_yticks(range(len(loops)), [l[1] for l in loops])
ax.set_xticks(range(3), ["value 1", "value 2", "value 3"])
ax.set_xlim(-0.5, 3.9)
ax.axhline(3.5, color="k", lw=1)
ax.set_title("EX-GC13 (EXPLORATORY): loop class of each motion of lattice knots\n"
             "1 = non-contractible loop (topological class), 0 = contractible; three independent regular values", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "wp_loop_classes.png"), dpi=150)
assert all(M[i, s] == loops[i][2] for i in range(len(loops)) for s in range(3))

# --- coasting times (panel a: Model Gc; panel b: EX-GC25 kinks, separate units and measure)
rows = []
for lab, retained in (("EX-GC23 strong, transport off", .739108), ("EX-GC23 strong, transport on", .970598),
                      ("EX-GC23 weak, transport on", .998139)):       # V(20)/V(0), results/EX-GC23 note table
    rows.append((lab + "\n(measured speed retention)", -20 / np.log(retained), "C0"))
d = json.load(open(os.path.join(ROOT, "results", "EX-GC24", "diagnostics.json")))
st = {r["S"]: r for r in d["steady"] if r["v"] == 0.01}
wh = {r["S"]: r for r in d["weighted_harmonic"]}
for S, lab in ((5.32, "EX-GC24A strong"), (0.266, "EX-GC24A weak")):
    b_total = -st[S]["total"] / 0.01
    rows.append((lab + "\n(derived m / b_total)", wh[S]["m"] / b_total, "C3"))
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4), gridspec_kw=dict(width_ratios=[1.4, 1]))
y = np.arange(len(rows))[::-1]
ax.barh(y, [r[1] for r in rows], color=[r[2] for r in rows])
for yi, r in zip(y, rows):
    ax.text(r[1] * 1.3, yi, f"{r[1]:.3g}", va="center", fontsize=8)
ax.set_yticks(y, [r[0] for r in rows], fontsize=7)
ax.set_xscale("log")
ax.set_xlim(1e-3, 1e6)
ax.axvline(20, color="k", ls=":", lw=1)
ax.text(22, y[0] + 0.45, "length of the EX-GC23 runs (20)", fontsize=7)
ax.set_xlabel("characteristic coasting time (Model Gc time units)")
ax.set_title("(a) Model Gc: imposed inertia, no head-on drag (blue)\nvs medium inertia against total drag (red)", fontsize=9)
r25 = json.load(open(os.path.join(ROOT, "results", "EX-GC25", "EX-GC25_runs.json")))
for k, mk in (("A", "o-"), ("B", "s-")):
    rr = [r for r in r25["coast"]["rows"] if r["kink"] == k]
    bx.plot([r["h"] for r in rr], [abs(r["loss"]) for r in rr], mk, label=f"kink {k}")
bx.set_xscale("log"); bx.set_yscale("log"); bx.invert_xaxis()
bx.set_xlabel("lattice spacing h (finer to the right)")
bx.set_ylabel("|fractional speed change|, early vs late fit")
bx.set_title("(b) EX-GC25, 1+1 D kinks at u = 0.5, t = 0-200\n(different model and units; not comparable to (a))", fontsize=9)
bx.legend(fontsize=8)
fig.suptitle("EXPLORATORY", fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "wp_coasting_times.png"), dpi=150)
for r in rows:
    print(r[0].replace("\n", " "), f"{r[1]:.4g}")
