"""EX-GC19 figure: core velocity, mass and knot/core positions for the flow-field runs."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
D = os.path.join(ROOT, "results", "EX-GC19")
R = {c: json.load(open(os.path.join(D, c + ".json"))) for c in
     ("a_drag", "b_nodrag", "c_kick_nodrag", "d_kick_drag", "e_nodrag_kappa0", "f_K04_nodrag", "g_knot_nodrag")}
lab = dict(a_drag="fall, drag on", b_nodrag="fall, drag off", c_kick_nodrag="kick, drag off", d_kick_drag="kick, drag on",
           e_nodrag_kappa0="fall, drag off, no cohesion", f_K04_nodrag="fall, drag off, K 0.4",
           g_knot_nodrag="knot particle, drag off")
col = dict(a_drag="#1f77b4", b_nodrag="#d62728", c_kick_nodrag="#2ca02c", d_kick_drag="#9467bd",
           e_nodrag_kappa0="#7f7f7f", f_K04_nodrag="#ff7f0e", g_knot_nodrag="#8c564b")
fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
t = np.linspace(0, 24, 100)
for c in ("a_drag", "b_nodrag", "f_K04_nodrag", "c_kick_nodrag", "d_kick_drag"):
    r = R[c]["rows"]
    ax[0].plot([x["t"] for x in r], [-x["Vcore"][2] for x in r], color=col[c], label=lab[c])
ax[0].plot(t, 0.049 / 20 * t, "k:", lw=.8, label="no wake: push/mass, 0.00245 t")
ax[0].plot(t, 0.05 * np.ones_like(t), "k--", lw=.8, label="perfect coasting")
ax[0].set_xlabel("t")
ax[0].set_ylabel("downward core velocity")
ax[0].set_title("Velocity of the particle")
ax[0].legend(fontsize=7)
for c in R:
    r = R[c]["rows"]
    ax[1].plot([x["t"] for x in r], [x["M"] / r[0]["M"] for x in r], color=col[c], label=lab[c])
ax[1].set_xlabel("t")
ax[1].set_ylabel("M / M(0)")
ax[1].set_title("Particle strength")
ax[1].legend(fontsize=7)
r = R["g_knot_nodrag"]["rows"]
ax[2].plot([x["t"] for x in r], [x["Wcen"][2] for x in r], color="#8c564b", label="core (W)")
ax[2].plot([x["t"] for x in r], [x["Kcen"][2] for x in r], color="#8c564b", alpha=.4, lw=3, label="knot")
ax[2].set_xlabel("t")
ax[2].set_ylabel("z position")
ax[2].set_title("Knot particle: does the knot keep up?")
ax[2].legend(fontsize=7)
fig.suptitle("EX-GC19 (EXPLORATORY): flow (momentum) field for the Gc particle, E = 0.02, tau = 20 — our construction",
             fontsize=10)
fig.tight_layout()
out = os.path.join(ROOT, "paper", "figures", "exgc19_flow.png")
fig.savefig(out, dpi=130)
print(out)
