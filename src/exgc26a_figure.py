"""EX-GC26A figure from results/EX-GC26A/*.json."""
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
R = os.path.join(ROOT, "results", "EX-GC26A")
ld = lambda f: json.load(open(os.path.join(R, f)))
rad = ld("radial.json")
st = {h: ld(f"static_{h}_5_hedgehog.json") for h in ("0.25", "0.2", "0.15")}
kick = ld("kick_0.2_5_6.json")
fig, ax = plt.subplots(1, 3, figsize=(14, 4))
hs = np.array([float(h) for h in st])
E = np.array([st[h]["E"] for h in st])
B = np.array([abs(st[h]["B"]) for h in st])
ax[0].plot(hs ** 2, E, "o", label="3D relaxed energy (box ±5)")
a, b = np.polyfit(hs[1:] ** 2, E[1:], 1)
ax[0].plot([0, 0.07], [b, a * 0.07 + b], "k--", lw=1, label=f"h² extrapolation (0.2, 0.15): {b:.4f}")
ax[0].axhline(rad["E"], color="C3", label=f"radial reference {rad['E']:.4f}")
ax[0].set_xlabel("h²"); ax[0].set_ylabel("E / 12π²"); ax[0].legend(fontsize=7)
ax[0].set_title("Rest energy (Skyrme units)", fontsize=9)
a2 = ax[0].twinx(); a2.plot(hs ** 2, B, "s", color="C2", ms=4); a2.set_ylabel("|degree| on the lattice", color="C2")
t = np.array([r["t"] for r in kick["rows"]])
H = np.array([r["H"] for r in kick["rows"]])
P = np.array([r["P"] for r in kick["rows"]])
Bk = np.array([abs(r["B"]) for r in kick["rows"]])
ax[1].plot(t, np.abs(H / H[0] - 1) + 1e-12, label="|H(t)/H(0) − 1|")
ax[1].plot(t, np.linalg.norm(P - P[0], axis=1) / kick["E0"] + 1e-12, label="|P(t) − P(0)| / E₀")
ax[1].axhline(1e-3, color="k", ls=":", lw=1); ax[1].axvline(3, color="grey", ls=":", lw=1)
ax[1].set_yscale("log"); ax[1].set_xlabel("t"); ax[1].legend(fontsize=7)
ax[1].set_title("Kick test (h = 0.2): conservation ledgers", fontsize=9)
ax[2].plot(t, Bk, "o-", ms=3)
ax[2].axhline(abs(kick["B0"]), color="k", ls=":", lw=1)
ax[2].set_ylim(abs(kick["B0"]) - 0.03, abs(kick["B0"]) + 0.03)
ax[2].set_xlabel("t"); ax[2].set_ylabel("|degree|")
ax[2].set_title(f"Degree during the kick; relaxed afterwards to E = {kick['recovered_E']:.5f} (E₀ = {kick['E0']:.5f})",
                fontsize=8)
fig.suptitle("EX-GC26A (EXPLORATORY): plain Skyrme model on a lattice, our construction of the tests", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "exgc26a_skyrmion.png"), dpi=150)
