"""EX-GC02 figure: drag, effective mass and mobility vs tau_W; time response with and without inertia."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
o = json.load(open(os.path.join(ROOT, "results", "EX-GC02", "results.json")))
col = {"ref": "0.2", "dw0.2": "#1f77b4", "k0.4": "#d62728"}
lab = {"ref": "reference", "dw0.2": "D_W 0.2", "k0.4": "k 0.4"}
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), dpi=140)
for n, rows in o["scan"].items():
    t = [r["tauW"] for r in rows]
    ax[0].plot(t, [r["gamma"] for r in rows], "o-", color=col[n], label=f"drag γ, {lab[n]}")
    ax[0].plot(t, [r["m"] / 40 for r in rows], "s--", color=col[n], mfc="white", label=f"mass m/40, {lab[n]}")
    ok = [r for r in rows if r["gamma"] > 0]           # gamma < 0: resting particle unstable (self-propelled)
    ax[1].plot([r["tauW"] for r in ok], [r["mu"] for r in ok], "o-", color=col[n], label=lab[n])
ax[0].axhline(0, color="0.6", lw=0.6)
ax[0].set(xlabel="core response time τ_W", ylabel="drag γ (γ = 1 at τ = 1), mass m/40", title="Drag falls, inertia grows")
ax[0].legend(fontsize=7)
ax[1].axhline(0, color="0.6", lw=0.6)
ax[1].set(xlabel="core response time τ_W", ylabel="G mobility μ", title="Mobility stops being common near threshold")
ax[1].set_ylim(-0.6, 0.02)
ax[1].text(16.5, -0.45, "k 0.4 self-propels\nbeyond τ_W ≈ 17.5", color="#d62728", fontsize=8)
ax[1].legend(fontsize=8)
r = o["response"]["25"]
ax[2].plot(r["t"], r["X"], color="k", lw=2, label="linear simulation (l = 1)")
ax[2].plot(r["t"], r["pred_with_mass"], "--", color="#d62728", label=f"m X'' + γ X' = f  (m = {r['m']:.1f})")
ax[2].plot(r["t"], r["pred_no_mass"], ":", color="#1f77b4", label="γ X' = f  (no inertia)")
ax[2].set(xlabel="time after the G gradient is switched on", ylabel="displacement X (E = 1e-3)",
          title="Response at τ_W = 25: velocity builds up")
ax[2].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC02_inertia.png"))
