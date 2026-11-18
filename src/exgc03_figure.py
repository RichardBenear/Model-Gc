"""EX-GC03 figure: gravitational push per unit mass (p/M) is nearly common; mobility is not; effect of slow G."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
o = json.load(open(os.path.join(ROOT, "results", "EX-GC03", "results.json")))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), dpi=140)
for cp, mk, col in (("GX", "o", "#1f77b4"), ("GY", "s", "#2ca02c")):
    T = o["per_type"][cp]; names = list(T)
    C = np.array([T[n]["c"] for n in names]); P = np.array([T[n]["p"] for n in names]); M = np.array([T[n]["M"] for n in names])
    pm = P / M; mu = -P / C.sum(1)
    ax[0].plot(M, pm / pm.mean(), mk, color=col, label=f"{cp}: push per unit M (p/M)")
    ax[0].plot(M, mu / mu.mean(), mk, color=col, mfc="white", label=f"{cp}: mobility μ (fall speed)")
    tg = np.logspace(0, 4.3, 40)
    sp = [np.std(-P / (C @ np.array([1, 1, 1, t]))) / abs(np.mean(-P / (C @ np.array([1, 1, 1, t])))) for t in tg]
    ax[1].semilogx(tg, sp, color=col, label=cp)
for i, n in enumerate(o["per_type"]["GX"]):
    T = o["per_type"]["GX"][n]
    ax[0].annotate(n, (T["M"], 0.45 + 0.0 * i), fontsize=7, rotation=90, color="0.4", ha="center")
ax[0].axhline(1, color="0.6", lw=0.6)
ax[0].set(xlabel="gravitational mass M (G consumption)", ylabel="relative to the mean over 7 types",
          title="Push ∝ M (within ~5%), but fall speed varies 35%")
ax[0].legend(fontsize=7)
ax[1].set(xlabel="G response time τ_G (others 1)", ylabel="spread of μ across types (std/mean)",
          title="Slow G (field drag) narrows the spread only when tuned")
ax[1].legend(fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC03_coupling.png"))
