"""Schematic of the EX-GC51 declared model (illustration only, no data). Our construction; not LaViolette's Model G."""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

fig = plt.figure(figsize=(14, 7.2), dpi=150)
ax = fig.add_axes([0.02, 0.06, 0.52, 0.86]); ax.set_aspect("equal"); ax.axis("off")
n = 13; c = n // 2; xs, ys = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
rr = np.hypot(xs - c, ys - c) + 0.8; N = np.exp(-1.2 / rr)
# bonds coloured by coincidence weight N_v N_w
for i in range(n):
    for j in range(n):
        for di, dj in ((1, 0), (0, 1)):
            k, l = i + di, j + dj
            if k < n and l < n:
                w = N[i, j] * N[k, l]; ax.plot([i, k], [j, l], color=plt.cm.viridis((w - 0.4) / 0.6), lw=1.6, zorder=1)
sc = ax.scatter(xs, ys, c=N, cmap="viridis", vmin=0.4, vmax=1.0, s=70, zorder=2, edgecolors="k", linewidths=0.3)
# source (frozen fermion) and probe packet
ax.add_patch(Circle((c, c), 1.3, color="#d62728", alpha=0.35, zorder=3)); ax.text(c, c - 2.0, "source fermion\n(charge q = ∂H/∂N)", ha="center", va="top", fontsize=9, color="#a00", bbox=dict(fc="white", ec="none", alpha=0.8), zorder=5)
px, py = c + 4.5, c + 3.5
ax.add_patch(Circle((px, py), 0.9, color="#1f77b4", alpha=0.45, zorder=3)); ax.annotate("", xy=(px - 1.3, py - 1.0), xytext=(px, py), arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=2), zorder=4)
ax.text(px - 1.5, py + 1.3, "probe fermion falls", fontsize=9, color="#1f77b4", bbox=dict(fc="white", ec="none", alpha=0.8), zorder=5)
# light ray bending past the source
t = np.linspace(-0.5, n - 0.5, 200); yb = c - 3.6 + 0.9 / (1 + ((t - c) / 2.5) ** 2) * 0 ; yb = c - 3.6 + 0.45 * (1 + np.tanh((t - c) / 2.0))
ax.plot(t, yb, color="orange", lw=2.5, zorder=4); ax.text(0, c - 4.3, "light ray (bent toward the source, delayed)", fontsize=9, color="#c60", bbox=dict(fc="white", ec="none", alpha=0.8))
cb = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.04, pad=0.01, shrink=0.7); cb.set_label("local firing rate N (clock rate)", fontsize=9)
ax.set_title("One cubic lattice, one field N(x)  (2D slice; well depth exaggerated)", fontsize=12)

def box(y, h, title, lines, col):
    bx = fig.add_axes([0.57, y, 0.41, h]); bx.axis("off")
    bx.add_patch(FancyBboxPatch((0, 0), 1, 1, boxstyle="round,pad=0.02", fc=col, ec="0.3", transform=bx.transAxes))
    bx.text(0.03, 0.88, title, fontsize=11, weight="bold", va="top", transform=bx.transAxes)
    bx.text(0.03, 0.62, "\n".join(lines), fontsize=9.5, va="top", transform=bx.transAxes, linespacing=1.5)
box(0.70, 0.22, "1. Clocks and signals (37B event counting)",
    ["Each site ticks at rate N:  dτ = N dt   (a one-site event)",
     "A signal between neighbours needs both to fire: weight N_v N_w",
     "Light:  ∂²u/∂τ² = Σ N_x N_y (u_y − u_x)"], "#fff4e0")
box(0.47, 0.21, "2. Matter: one Wilson–Dirac fermion (48)",
    ["Hops and Wilson bonds weighted N_v N_w;  mass m β on a site weighted N_v",
     "⇒ rest energy ∝ N, limiting speed N² (same as light)",
     "Gravitational charge  q_x = ∂⟨H⟩/∂N_x  (≈ energy for a packet at rest)"], "#e8f0fb")
box(0.24, 0.21, "3. Field: coincidences pass on light's slowness (49)",
    ["Transported quantity  f(N) = −1/(2N²)  (coordinate delay per hop)",
     "(1/N)∂_t[(1/N)∂_t f] = Σ N_x N_y (f_y − f_x) − 4πG q",
     "⇒ same light cone as light; static exterior N = e^(−U), β = 1"], "#e9f7e9")
box(0.03, 0.19, "4. Back-reaction: nothing else put in",
    ["Matter moves in the evolving N; N is sourced by matter's q",
     "Tests: 51A cone, 51B Newton and β, 51C fall, 51D action = reaction, 51E bending"], "#f3e9f7")
fig.text(0.775, 0.005, "EX-GC51 declared model (schematic, EXPLORATORY). Our construction; not LaViolette's Model G.", ha="center", fontsize=8.5, color="0.35")
fig.savefig("paper/figures/exgc51_model_schematic.png", bbox_inches="tight", pad_inches=0.15); print("ok")
