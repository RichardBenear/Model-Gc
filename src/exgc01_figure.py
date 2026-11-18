"""EX-GC01 figure: weighted action/reaction M_a v_a vs M_b v_b, and speeds vs prediction."""
import json, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
o = json.load(open(os.path.join(ROOT, "results", "EX-GC01", "pairs.json")))
lab = {("ref", "ref"): "equal particles", ("ref", "dw0.2"): "b: D_W 0.2 (M +13%, μ −0.7%)",
       ("ref", "k0.4"): "b: k 0.4 (M +46%, μ +81%)"}
col = {("ref", "ref"): "0.4", ("ref", "dw0.2"): "#1f77b4", ("ref", "k0.4"): "#d62728"}
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), dpi=140)
for key in lab:
    P = [p for p in o["pairs"] if (p["a"], p["b"]) == key]
    d = [p["d"] for p in P]
    ax[0].plot(d, [p["momentum_ratio"] for p in P], "o-", color=col[key], label=lab[key])
    ax[0].axhline(P[0]["momentum_ratio_pred"], color=col[key], ls="--", lw=0.8)
    ax[1].plot(d, [p["va"] / p["pred_va"] for p in P], "o-", color=col[key], label=f"{lab[key]}: particle a")
    ax[1].plot(d, [p["vb"] / p["pred_vb"] for p in P], "s--", color=col[key], mfc="white", label="particle b")
ax[0].set(xlabel="separation d", ylabel="(M_a v_a) / (M_b v_b)", title="Weighted action and reaction (dashed: μ_a/μ_b)")
ax[0].legend(fontsize=8)
ax[1].set(xlabel="separation d", ylabel="v / (μ_i M_j / 4π d²)", title="Each particle's speed vs linear prediction")
ax[1].legend(fontsize=6.5, ncol=1)
fig.tight_layout()
fig.savefig(os.path.join(ROOT, "paper", "figures", "EX-GC01_action_reaction.png"))
