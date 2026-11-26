"""EX-GC38B descriptive diagnostic: local vs row-summed vs link readouts of the lapse response (L = 40, m = 0.1, s = 1)."""
import json
import numpy as np, flint
from flint import arb_mat
import exgc38b as B

L, m = 40, 0.1
flint.ctx.prec = 1000
xs = np.arange(-L + 1, L + 1) - 0.5
f = np.where(xs > 0, np.exp(-(xs - 10) ** 2 / 32), 0.0)
mats = []
for lap in (np.ones(2 * L), 1 + 0.1 * f):
    X, P = B.ground_XP(L, m, lap); idx = range(L, 2 * L)
    Xr = arb_mat([[X[i, k] for k in idx] for i in idx]); Pr = arb_mat([[P[i, k] for k in idx] for i in idx])
    Nk, Mk, _, _ = B.ent_hamiltonian(Xr, Pr)
    mats.append((np.array([[float(Nk[i, k].mid()) for k in range(L)] for i in range(20)]),
                 np.array([[float(Mk[i, k].mid()) for k in range(L)] for i in range(21)])))
(N0, M0), (N1, M1) = mats
x = np.arange(1, 21) - 0.5
xx = np.linspace(0, 20, 20001); ff = np.exp(-(xx - 10) ** 2 / 32)
If = np.concatenate([[0], np.cumsum((ff[1:] + ff[:-1]) / 2 * np.diff(xx))])
pred = 2 * np.pi * 0.1 * (x * np.exp(-(x - 10) ** 2 / 32) - np.interp(x, xx, If))
out = dict(x=x.tolist(), predicted=pred.tolist(),
           diag=(np.diag(N1[:, :20]) - np.diag(N0[:, :20])).tolist(),
           rowsum=(N1.sum(1) - N0.sum(1)).tolist(),
           link_x=(x + 0.5).tolist(),
           link=(-(np.array([M1[i, i + 1] for i in range(20)]) - np.array([M0[i, i + 1] for i in range(20)]))).tolist())
(B.OUT / "rowsum.json").write_text(json.dumps(out, indent=1) + "\n")
print({k: np.round(v, 3).tolist() for k, v in out.items() if k != "link_x"})

# step lapse (N = 1 -> 1.1 for x > 5) and uniform-lapse consistency, L = 40, 700 bits
def diagN(L, m, lap, mode="physical", prec=700, n=15):
    flint.ctx.prec = prec
    X, P = B.ground_XP(L, m, lap, mode); idx = range(L, 2 * L)
    Xr = arb_mat([[X[i, k] for k in idx] for i in idx]); Pr = arb_mat([[P[i, k] for k in idx] for i in idx])
    Nk, _, _, _ = B.ent_hamiltonian(Xr, Pr)
    return np.array([float(Nk[i, i].mid()) for i in range(n)])
xs = np.arange(-40 + 1, 41) - 0.5
step = {}
for mm, mode in ((0.005, "conformal"), (0.1, "physical")):
    b0 = diagN(40, mm, np.ones(80), mode)
    step[f"m{mm}_{mode}"] = (diagN(40, mm, np.where(xs > 5, 1.1, 1.0), mode) - b0).tolist()
    if mm == 0.1:
        step["uniform_1p1_maxabs"] = float(np.max(np.abs(diagN(40, mm, 1.1 * np.ones(80)) - b0)))
step["predicted_beyond_step"] = np.pi
(B.OUT / "step.json").write_text(json.dumps(step, indent=1) + "\n")
print(step)
