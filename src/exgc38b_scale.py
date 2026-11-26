"""EX-GC38B descriptive diagnostic (after the registered P2 failed): does the Tolman miss shrink when the lapse bump
is made smoother in lattice units? Bump centre 10s, width 4s, eps 0.1, m = 0.1/s, L = 40s, compare delta beta with
2 pi eps [x f - int_0^x f] on 2.5s <= x <= 19.5s."""
import json, sys, time
import numpy as np, flint
from flint import arb_mat
import exgc38b as B

def beta(L, m, lapse, prec, nrow):
    flint.ctx.prec = prec
    X, P = B.ground_XP(L, m, lapse)
    idx = range(L, 2 * L)
    Xr = arb_mat([[X[i, k] for k in idx] for i in idx]); Pr = arb_mat([[P[i, k] for k in idx] for i in idx])
    Nk, _, _, _ = B.ent_hamiltonian(Xr, Pr)
    return np.array([float(Nk[i, i].mid()) for i in range(nrow)])

out = {}
for s in [int(a) for a in sys.argv[1:]]:
    t = time.process_time()
    L, m, nrow = 40 * s, 0.1 / s, 20 * s
    xs = np.arange(-L + 1, L + 1) - 0.5
    c, w = 10.0 * s, 4.0 * s
    f = np.where(xs > 0, np.exp(-(xs - c) ** 2 / (2 * w * w)), 0.0)
    b0 = beta(L, m, np.ones(2 * L), 1000 * s, nrow)
    b1 = beta(L, m, 1 + 0.1 * f, 1000 * s, nrow)
    x = np.arange(1, nrow + 1) - 0.5
    xx = np.linspace(0, x.max(), 40001); ff = np.exp(-(xx - c) ** 2 / (2 * w * w))
    If = np.concatenate([[0], np.cumsum((ff[1:] + ff[:-1]) / 2 * np.diff(xx))])
    pred = 2 * np.pi * 0.1 * (x * np.exp(-(x - c) ** 2 / (2 * w * w)) - np.interp(x, xx, If))
    d = b1 - b0
    sel = (x >= 2.5 * s) & (x <= 19.5 * s)
    rel = float(np.sqrt(np.mean((d - pred)[sel] ** 2)) / np.sqrt(np.mean(pred[sel] ** 2)))
    slope = float(np.sum(d[sel] * pred[sel]) / np.sum(pred[sel] ** 2))
    out[s] = dict(L=L, m=m, x=x.tolist(), delta_beta=d.tolist(), predicted=pred.tolist(), rel_rms_err=rel,
                  best_scale=slope, cpu_s=time.process_time() - t)
    print(s, rel, slope, out[s]["cpu_s"], flush=True)
    B.OUT.mkdir(parents=True, exist_ok=True)
    (B.OUT / f"scale_s{s}.json").write_text(json.dumps(out[s], indent=1) + "\n")
