"""EX-GC46 separate recomputation: independently built 7-band H(k) (opposite phase convention, own loops), a dense
four-band two-body matrix at L = 6 (own pair basis, no rank-1 solver), eta by central difference in the hopping scale,
compared with the driver's rank-1 solver at the same L and U. Shared: the model definition and the U values from 45's C.json.
EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
from exgc45_recompute import H7
import exgc45 as X, exgc46 as Y

R = Path(__file__).resolve().parents[1] / "results" / "EX-GC46"


def dense4(L, lam0, m0, g, U, t):
    ks = 2 * np.pi * np.arange(L) / L; N = L ** 3
    E = np.zeros((N, 4)); H6 = np.zeros((N, 4), complex); idx = {}
    n = 0
    for a in range(L):
        for b in range(L):
            for c in range(L):
                w, v = np.linalg.eigh(H7(ks[a], ks[b], ks[c], lam0, m0, g, t)); E[n] = w[3:]; H6[n] = v[6, 3:]; idx[(a, b, c)] = n; n += 1
    mi = np.array([idx[((-a) % L, (-b) % L, (-c) % L)] for a, b, c in sorted(idx, key=idx.get)])
    D = (E[:, :, None] + E[mi][:, None, :]).reshape(-1); P = (np.conj(H6)[:, :, None] * np.conj(H6[mi])[:, None, :]).reshape(-1)
    Hm = np.diag(D).astype(complex) - (U / N) * np.outer(np.conj(P), P)
    return np.linalg.eigvalsh(Hm)[0], D.min()


def main(L=6, lam0=0.1, h=1e-5):
    C = json.loads((R.parent / "EX-GC45" / "C.json").read_text())["result"]; out = {}
    for key in ("0.3_0.2_0.15", "0.5_0.1_0.15", "0.7_0.05_0.15"):
        c = C[key]; m0, g, U = c["m0"], c["g"], c["U"]
        p = dense4(L, lam0, m0, g, U, 1 + h); m = dense4(L, lam0, m0, g, U, 1 - h); z = dense4(L, lam0, m0, g, U, 1.0)
        eta = ((p[0] - m[0]) - (p[1] - m[1])) / (2 * h) / (z[1] - z[0])
        e, u, _ = X.eig(L, lam0, m0=m0, g=g); r0 = Y.solve(e, u, U)
        rr = {}
        for s in (1, -1):
            es, us, _ = X.eig(L, lam0, m0=m0, g=g, t=(1 + s * h,) * 3); rr[s] = Y.solve(es, us, U)
        etaY = ((rr[1]["E0"] - rr[-1]["E0"]) - (rr[1]["Dmin"] - rr[-1]["Dmin"])) / (2 * h) / r0["Eb"]
        out[key] = dict(L=L, eta_dense=eta, eta_driver_sameL=etaY, diff=abs(eta - etaY), Eb_dense=z[1] - z[0], Eb_driver=r0["Eb"])
        print(key, out[key], flush=True)
    (R / "recompute.json").write_text(json.dumps(dict(result=out, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}"), indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
