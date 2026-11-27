"""EX-GC45 separate recomputation. (1) 45B Q_a/dE at the heavy-band minimum from eigenvalues only of an independently
built 7-band Bloch H (opposite phase convention, 4-point stencil in the hopping scale, no eigenvectors). (2) eta at
45C points from a dense projected two-body matrix at L = 10 (own heavy-band eigenvectors and pair basis), compared with
the driver's separable solver at the same L. Shared with the driver: the model definition, the U values in C.json.
EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys
from pathlib import Path
import numpy as np

R = Path(__file__).resolve().parents[1] / "results" / "EX-GC45"
sys.path.insert(0, str(Path(__file__).resolve().parent))


def H7(kx, ky, kz, lam0, m0, g, t):
    H = np.zeros((7, 7), complex)
    for a, k in enumerate((kx, ky, kz)):            # link a: mode 2a at v, mode 2a+1 at v+e_a; opposite phase convention
        H[2 * a + 1, 2 * a] = -t * np.exp(1j * k); H[2 * a, 2 * a + 1] = np.conj(H[2 * a + 1, 2 * a])
    for i in range(6):
        for j in range(6):
            if i != j:
                H[i, j] += lam0
    H[6, 6] = m0
    for i in range(6):
        H[6, i] = H[i, 6] = g
    return H


def heavy_eigs(L, lam0, m0, g, t, vecs=False):
    ks = 2 * np.pi * np.arange(L) / L
    E = np.zeros((L, L, L)); U6 = np.zeros((L, L, L), complex)
    for a, kx in enumerate(ks):
        for b, ky in enumerate(ks):
            for c, kz in enumerate(ks):
                w, v = np.linalg.eigh(H7(kx, ky, kz, lam0, m0, g, t)); E[a, b, c] = w[3]; U6[a, b, c] = v[6, 3]
    return (E, U6) if vecs else E


def part1(L=16, lam0=0.1, h=1e-4):
    out = {}; drv = json.loads((R / "B.json").read_text())["result"]
    for m0 in (0.3, 0.5, 0.7):
        for g in (0.025, 0.05, 0.1, 0.2):
            Es = {s: heavy_eigs(L, lam0, m0, g, 1 + s * h).min() for s in (-2, -1, 0, 1, 2)}
            dEdt = (-Es[2] + 8 * Es[1] - 8 * Es[-1] + Es[-2]) / (12 * h)
            r = dEdt / Es[0]; d = drv[f"{m0}_{g}"]["Qa_over_dE_min"]
            out[f"{m0}_{g}"] = dict(Qa_over_dE=r, driver=d, rel=abs(r / d - 1))
            print("R1", m0, g, out[f"{m0}_{g}"], flush=True)
    return out


def pair_dense(L, lam0, m0, g, U, t):
    E, U6 = heavy_eigs(L, lam0, m0, g, t, vecs=True); N = L ** 3
    idx = np.arange(L); neg = (-idx) % L
    rows = []; D = []
    for a in range(L):
        for b in range(L):
            for c in range(L):
                D.append(E[a, b, c] + E[neg[a], neg[b], neg[c]]); rows.append(U6[a, b, c] * U6[neg[a], neg[b], neg[c]])
    D = np.array(D); f = np.array(rows)
    H = np.diag(D) - (U / N) * np.outer(np.conj(f), f)      # <k'|V|k> = -(U/N) conj(f_k') f_k
    return np.linalg.eigvalsh(H)[0], D.min()


def part2(L=10, lam0=0.1, h=1e-5, keys=("0.3_0.2_0.15", "0.5_0.1_0.15", "0.7_0.05_0.15")):
    import exgc45 as X
    C = json.loads((R / "C.json").read_text())["result"]; out = {}
    for key in keys:
        c = C[key]; m0, g, U = c["m0"], c["g"], c["U"]
        p = pair_dense(L, lam0, m0, g, U, 1 + h); m = pair_dense(L, lam0, m0, g, U, 1 - h); z = pair_dense(L, lam0, m0, g, U, 1.0)
        dQ = (p[0] - m[0]) / (2 * h) - (p[1] - m[1]) / (2 * h); eta = dQ / (z[1] - z[0])
        sep = {}
        for s in (1, -1):
            e, u, _ = X.eig(L, lam0, m0=m0, g=g, t=(1 + s * h,) * 3); sep[s] = X.pair_solve(e[..., 3], u[..., :, 3], U)
        e, u, _ = X.eig(L, lam0, m0=m0, g=g); s0 = X.pair_solve(e[..., 3], u[..., :, 3], U)
        dQs = (sep[1]["E0"] - sep[-1]["E0"]) / (2 * h) - (sep[1]["Dmin"] - sep[-1]["Dmin"]) / (2 * h)
        out[key] = dict(L=L, eta_dense=eta, eta_separable_sameL=dQs / s0["Eb"], eta_driver_L32=c["eta"],
                        diff_sameL=abs(eta - dQs / s0["Eb"]))
        print("R2", key, out[key], flush=True)
    return out


if __name__ == "__main__":
    res = dict(part1=part1(), part2=part2(), platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (R / "recompute.json").write_text(json.dumps(res, indent=1, default=float) + "\n"); print("done")
