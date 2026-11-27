"""EX-GC44 separate recomputation of the exciton's von Laue exponent w* = -d_link/<V>. Independent of the separable
solver: builds the projected two-body matrix densely at L = 12 from its own Bloch construction and takes d_link and <V>
from eigenvalue finite differences only (hopping scale t, coupling U). Uses the U values stored in D.json. Shared
dependencies: the projected model definition (bands 5 and 2, on-site density attraction) and the stored U."""
import json
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1] / "results" / "EX-GC44"
M = 6
def bands(L, lam0, t):
    k = 2 * np.pi * np.fft.fftfreq(L); KK = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    H = np.zeros((L, L, L, M, M), complex)
    for a in range(3):
        H[..., 2 * a, 2 * a + 1] = -t * np.exp(1j * KK[..., a]); H[..., 2 * a + 1, 2 * a] = -t * np.exp(-1j * KK[..., a])
    H += lam0 * (np.ones((M, M)) - np.eye(M))
    return np.linalg.eigh(H)
def E0(L, lam0, t, U):
    e, u = bands(L, lam0, t)
    ep, ev = e[..., 5].reshape(-1), e[..., 2].reshape(-1); up, uv = u[..., :, 5].reshape(-1, M), u[..., :, 2].reshape(-1, M)
    H = np.diag(ep - ev).astype(complex) - (U / L ** 3) * (up.conj() @ up.T) * (uv @ uv.conj().T)
    return np.linalg.eigvalsh(H)[0], (ep - ev).min()
D = json.loads((R / "D.json").read_text())["result"]; out = {}; L = 12; h = 1e-5
for lam0, r in D.items():
    for fr, row in r["rows"].items():
        U = row["U"]; e0, edge = E0(L, float(lam0), 1.0, U)
        dl = (E0(L, float(lam0), 1 + h, U)[0] - E0(L, float(lam0), 1 - h, U)[0]) / (2 * h) - \
             (E0(L, float(lam0), 1 + h, U)[1] - E0(L, float(lam0), 1 - h, U)[1]) / (2 * h)
        V = U * (E0(L, float(lam0), 1.0, U * (1 + h))[0] - E0(L, float(lam0), 1.0, U * (1 - h))[0]) / (2 * h * U)
        out[f"{lam0}_{fr}"] = dict(Eb_L12=float(edge - e0), Eb_L32=row["Eb"], w_star_L12=float(-dl / V), w_star_L32=row["w_star"])
        print(lam0, fr, out[f"{lam0}_{fr}"], flush=True)
(R / "recompute.json").write_text(json.dumps(out, indent=1) + "\n")
