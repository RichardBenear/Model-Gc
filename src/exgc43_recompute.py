"""EX-GC43 separate recomputation of the packet charges Q = 2 E_link + E_site. Independent of band_data/packet_QE:
builds the Bloch Hamiltonian H(t) = t*H_link + lam0*J from scratch and gets E_link = t dE/dt by a finite difference of
band eigenvalues (no eigenvectors, no expectation values). Packets use sorted band energies, so this path is exact
for non-degenerate or flat bands. Shared dependency: the packet definition (Gaussian weights, centres) and the plan."""
import json
from pathlib import Path
import numpy as np
R = Path(__file__).resolve().parents[1] / "results" / "EX-GC43"
M = 6
def bands(L, lam0, t):
    k = 2 * np.pi * np.fft.fftfreq(L); KK = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    H = np.zeros((L, L, L, M, M), complex)
    for a in range(3):
        H[..., 2 * a, 2 * a + 1] = -t * np.exp(1j * KK[..., a]); H[..., 2 * a + 1, 2 * a] = -t * np.exp(-1j * KK[..., a])
    H += lam0 * (np.ones((M, M)) - np.eye(M))
    return np.linalg.eigvalsh(H), KK
def weights(KK, k0, s):
    dk = (KK - np.array(k0) + np.pi) % (2 * np.pi) - np.pi; w = np.exp(-s ** 2 * np.sum(dk ** 2, -1)); return w / w.sum()
saved = json.loads((R / "packets.json").read_text())["result"]; out = {}; L = 32; h = 1e-5
for lam0 in ("0.1", "0.2"):
    e0, KK = bands(L, float(lam0), 1.0); ep, _ = bands(L, float(lam0), 1 + h); em, _ = bands(L, float(lam0), 1 - h)
    el = (ep - em) / (2 * h)                                     # E_link per state = t dE/dt at t = 1
    sv = saved[lam0]; kc = sv["k_cond_min"]
    i_top = np.unravel_index(np.argmin(e0[..., 5]), e0[..., 5].shape); kt = [float(KK[i_top][a]) for a in range(3)]
    i_v = np.unravel_index(np.argmax(e0[..., 2]), e0[..., 2].shape); kv = [float(KK[i_v][a]) for a in range(3)]
    def pk(b, k0, s):
        w = weights(KK, k0, s); dE = np.sum(w * e0[..., b]); Q = np.sum(w * (e0[..., b] + el[..., b])); return dE, Q
    rows = {}
    for name, b, k0 in (("cond_k0", 3, [0, 0, 0]), ("cond_min", 3, kc), ("top_k0", 5, [0, 0, 0]), ("top_min", 5, kt)):
        for s in (1.0, 2.0, 4.0):
            dE, Q = pk(b, k0, s); key = f"{name}_s{s:g}"
            rows[key] = dict(ratio=Q / dE, saved=sv["packets"][key]["ratio"], diff=abs(Q / dE - sv["packets"][key]["ratio"]))
    (dc, qc), (dv, qv) = pk(3, kc, 2.0), pk(2, kv, 2.0)
    rows["pair_s2"] = dict(ratio=(qc - qv) / (dc - dv), saved=sv["packets"]["pair_s2"]["ratio"],
                           diff=abs((qc - qv) / (dc - dv) - sv["packets"]["pair_s2"]["ratio"]))
    out[lam0] = rows
    print(lam0, "max |diff|", max(r["diff"] for r in rows.values()))
(R / "recompute.json").write_text(json.dumps(out, indent=1) + "\n")
