"""EX-GC47 separate recomputation at L = 6. Own 7-band Bloch build (exgc45_recompute.H7, opposite bond convention:
mode 2a+1 at v <-> mode 2a at v+e_a), own dense heavy-band pair problem at every total momentum K, own gauge rule
(a_K = <on-site pair Bloch sum|Psi_K> real positive), real-space two-body amplitudes, real-space lapse/strain operators
and contact term, transition densities F(q) = <Psi_q|o(q)|Psi_0> for all 216 q (no symmetry fill), and the cross
energy V12 for B1/B1 at R = (3,0,0) and (2,1,0). Compared with the driver's form-factor result at the same L and U.
Shared: model definition, U, the kernel formula. EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, itertools
from pathlib import Path
import numpy as np
import scipy.sparse as sp
sys.path.insert(0, str(Path(__file__).resolve().parent))
from exgc45_recompute import H7

R = Path(__file__).resolve().parents[1] / "results" / "EX-GC47"
LAM0, C2 = 0.1, 1.001


def main(L=6, G4pi=0.01):
    A = json.loads((R / "A.json").read_text())["result"]["B1"]; m0, g, U = A["m0"], A["g"], A["U"]; N = L ** 3
    grid = list(itertools.product(range(L), repeat=3)); ks = [2 * np.pi * np.array(n) / L for n in grid]
    E = np.zeros(N); Uv = np.zeros((N, 7), complex)
    for i, k in enumerate(ks):
        w, v = np.linalg.eigh(H7(*k, LAM0, m0, g, 1.0)); E[i] = w[3]; Uv[i] = v[:, 3]
    ix = {n: i for i, n in enumerate(grid)}; add = lambda a, b: ix[tuple((np.array(grid[a]) + np.array(grid[b])) % L)]
    mk = [ix[tuple((-np.array(n)) % L)] for n in grid]
    X = np.array(grid)
    def state(Kn):
        K = ix[Kn]; i1 = [add(i, K) for i in range(N)]
        D = E[i1] + E[mk]; Af = np.conj(Uv[i1, 6]) * np.conj(Uv[mk, 6])
        H = np.diag(D).astype(complex) - (U / N) * np.outer(Af, np.conj(Af))
        w, v = np.linalg.eigh(H); phi = v[:, 0]
        a = np.sum(phi * Uv[i1, 6] * Uv[mk, 6]); phi = phi * np.conj(a) / abs(a)                 # gauge: a_K > 0
        # real-space amplitude A[(x1,i),(x2,j)] = (1/N) sum_k phi(k) u_i(k+K) e^{i(k+K)x1} u_j(-k) e^{-ik x2}
        kK = np.array([ks[i] for i in i1]); km = np.array([ks[i] for i in mk])
        P1 = (np.exp(1j * X @ kK.T)[:, None, :] * Uv[i1].T[None, :, :]).reshape(N * 7, N) / np.sqrt(N)
        P2 = (np.exp(1j * X @ km.T)[:, None, :] * Uv[mk].T[None, :, :]).reshape(N * 7, N) / np.sqrt(N)
        return w[0], (P1 * phi[None, :]) @ P2.T
    def ops(q):
        rows, cols, vN = [], [], []; site = LAM0 * (np.ones((7, 7)) - np.eye(7)); site[6, :] = g; site[:, 6] = g; site[6, 6] = m0
        vS = []
        for v, x in enumerate(grid):
            ph = np.exp(1j * np.dot(q, x))
            for i in range(7):
                for j in range(7):
                    if site[i, j] != 0:
                        rows.append(v * 7 + i); cols.append(v * 7 + j); vN.append(ph * site[i, j]); vS.append(0.0)
            for a in range(3):
                y = list(x); y[a] = (y[a] + 1) % L; w = ix[tuple(y)]; wt = (ph + np.exp(1j * np.dot(q, y))) / 2
                for (r, c) in ((v * 7 + 2 * a + 1, w * 7 + 2 * a), (w * 7 + 2 * a, v * 7 + 2 * a + 1)):
                    rows.append(r); cols.append(c); vN.append(-wt); vS.append(-wt)
        return sp.csr_matrix((vN, (rows, cols)), shape=(N * 7, N * 7)), sp.csr_matrix((vS, (rows, cols)), shape=(N * 7, N * 7))
    E0, A0 = state((0, 0, 0)); F = {}
    hv = [v * 7 + 6 for v in range(N)]
    for n in grid:
        q = 2 * np.pi * np.array(n) / L; Eq, Aq = state(n); ON, Oa = ops(q); out = {}
        for lab, O in (("N", ON), ("a", Oa)):
            out[lab] = np.sum(np.conj(Aq) * (O @ A0)) + np.sum(np.conj(Aq) * (O @ A0.T).T)
        ph = np.exp(1j * X @ q)
        out["V"] = -U * np.sum(ph * np.conj(Aq[hv, hv]) * A0[hv, hv]); out["N"] += out["V"]; F[n] = out
    FN = np.zeros((L, L, L), complex); Fa = np.zeros_like(FN)
    for n in grid:
        FN[n], Fa[n] = F[n]["N"], F[n]["a"]
    qq = 2 * np.pi * np.fft.fftfreq(L); Q = np.meshgrid(qq, qq, qq, indexing="ij"); qt2 = sum(2 * (1 - np.cos(x)) for x in Q)
    Kq = np.zeros_like(qt2); Kq[qt2 > 0] = -G4pi / ((1 + G4pi * C2) * qt2[qt2 > 0])
    Vt = Kq * ((FN + Fa) * FN + FN * Fa)
    res = dict(E0=E0, F0_N=FN[0, 0, 0].real, F0_a=Fa[0, 0, 0].real, Fimag_max=float(max(abs(FN.imag).max(), abs(Fa.imag).max())))
    import exgc47 as D
    S = D.setup(L, m0, g, U=U); b0, FNd, Fad, _ = D.pair_tables(S)
    res["dF_N_max"] = float(np.abs(FN - FNd).max()); res["dF_a_max"] = float(np.abs(Fa - Fad).max())
    for Rv in ((3, 0, 0), (2, 1, 0)):
        ph = np.exp(1j * (Q[0] * Rv[0] + Q[1] * Rv[1] + Q[2] * Rv[2]))
        Vr = float(np.sum(Vt * ph).real / N); Vd = float(D.V12(dict(N=FNd, a=Fad), dict(N=FNd, a=Fad), D.kernel(L, G4pi), Rv).real)
        res[f"V12_{Rv}"] = dict(recompute=Vr, driver=Vd, rel=abs(Vr / Vd - 1))
    print(res, flush=True)
    (R / "recompute.json").write_text(json.dumps(dict(result=res, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}"), indent=1, default=float) + "\n")


if __name__ == "__main__":
    main()
