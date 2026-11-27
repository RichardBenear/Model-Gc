"""EX-GC50 separate paths. C1: (L,1,1) supercell exact diagonalisation of the 40A network with dN = eps cos qx and/or
da = eta cos qx (homogeneous and 37B weights), chi from energy second differences with Richardson in eps, vs the
perturbative kernel. C2: boson-fermion partner check on an L = 4 real-space Wilson-Dirac lattice with random N(x), a(bond):
filled-sea energy + zero-point energy of one real oscillator per mode with stiffness h^2. Shared: numpy, exgc40a.bloch_h
(only for the uniform k-space reference in C1's comparison). EXPLORATORY. Our construction; not LaViolette's Model G.
"""
import json, time, platform
import numpy as np
import exgc50 as X
import exgc48 as D48

M6 = 6; J6 = np.ones((6, 6)) - np.eye(6)


def supercell_E(lam0, Nx, ax, weights, L=32):
    """Energy per site of the half-filled network; N and a vary along x only (supercell (L,1,1)), transverse k grid L x L."""
    S = len(Nx); dim = S * M6; base = np.zeros((dim, dim), complex)
    w = lambda i, j: (Nx[i] * Nx[j]) if weights == "37B" else (Nx[i] + Nx[j]) / 2
    for i in range(S):
        base[i * 6:(i + 1) * 6, i * 6:(i + 1) * 6] += Nx[i] * lam0 * J6
        j = (i + 1) % S; abar = (ax[i] + ax[j]) / 2; t = w(i, j) / abar
        base[i * 6 + 0, j * 6 + 1] += -t; base[j * 6 + 1, i * 6 + 0] += -t
    ks = 2 * np.pi * np.fft.fftfreq(L); E = 0.0
    yz = []
    for i in range(S):
        t = w(i, i) / ax[i]; yz.append(t)
    for ky in ks:
        for kz in ks:
            h = base.copy()
            for i in range(S):
                t = yz[i]
                for a, k in ((1, ky), (2, kz)):
                    p, m = i * 6 + 2 * a, i * 6 + 2 * a + 1
                    h[p, m] += -t * np.exp(1j * k); h[m, p] += -t * np.exp(-1j * k)
            e = np.linalg.eigvalsh(h); E += np.sum(e[:dim // 2])
    return E / (L * L * S)


def c1(lam0=0.1, L=32, ms=(1, 2, 4)):
    S0 = X.net_setup(L, lam0); out = {}; x = np.arange(L); one = np.ones(L)
    E0H = supercell_E(lam0, one, one, "H", L); E037 = supercell_E(lam0, one, one, "37B", L)
    for mm in ms:
        c = np.cos(2 * np.pi * mm * x / L); row = {}
        def chis(eps):
            r = {}
            EpN = supercell_E(lam0, 1 + eps * c, one, "H", L); EmN = supercell_E(lam0, 1 - eps * c, one, "H", L)
            r["NN"] = 2 * (EpN + EmN - 2 * E0H) / eps ** 2
            Epa = supercell_E(lam0, one, 1 + eps * c, "H", L); Ema = supercell_E(lam0, one, 1 - eps * c, "H", L)
            r["aa"] = 2 * (Epa + Ema - 2 * E0H) / eps ** 2
            Epp = supercell_E(lam0, 1 + eps * c, 1 + eps * c, "H", L); Epm = supercell_E(lam0, 1 + eps * c, 1 - eps * c, "H", L)
            Emp = supercell_E(lam0, 1 - eps * c, 1 + eps * c, "H", L); Emm = supercell_E(lam0, 1 - eps * c, 1 - eps * c, "H", L)
            r["Na"] = 2 * (Epp - Epm - Emp + Emm) / (4 * eps ** 2)
            Ep3 = supercell_E(lam0, 1 + eps * c, one, "37B", L); Em3 = supercell_E(lam0, 1 - eps * c, one, "37B", L)
            r["NN37"] = 2 * (Ep3 + Em3 - 2 * E037) / eps ** 2
            return r
        r1, r2 = chis(0.01), chis(0.02); rich = {k: (4 * r1[k] - r2[k]) / 3 for k in r1}
        pH = X.net_kernel(S0, (mm, 0, 0), "H"); p37 = X.net_kernel(S0, (mm, 0, 0), "37B")
        pert = dict(NN=pH["NN"], aa=pH["aa"], Na=pH["Na"], NN37=p37["NN"])
        for k in rich:
            d = rich[k] - pert[k]; row[k] = dict(supercell=rich[k], pert=pert[k], abs=float(d), rel=float(d / pert[k]))
        out[str(mm)] = row; print("C1 m", mm, {k: "%.3e/%.2e" % (v["abs"], v["rel"]) for k, v in row.items()}, flush=True)
    return out


def c2(L=4, m=0.3, r=1.0, draws=5, seed=20261127):
    rng = np.random.default_rng(seed); Ns = L ** 3; G = D48.GAM; al = G[1:4]; be = G[4]
    idx = np.arange(Ns).reshape(L, L, L); X3 = np.stack(np.meshgrid(*[np.arange(L)] * 3, indexing="ij"), -1).reshape(-1, 3)
    def smooth():
        f = np.zeros(Ns)
        for _ in range(4):
            kv = rng.integers(-1, 2, 3) * 2 * np.pi / L; f += rng.normal() * np.cos(X3 @ kv + rng.uniform(0, 2 * np.pi))
        return 1 + 0.2 * f / max(1e-12, np.abs(f).max())
    out = []
    for d in range(draws):
        N = smooth(); A = [smooth() for _ in range(3)]       # a on bonds along axis a, indexed by the bond's lower site
        h = np.zeros((4 * Ns, 4 * Ns), complex)
        for v in range(Ns):
            h[4 * v:4 * v + 4, 4 * v:4 * v + 4] += N[v] * m * be
            for a in range(3):
                c = X3[v].copy(); c[a] = (c[a] + 1) % L; w = idx[tuple(c)]
                wt = (N[v] + N[w]) / 2 / A[a][v]
                blk = lambda i, j, B: h.__setitem__((slice(4 * i, 4 * i + 4), slice(4 * j, 4 * j + 4)), h[4 * i:4 * i + 4, 4 * j:4 * j + 4] + wt * B)
                blk(v, w, al[a] / 2j); blk(w, v, -al[a] / 2j)
                blk(v, v, r / 2 * be); blk(w, w, r / 2 * be); blk(v, w, -r / 2 * be); blk(w, v, -r / 2 * be)
        e = np.linalg.eigvalsh(h); EF = float(np.sum(e[e < 0])); EB = float(0.5 * np.sum(np.abs(e)))
        sym = float(np.max(np.abs(np.sort(e) + np.sort(e)[::-1])))
        out.append(dict(E_F=EF, E_B=EB, sum=EF + EB, rel=(EF + EB) / abs(EF), spectral_asym=sym, herm=float(np.abs(h - h.conj().T).max())))
        print("C2", d, out[-1], flush=True)
    return out


if __name__ == "__main__":
    t0 = time.time(); r2 = c2(); r1 = c1()
    rec = dict(study="EX-GC50", part="checks", label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(),
               wall_s=time.time() - t0, result=dict(C1=r1, C2=r2))
    (X.OUT / "checks.json").write_text(json.dumps(rec, indent=1, default=float))
