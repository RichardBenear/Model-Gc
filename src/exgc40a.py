"""EX-GC40A: entanglement equilibrium on a fermionic pair network.

Each site of a periodic cubic lattice carries one fermion mode per link end (kpl modes per end). The pair term joins
the two end modes of each link; at half filling and lam = 0 every link is one maximally entangled pair. lam couples
all modes of a site, spreading entanglement beyond single links.
Plan: spec/EX-GC40_plan.md. Our construction; not LaViolette's Model G and not the Chain Theory construction.
"""
import json, sys, time
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC40A"
LAMS = (0.0, 0.1, 0.2, 0.3, 0.5)


# ---------------------------------------------------------------- model
def nmodes(kpl):
    return 6 * kpl


def mode(d, s, kpl):
    """d = 0..5 for +x,-x,+y,-y,+z,-z; s = copy index."""
    return d * kpl + s


def bloch_h(L, lam, kpl=1):
    M = nmodes(kpl)
    k = 2 * np.pi * np.fft.fftfreq(L)
    K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    H = np.zeros((L, L, L, M, M), complex)
    for a in range(3):
        ph = np.exp(1j * K[..., a])
        for s in range(kpl):
            p, m = mode(2 * a, s, kpl), mode(2 * a + 1, s, kpl)
            H[..., p, m] = -ph          # (v,+a) <- (v+a,-a)
            H[..., m, p] = -np.conj(ph)
    H += lam * (np.ones((M, M)) - np.eye(M))
    return H, K


def bands(L, lam, kpl=1):
    H, K = bloch_h(L, lam, kpl)
    e, u = np.linalg.eigh(H)
    return e, u, K


def gap(e):
    M = e.shape[-1]
    return float(e[..., M // 2].min() - e[..., M // 2 - 1].max())


def P_real(e, u):
    """Pr[r] with P_{(x,m),(y,m')} = Pr[x-y][m,m'] (occupied projector, half filling)."""
    M = e.shape[-1]
    uo = u[..., :M // 2]
    Pk = np.einsum("...mi,...ni->...mn", uo, uo.conj())
    return np.fft.ifftn(Pk, axes=(0, 1, 2))


def block_coords(l, origin):
    r = np.arange(l)
    g = np.stack(np.meshgrid(r, r, r, indexing="ij"), -1).reshape(-1, 3)
    return g + np.asarray(origin)


def block_P(Pr, coords, L):
    M = Pr.shape[-1]
    d = (coords[:, None, :] - coords[None, :, :]) % L
    blk = Pr[d[..., 0], d[..., 1], d[..., 2]]            # (n, n, M, M)
    n = len(coords)
    return blk.transpose(0, 2, 1, 3).reshape(n * M, n * M)


def entropy(PB):
    p = np.linalg.eigvalsh(0.5 * (PB + PB.conj().T))
    p = np.clip(p, 1e-300, 1 - 1e-16)
    q = 1 - p
    return float(-np.sum(p * np.log(p) + np.where(q > 0, q * np.log(np.maximum(q, 1e-300)), 0.0)))


# ---------------------------------------------------------------- real space (dense)
def real_h(L, lam, kpl=1, rewire=None):
    M = nmodes(kpl)
    N = L ** 3
    idx = lambda x, m: ((x[0] % L) * L * L + (x[1] % L) * L + (x[2] % L)) * M + m
    h = np.zeros((N * M, N * M))
    for x0 in range(L):
        for x1 in range(L):
            for x2 in range(L):
                x = (x0, x1, x2)
                for a in range(3):
                    y = list(x); y[a] += 1
                    for s in range(kpl):
                        i, j = idx(x, mode(2 * a, s, kpl)), idx(y, mode(2 * a + 1, s, kpl))
                        h[i, j] = h[j, i] = -1.0
                base = idx(x, 0)
                h[base:base + M, base:base + M] += lam * (np.ones((M, M)) - np.eye(M))
    if rewire is not None:
        for (i, j) in rewire["remove"]:
            ii, jj = idx(*i), idx(*j)
            h[ii, jj] = h[jj, ii] = 0.0
        for (i, j) in rewire["add"]:
            ii, jj = idx(*i), idx(*j)
            h[ii, jj] = h[jj, ii] = -1.0
    return h, idx


def dense_P(h):
    e, v = np.linalg.eigh(h)
    n = h.shape[0] // 2
    return v[:, :n] @ v[:, :n].T, float(e[n] - e[n - 1])


def dense_block(P, coords, idx, M):
    ids = np.array([idx(tuple(c), m) for c in coords for m in range(M)])
    return P[np.ix_(ids, ids)]


# ---------------------------------------------------------------- rewiring
def rewire_spec(a0, kpl=1):
    """Two adjacent sites a, c = a + x on a block's top (+z) face; their +z links to b = a+z, d = c+z are rewired
    to (a,+z)-(c,+z) and (b,-z)-(d,-z). Degrees and modes are unchanged; two cut links are removed."""
    a = np.array(a0); c = a + [1, 0, 0]; b = a + [0, 0, 1]; d = c + [0, 0, 1]
    pz, mz = mode(4, 0, kpl), mode(5, 0, kpl)
    t = lambda v: tuple(int(x) for x in v)
    return dict(remove=[((t(a), pz), (t(b), mz)), ((t(c), pz), (t(d), mz))],
                add=[((t(a), pz), (t(c), pz)), ((t(b), mz), (t(d), mz))],
                modes=[(t(a), pz), (t(b), mz), (t(c), pz), (t(d), mz)])


def block_origin_for(a0, l, where):
    """Block of side l whose top face (z max) contains a = a0 and c = a0 + x at the given position."""
    a = np.array(a0)
    zlo = a[2] - (l - 1)
    if where == "face":
        xlo = a[0] - (l - 2) // 2; ylo = a[1] - (l - 1) // 2
    elif where == "edge":                    # a, c on the y = max edge of the top face
        xlo = a[0] - (l - 2) // 2; ylo = a[1] - (l - 1)
    elif where == "corner":                  # c at the (x max, y max) corner
        xlo = a[0] + 1 - (l - 1); ylo = a[1] - (l - 1)
    else:
        raise ValueError(where)
    return np.array([xlo, ylo, zlo])


def woodbury_block_P(L, lam, spec, coords, kpl=1, nq=160, scale=1.0):
    """Separate path (G4): P = 1/2 + (1/pi) int_0^inf Herm G(iy) dy, with the rewired G from the translation-invariant
    G0 by a rank-4 Woodbury update. Returns P on the block (rewired)."""
    M = nmodes(kpl)
    H, _ = bloch_h(L, lam, kpl)
    e = np.linalg.eigvalsh(H)
    mu = 0.5 * (e[..., M // 2].min() + e[..., M // 2 - 1].max())     # Fermi level in the gap
    R = spec["modes"]
    # D: change of h on the four modes, ordering R = [a+, b-, c+, d-]
    D = np.zeros((4, 4))
    D[0, 1] = D[1, 0] = 1.0; D[2, 3] = D[3, 2] = 1.0       # remove -1 hoppings
    D[0, 2] = D[2, 0] = -1.0; D[1, 3] = D[3, 1] = -1.0     # add -1 hoppings
    t, w = np.polynomial.legendre.leggauss(nq)
    t = 0.5 * (t + 1); w = 0.5 * w
    ys = scale * t / (1 - t); wy = w * scale / (1 - t) ** 2
    n = len(coords)
    Rx = np.array([r[0] for r in R]); Rm = np.array([r[1] for r in R])
    dP = np.zeros((n * M, n * M), complex)
    I = np.eye(M)
    for y, wt in zip(ys, wy):
        Gk = np.linalg.inv((mu + 1j * y) * I - H)
        Gr = np.fft.ifftn(Gk, axes=(0, 1, 2))
        dBR = (coords[:, None, :] - Rx[None, :, :]) % L         # (n,4,3)
        GBR = np.stack([Gr[dBR[:, j, 0], dBR[:, j, 1], dBR[:, j, 2]][:, :, Rm[j]] for j in range(4)], 1)  # (n,4,M)
        GBR = GBR.transpose(0, 2, 1).reshape(n * M, 4)
        dRB = (Rx[:, None, :] - coords[None, :, :]) % L
        GRB = np.stack([Gr[dRB[j, :, 0], dRB[j, :, 1], dRB[j, :, 2]][:, Rm[j], :] for j in range(4)], 0)  # (4,n,M)
        GRB = GRB.reshape(4, n * M)
        dRR = (Rx[:, None, :] - Rx[None, :, :]) % L
        GRR = np.array([[Gr[tuple(dRR[i, j])][Rm[i], Rm[j]] for j in range(4)] for i in range(4)])
        dG = GBR @ D @ np.linalg.solve(np.eye(4) - GRR @ D, GRB)
        dP += wt * 0.5 * (dG + dG.conj().T) / np.pi
    return dP


# ---------------------------------------------------------------- matter deformation
def matter_pairs(L, lam, kpl, c0, k0_index=2, qz=1):
    """Pairs (v,k) -> (c,k+q), |k| <= 2 pi k0_index / L, with q = 2 pi qz / L z. The relative phase of each pair is
    fixed so that its first-order energy density at the centre site c0 is positive and maximal."""
    e, u, K = bands(L, lam, kpl)
    M = e.shape[-1]; v, c = M // 2 - 1, M // 2
    n = np.fft.fftfreq(L) * L
    pairs = []
    for i in range(L):
        for j in range(L):
            for k in range(L):
                if n[i] ** 2 + n[j] ** 2 + n[k] ** 2 > k0_index ** 2 + 1e-9:
                    continue
                kv = K[i, j, k]; k2 = (k + qz) % L
                kc = K[i, j, k2]
                uv = u[i, j, k, :, v]; uc = u[i, j, k2, :, c]
                ev, ec = e[i, j, k, v], e[i, j, k2, c]
                # first-order delta rho at c0 for phase phi: eps (ec+ev) Im(e^{i phi} psi_c psi_v^*) summed over site modes
                z = np.sum(uc * uv.conj()) * np.exp(1j * (kc - kv) @ np.asarray(c0, float)) / L ** 3
                amp = (ec + ev) * z
                phi = float(np.angle(1j * np.conj(amp))) if abs(amp) > 0 else 0.0   # makes Im(e^{i phi} amp) = |amp|
                pairs.append(dict(kv=kv, kc=kc, uv=uv, uc=uc * np.exp(1j * phi), ev=ev, ec=ec))
    return pairs


def plane(kvec, uvec, coords, L):
    ph = np.exp(1j * coords @ kvec) / np.sqrt(L ** 3)
    return (ph[:, None] * uvec[None, :]).reshape(-1)


def deform_block(PB, pairs, coords, L, eps):
    """Exact rotation of each pair: psi_v -> cos(eps) psi_v - i sin(eps) psi_c."""
    out = PB.copy()
    for p in pairs:
        pv = plane(p["kv"], p["uv"], coords, L); pc = plane(p["kc"], p["uc"], coords, L)
        w = np.cos(eps) * pv - 1j * np.sin(eps) * pc
        out += np.outer(w, w.conj()) - np.outer(pv, pv.conj())
    return out


def drho_sites(pairs, coords, L, M):
    """First-order energy density per unit eps at each site: (ec + ev) Im sum_m psi_c psi_v^*."""
    out = np.zeros(len(coords))
    for p in pairs:
        pv = plane(p["kv"], p["uv"], coords, L).reshape(len(coords), M)
        pc = plane(p["kc"], p["uc"], coords, L).reshape(len(coords), M)
        out += (p["ec"] + p["ev"]) * np.imag(np.sum(pc * pv.conj(), 1))
    return out


# ---------------------------------------------------------------- parts
def part_area(L=16, ells=(2, 3, 4, 5, 6)):
    res = {}
    for kpl, lams in ((1, LAMS), (2, (0.0,))):
        for lam in lams:
            e, u, _ = bands(L, lam, kpl)
            Pr = P_real(e, u)
            rows = {}
            for l in ells:
                if kpl == 2 and l > 4:
                    continue
                S = entropy(block_P(Pr, block_coords(l, (0, 0, 0)), L))
                rows[l] = dict(S=S, cut=6 * l * l, s_link=S / (6 * l * l))
            ls = np.array(sorted(rows)); Ss = np.array([rows[l]["S"] for l in ls])
            fit = np.polyfit(ls, Ss, 2) if len(ls) >= 3 else None
            res[f"kpl{kpl}_lam{lam}"] = dict(kpl=kpl, lam=lam, gap=gap(e), rows=rows,
                                            fit_l2_coeff_per_link=(float(fit[0]) / 6 if fit is not None else None),
                                            ln_d=float(np.log(2 ** kpl)))
    return res


def part_g1(L=6):
    out = {}
    for lam in (0.0, 0.3):
        e, u, _ = bands(L, lam)
        Pr = P_real(e, u)
        h, idx = real_h(L, lam)
        P, g = dense_P(h)
        coords = block_coords(3, (1, 2, 0))
        A = block_P(Pr, coords, L); B = dense_block(P, coords, idx, 6)
        out[str(lam)] = dict(max_dev=float(np.abs(A - B).max()), gap_dense=g, gap_bloch=gap(e))
        if lam == 0.0:
            ev = np.linalg.eigvalsh(A)
            out[str(lam)]["eig_half_count"] = int(np.sum(np.abs(ev - 0.5) < 1e-10))
            out[str(lam)]["cut_links"] = 6 * 9
    return out


def part_rewire(L=10, ells=(3, 4, 5), check_g4=True):
    res = {}
    a0 = (4, 4, 6)
    spec = rewire_spec(a0)
    for lam in LAMS:
        t0 = time.time()
        e, u, _ = bands(L, lam)
        Pr = P_real(e, u)
        h, idx = real_h(L, lam, rewire=spec)
        P, g = dense_P(h)
        rows = {}
        for l in ells:
            for where in ("face", "edge", "corner"):
                o = block_origin_for(a0, l, where)
                coords = block_coords(l, o) % L
                S0 = entropy(block_P(Pr, coords, L))
                S1 = entropy(dense_block(P, coords, idx, 6))
                r = dict(S0=S0, S1=S1, dS=S1 - S0, eta_count=-(S1 - S0) / 2, s_link=S0 / (6 * l * l))
                r["ratio"] = r["eta_count"] / r["s_link"]
                if check_g4 and l == 3:
                    dPB = woodbury_block_P(L, lam, spec, coords)
                    S1w = entropy(block_P(Pr, coords, L) + dPB)
                    r["G4_S1_woodbury"] = S1w; r["G4_dev"] = abs(S1w - S1)
                rows[f"{l}_{where}"] = r
        res[str(lam)] = dict(gap_rewired=g, gap_bloch=gap(e), rows=rows, wall_s=time.time() - t0)
        print(lam, json.dumps({k: round(v["ratio"], 4) for k, v in rows.items()}), flush=True)
    return res


def part_matter(L=16, ells=(2, 3, 4, 5, 6), eps_list=(1e-3, 1e-4)):
    res = {}
    c0 = np.array([8, 8, 8])
    for lam in LAMS:
        e, u, _ = bands(L, lam)
        Pr = P_real(e, u)
        pairs = matter_pairs(L, lam, 1, c0)
        rows = {}
        for l in ells:
            coords = block_coords(l, c0 - (l - 1) // 2)
            PB = block_P(Pr, coords, L)
            S0 = entropy(PB)
            drho = drho_sites(pairs, coords, L, 6)
            r = dict(S0=S0, drho_block_per_eps=float(drho.mean()), drho_centre_per_eps=float(drho[np.all(coords == c0, 1)][0]))
            for eps in eps_list:
                Sp = entropy(deform_block(PB, pairs, coords, L, eps)); Sm = entropy(deform_block(PB, pairs, coords, L, -eps))
                r[f"dS_plus_{eps:g}"] = Sp - S0; r[f"dS_minus_{eps:g}"] = Sm - S0
                r[f"dS_central_per_eps_{eps:g}"] = (Sp - Sm) / (2 * eps)
            eps = eps_list[0]
            r["G2_lin"] = abs(abs(r[f"dS_plus_{eps:g}"]) - abs(r[f"dS_minus_{eps:g}"])) / max(abs(r[f"dS_plus_{eps:g}"]), 1e-300)
            s_link = S0 / (6 * l * l)
            r["s_link"] = s_link
            if abs(r["drho_block_per_eps"]) > 1e-14:
                kB = r[f"dS_central_per_eps_{eps:g}"] / r["drho_block_per_eps"]
                r["k_B"] = kB
                r["G_eq"] = 3 * kB / (4 * np.pi * s_link * l ** 4)
            rows[l] = r
        res[str(lam)] = dict(gap=gap(e), npairs=len(pairs), rows=rows)
        print(lam, {l: (round(v.get("k_B", float("nan")), 4), v.get("G_eq")) for l, v in rows.items()}, flush=True)
    return res


def part_g3(L=16, ells=(2, 3), lams=(0.1, 0.3), dps=60):
    """First law in arbitrary precision: dS = Tr(K_B dP_B), K_B = ln((1-P_B)/P_B), with P_B built from mpmath Bloch
    projectors (P(k) = (1 - h/|h|)/2 via mp.eighe on each 6x6 h(k))."""
    import mpmath as mp
    mp.mp.dps = dps
    out = {}
    c0 = np.array([8, 8, 8])
    for lam in lams:
        H, _ = bloch_h(L, lam)
        Pk = {}
        for i in range(L):
            for j in range(L):
                for k in range(L):
                    A = mp.matrix(6, 6)
                    for a in range(6):
                        for b in range(6):
                            A[a, b] = mp.mpc(H[i, j, k, a, b].real, H[i, j, k, a, b].imag)
                    # exact phases at high precision
                    kk = [2 * mp.pi * mp.mpf(((x + L // 2) % L) - L // 2) / L for x in (i, j, k)]
                    for ax in range(3):
                        ph = mp.expj(kk[ax])
                        A[2 * ax, 2 * ax + 1] = -ph + mp.mpf(lam); A[2 * ax + 1, 2 * ax] = -mp.conj(ph) + mp.mpf(lam)   # link + on-site term
                    E, Q = mp.eighe(A)
                    occ = sorted(range(6), key=lambda n: E[n])[:3]     # mp.eighe does not sort its eigenvalues
                    P = mp.matrix(6, 6)
                    for a in range(6):
                        for b in range(6):
                            P[a, b] = sum(Q[a, n] * mp.conj(Q[b, n]) for n in occ)
                    Pk[(i, j, k)] = (P, kk)
        pairs = matter_pairs(L, lam, 1, c0)
        res = {}
        for l in ells:
            coords = block_coords(l, c0 - (l - 1) // 2)
            n = len(coords) * 6
            cache = {}
            def Pr(d):
                key = tuple(int(x) % L for x in d)
                if key not in cache:
                    acc = mp.matrix(6, 6)
                    for (P, kk) in Pk.values():
                        ph = mp.expj(sum(kk[a] * key[a] for a in range(3)))
                        acc += P * ph
                    cache[key] = acc / L ** 3
                return cache[key]
            PB = mp.matrix(n, n)
            for x in range(len(coords)):
                for y in range(len(coords)):
                    blk = Pr(coords[x] - coords[y])
                    for a in range(6):
                        for b in range(6):
                            PB[6 * x + a, 6 * y + b] = blk[a, b]
            E, Q = mp.eighe(PB)
            E = [mp.re(p) for p in E]
            tol = mp.mpf(10) ** (-dps // 2)
            pure = [bool(p < tol or 1 - p < tol) for p in E]        # exactly pure modes (structural zeros/ones)
            Kd = [mp.mpf(0) if pu else mp.log((1 - p) / p) for p, pu in zip(E, pure)]
            # first-order dP_B per eps (double precision is enough for dP): -i (psi_c psi_v^+ - psi_v psi_c^+)
            dP = np.zeros((n, n), complex)
            for p in pairs:
                pv = plane(p["kv"], p["uv"], coords, L); pc = plane(p["kc"], p["uc"], coords, L)
                dP += -1j * (np.outer(pc, pv.conj()) - np.outer(pv, pc.conj()))
            # Tr(K dP) = sum_n K_n <q_n| dP |q_n>
            Qn = np.array([[complex(Q[a, b]) for b in range(n)] for a in range(n)])
            diag = np.real(np.einsum("an,ab,bn->n", Qn.conj(), dP, Qn))
            dK = float(sum(float(Kd[m]) * diag[m] for m in range(n)))
            imp = [m for m in range(n) if not pure[m]]
            res[l] = dict(dK_per_eps=dK, n_pure=int(sum(pure)), max_abs_dP_on_pure=float(max([abs(diag[m]) for m in range(n) if pure[m]], default=0.0)),
                          min_eig_mixed=float(min(E[m] for m in imp)), max_eig_mixed=float(max(E[m] for m in imp)))
        out[str(lam)] = res
        print(lam, res, flush=True)
    return out


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.process_time()
    res = dict(area=part_area, g1=part_g1, rewire=part_rewire, matter=part_matter, g3=part_g3)[part]()
    res = dict(result=res, cpu_seconds=time.process_time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("done", part, res["cpu_seconds"])
