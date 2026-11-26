"""EX-GC40C: is a shared-entanglement hub on the 40A pair network a mass, or only a local defect?

Hub: lam0 + dlam on one site (or a 7-site cluster). Everything follows from the background Green's function by
Woodbury updates on the imaginary axis z = mu + iy:
  ground-state energy shift  dE0 = -(1/pi) int_0^inf ln|det(I - G0_RR D)| dy      (fixed N, mu in the gap)
  projector change           dP  = (1/pi) int_0^inf Herm[G0_TR D (I - G0_RR D)^-1 G0_RT] dy
Bound states: zero eigenvalues of D^-1 - G0_RR(E) for real E in the gap.
Plan: spec/EX-GC40C_plan.md. Our construction; not LaViolette's Model G and not the Chain Theory construction.
"""
import json, sys, time
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc40a as A

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC40C"
M = 6
DLAMS = (0.1, 0.2, 0.4)
NEIGH = [(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]


def onsite(dl):
    return dl * (np.ones((M, M)) - np.eye(M))


def mu_gap(L, lam):
    e = np.linalg.eigvalsh(A.bloch_h(L, lam)[0])
    lo, hi = e[..., M // 2 - 1].max(), e[..., M // 2].min()
    return 0.5 * (lo + hi), float(hi - lo), float(lo), float(hi)


def Gr_at(H, z):
    I = np.eye(M)
    return np.fft.ifftn(np.linalg.inv(z * I - H), axes=(0, 1, 2))


def block(Gr, sitesA, sitesB, L):
    """G0 between mode lists: sites given as (n,3) int arrays, all 6 modes each -> (6nA, 6nB)."""
    d = (sitesA[:, None, :] - sitesB[None, :, :]) % L
    g = Gr[d[..., 0], d[..., 1], d[..., 2]]                 # (nA, nB, 6, 6)
    return g.transpose(0, 2, 1, 3).reshape(len(sitesA) * M, len(sitesB) * M)


def targets(L, rmax):
    """Probe sites along +x (r,0,0) and the face diagonal (r,r,0), r = 1..rmax."""
    ax = np.array([(r, 0, 0) for r in range(1, rmax + 1)])
    dg = np.array([(r, r, 0) for r in range(1, rmax + 1)])
    return ax, dg


def bond_energy_change(dP_site_link, a_end, b_end):
    """Change of the link energy -2 Re P_{a,b} for the link joining modes a (site s) and b (site s+1)."""
    return -2.0 * np.real(dP_site_link)


def site_entropy_first_order(P0ss, dPss):
    w, Q = np.linalg.eigh(0.5 * (P0ss + P0ss.conj().T))
    w = np.clip(w, 1e-15, 1 - 1e-15)
    K = Q @ np.diag(np.log((1 - w) / w)) @ Q.conj().T
    return float(np.real(np.trace(K @ dPss)))


def run_case(L, lam0, nq=160, rmax=None, scale=None, two_hub_d=range(2, 13), far=True):
    t0 = time.time()
    H, _ = A.bloch_h(L, lam0)
    mu, gap, lo, hi = mu_gap(L, lam0)
    scale = scale or max(gap, 0.2)
    rmax = rmax or L // 2 - 1
    t, w = np.polynomial.legendre.leggauss(nq)
    t = 0.5 * (t + 1); w = 0.5 * w
    ys = scale * t / (1 - t); wy = w * scale / (1 - t) ** 2
    hub1 = np.array([(0, 0, 0)]); hub7 = np.array(NEIGH)
    ax, dg = targets(L, rmax)
    # far-field probes: each probe site plus its +x neighbour (for the link energy along x)
    probes = np.concatenate([ax, dg]); probes_nx = probes + np.array([1, 0, 0])
    D1 = {dl: onsite(dl) for dl in DLAMS}
    D7 = {dl: np.kron(np.eye(7), onsite(dl)) for dl in DLAMS}
    dE = {f"1_{dl}": 0.0 for dl in DLAMS} | {f"7_{dl}": 0.0 for dl in DLAMS}
    dE2 = {d: 0.0 for d in two_hub_d}
    dl_main = 0.2
    acc_link = {k: np.zeros(len(probes)) for k in ("1", "7")}
    acc_site = {k: np.zeros((len(probes), M, M), complex) for k in ("1", "7")}
    for y, wt in zip(ys, wy):
        Gr = Gr_at(H, mu + 1j * y)
        G11 = block(Gr, hub1, hub1, L); G77 = block(Gr, hub7, hub7, L)
        for dl in DLAMS:
            dE[f"1_{dl}"] += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(M) - G11 @ D1[dl])))
            dE[f"7_{dl}"] += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(7 * M) - G77 @ D7[dl])))
        for d in two_hub_d:
            hubs = np.array([(0, 0, 0), (d, 0, 0)])
            G22 = block(Gr, hubs, hubs, L); D2 = np.kron(np.eye(2), D1[dl_main])
            dE2[d] += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(2 * M) - G22 @ D2)))
        if far:
            for key, hub, Dm, Ghh in (("1", hub1, D1[dl_main], G11), ("7", hub7, D7[dl_main], G77)):
                core = np.linalg.solve(np.eye(len(Dm)) - Ghh @ Dm, np.eye(len(Dm)))
                GPh = block(Gr, probes, hub, L); GhP = block(Gr, hub, probes, L)
                GXh = block(Gr, probes_nx, hub, L); GhX = block(Gr, hub, probes_nx, L)
                T = Dm @ core
                # site blocks of dG: (6n x 6n) only diagonal site blocks needed
                left = GPh @ T
                leftX = GXh @ T
                for p in range(len(probes)):
                    sl = slice(p * M, (p + 1) * M)
                    dGss = left[sl] @ GhP[:, sl]
                    acc_site[key][p] += wt / np.pi * 0.5 * (dGss + dGss.conj().T)
                    # link (p,+x) <-> (p+x,-x): modes 0 at p and 1 at p+x
                    a_row = left[p * M + 0]
                    dG_ab = a_row @ GhX[:, p * M + 1]
                    b_row = leftX[p * M + 1]
                    dG_ba = b_row @ GhP[:, p * M + 0]
                    # Herm part element (a,b): (dG_ab + conj(dG_ba))/2
                    acc_link[key][p] += wt / np.pi * np.real(0.5 * (dG_ab + np.conj(dG_ba)))
    out = dict(L=L, lam0=lam0, mu=float(mu), gap=gap, nq=nq, dE=dE,
               E_int={d: float(dE2[d] - 2 * dE[f"1_{dl_main}"]) for d in two_hub_d})
    if far:
        Pr = A.P_real(*A.bands(L, lam0)[:2])
        P0ss = Pr[0, 0, 0]
        prof = {}
        for key in ("1", "7"):
            link = -2.0 * acc_link[key]           # bond energy change: h_ab = -1, energy 2 Re(h_ab P_ba) = -2 Re P_ab
            site = [site_entropy_first_order(P0ss, acc_site[key][p]) for p in range(len(probes))]
            prof[key] = dict(axis=dict(r=list(range(1, rmax + 1)), d_link_energy=link[:rmax].tolist(), dS_site=site[:rmax]),
                             diag=dict(r=[float(np.sqrt(2) * r) for r in range(1, rmax + 1)], d_link_energy=link[rmax:].tolist(),
                                       dS_site=site[rmax:]))
        out["far_field_dlam0.2"] = prof
    out["bound_states"] = bound_states(H, L, mu, lo, hi)
    out["wall_s"] = time.time() - t0
    return out


def bound_states(H, L, mu, lo, hi, nE=400):
    """Count zero crossings of eigenvalues of D^-1 - G0_RR(E) for E in the gap (single-site hub)."""
    if hi - lo <= 1e-9:
        return None
    Es = np.linspace(lo + 1e-3 * (hi - lo), hi - 1e-3 * (hi - lo), nE)
    hub = np.array([(0, 0, 0)])
    Dinv = {dl: np.linalg.inv(onsite(dl)) for dl in DLAMS}
    prev = {dl: None for dl in DLAMS}; cross = {dl: [] for dl in DLAMS}
    for E in Es:
        G = block(Gr_at(H, E + 0j), hub, hub, L)
        for dl in DLAMS:
            Mx = Dinv[dl] - G
            ev = np.linalg.eigvalsh(0.5 * (Mx + Mx.conj().T))
            if prev[dl] is not None:
                for k in range(M):
                    if np.sign(ev[k]) != np.sign(prev[dl][k]):
                        cross[dl].append(float(E))
            prev[dl] = ev
    return {str(dl): dict(n=len(cross[dl]), energies=[round(x, 4) for x in cross[dl]], mu=float(mu)) for dl in DLAMS}


def fit_profile(r, f):
    """Fit |f| to A e^{-r/xi} / r^p (log-linear in r and ln r) and to B/r; return residuals in log space."""
    r = np.asarray(r, float); f = np.abs(np.asarray(f, float))
    m = f > 1e-300
    r, f = r[m], f[m]
    X = np.column_stack([np.ones_like(r), -r, -np.log(r)])
    c, *_ = np.linalg.lstsq(X, np.log(f), rcond=None)
    res_exp = float(np.sum((X @ c - np.log(f)) ** 2))
    b = np.mean(np.log(f) + np.log(r))
    res_r = float(np.sum((b - np.log(r) - np.log(f)) ** 2))
    return dict(xi=float(1 / c[1]) if c[1] != 0 else None, p=float(c[2]), res_exp=res_exp, res_1overr=res_r,
                ratio=res_r / max(res_exp, 1e-300))


def corr_length(L, lam0):
    Pr = A.P_real(*A.bands(L, lam0)[:2])
    rs = np.arange(1, L // 2)
    c = np.array([np.abs(Pr[r, 0, 0]).max() for r in rs])
    return fit_profile(rs[2:], c[2:]) | dict(C=c.tolist())


# ---------------------------------------------------------------- checks
def check_c1(L=8, lam0=0.1, dl=0.2):
    """Woodbury (contour) vs dense diagonalisation: dE0 and site-block dP at two probe sites, 1- and 7-site hubs."""
    out = {}
    h0, idx = A.real_h(L, lam0)
    e0 = np.linalg.eigvalsh(h0); n = h0.shape[0] // 2
    P0, _ = A.dense_P(h0)
    H, _ = A.bloch_h(L, lam0)
    mu, gap, lo, hi = mu_gap(L, lam0)
    t, w = np.polynomial.legendre.leggauss(160); t = 0.5 * (t + 1); w = 0.5 * w
    ys = gap * t / (1 - t); wy = w * gap / (1 - t) ** 2
    for key, hub in (("1", np.array([(0, 0, 0)])), ("7", np.array(NEIGH))):
        h = h0.copy()
        for s in hub:
            b = idx(tuple(s), 0); h[b:b + M, b:b + M] += onsite(dl)
        e = np.linalg.eigvalsh(h); P, _ = A.dense_P(h)
        dE_dense = float(e[:n].sum() - e0[:n].sum())
        Dm = np.kron(np.eye(len(hub)), onsite(dl))
        probes = np.array([(2, 0, 0), (3, 1, 0)])
        dE = 0.0; dP = np.zeros((2 * M, 2 * M), complex)
        for y, wt in zip(ys, wy):
            Gr = Gr_at(H, mu + 1j * y)
            Ghh = block(Gr, hub, hub, L)
            dE += -wt / np.pi * np.log(abs(np.linalg.det(np.eye(len(Dm)) - Ghh @ Dm)))
            dG = block(Gr, probes, hub, L) @ Dm @ np.linalg.solve(np.eye(len(Dm)) - Ghh @ Dm, block(Gr, hub, probes, L))
            dP += wt / np.pi * 0.5 * (dG + dG.conj().T)
        ids = [idx(tuple(p), m) for p in probes for m in range(M)]
        dP_dense = (P - P0)[np.ix_(ids, ids)]
        out[key] = dict(dE_dense=dE_dense, dE_woodbury=float(dE), dE_dev=abs(dE - dE_dense),
                        dP_dev=float(np.abs(dP - dP_dense).max()), dP_max=float(np.abs(dP_dense).max()))
    return out


def check_locality_dense(L=8, dl=0.2):
    """P1 at lam0 = 0 by dense diagonalisation: max |dP| on modes beyond distance 1 of the hub."""
    h0, idx = A.real_h(L, 0.0); P0, _ = A.dense_P(h0)
    h = h0.copy(); h[0:M, 0:M] += onsite(dl); P, _ = A.dense_P(h)
    dP = np.abs(P - P0)
    near = set([idx((0, 0, 0), m) for m in range(M)])
    for a, s in ((0, (1, 0, 0)), (1, (-1, 0, 0)), (2, (0, 1, 0)), (3, (0, -1, 0)), (4, (0, 0, 1)), (5, (0, 0, -1))):
        near.add(idx(s, a ^ 1))                 # the partner mode of hub mode a
    far = np.array([i for i in range(h.shape[0]) if i not in near])
    e0 = np.linalg.eigvalsh(h0); e = np.linalg.eigvalsh(h); n = h.shape[0] // 2
    return dict(max_dP_far=float(dP[np.ix_(far, far)].max()), max_dP_near_far=float(dP[np.ix_(far, sorted(near))].max()),
                max_dP_near=float(dP[np.ix_(sorted(near), sorted(near))].max()), dE0=float(e[:n].sum() - e0[:n].sum()))


def gapless_dense(L, dl=0.2, lam0=0.5):
    h0, idx = A.real_h(L, lam0)
    e0, v0 = np.linalg.eigh(h0); n = h0.shape[0] // 2
    P0 = v0[:, :n] @ v0[:, :n].T
    h = h0.copy(); h[0:M, 0:M] += onsite(dl)
    e, v = np.linalg.eigh(h); P = v[:, :n] @ v[:, :n].T
    rows = []
    for r in range(1, L // 2):
        a, b = idx((r, 0, 0), 0), idx((r + 1, 0, 0), 1)
        rows.append(dict(r=r, d_link_energy=float(-2 * (P - P0)[a, b])))
    return dict(L=L, dE0=float(e[:n].sum() - e0[:n].sum()), gap_at_half=float(e0[n] - e0[n - 1]),
                gap_at_half_hub=float(e[n] - e[n - 1]), axis=rows)


def local_excitation(L=24, lam0=0.1, rmax=10):
    """P6 control: move one electron from a localized valence packet to a localized conduction packet at the origin.
    dP = -|phi_v><phi_v| + |phi_c><phi_c|, phi_v = P0 e_m / |.|, phi_c = (1-P0) e_m / |.| (mode m = 0 of the origin)."""
    Pr = A.P_real(*A.bands(L, lam0)[:2])
    # column of P0 for mode 0 at the origin: P0_{(x,a),(0,0)} = Pr[x][a,0]
    colP = Pr[..., :, 0]                          # (L,L,L,6)
    e0 = np.zeros((L, L, L, M), complex); e0[0, 0, 0, 0] = 1
    pv = colP / np.sqrt(np.sum(np.abs(colP) ** 2)); pc = (e0 - colP); pc = pc / np.sqrt(np.sum(np.abs(pc) ** 2))
    rows = []
    for r in range(1, rmax + 1):
        a = pv[r, 0, 0, 0]; b = pv[(r + 1) % L, 0, 0, 1]; ca = pc[r, 0, 0, 0]; cb = pc[(r + 1) % L, 0, 0, 1]
        dPab = -a * np.conj(b) + ca * np.conj(cb)
        wpk = float(np.sum(np.abs(pv[r, 0, 0]) ** 2 + np.abs(pc[r, 0, 0]) ** 2))
        rows.append(dict(r=r, d_link_energy=float(-2 * np.real(dPab)), packet_weight=wpk))
    return rows


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if part == "checks":
        res = dict(C1=check_c1(), locality_dense=check_locality_dense())
        a = run_case(12, 0.1, nq=160, rmax=4, two_hub_d=(2, 4)); b = run_case(12, 0.1, nq=320, rmax=4, two_hub_d=(2, 4))
        res["C3"] = dict(dE_dev=max(abs(a["dE"][k] - b["dE"][k]) for k in a["dE"]),
                         far_dev=max(abs(x - y) for k in ("1", "7") for x, y in zip(a["far_field_dlam0.2"][k]["axis"]["d_link_energy"],
                                                                                   b["far_field_dlam0.2"][k]["axis"]["d_link_energy"])))
        res["C2_note"] = "dlam -> 0 limit: dE scales as dlam^2 (see main), exact zero at dlam = 0 by construction"
    elif part == "main":
        res = {}
        for lam0 in (0.0, 0.1, 0.2):
            res[str(lam0)] = {}
            for L in (16, 24, 32):
                res[str(lam0)][L] = run_case(L, lam0)
                print(lam0, L, res[str(lam0)][L]["dE"], res[str(lam0)][L]["wall_s"], flush=True)
            res[str(lam0)]["corr"] = corr_length(32, lam0) if lam0 > 0 else None
    elif part == "gapless":
        res = {L: gapless_dense(L) for L in (8, 10)}
    elif part == "control":
        res = dict(local_excitation=local_excitation())
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=float) + "\n")
    print("done", part, res["wall_seconds"])
