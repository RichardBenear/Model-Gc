"""EX-GC34M (exploratory): energy, momentum and stress sources of a scalar packet from one discrete matter action
varied with respect to J1's lapse, shift and edge variables; test J1's lattice conservation identities.
Plan: spec/EX-GC34M_plan.md.

  python src/exgc34m.py m0   -> results/EX-GC34M/m0.json          (derivation check)
  python src/exgc34m.py run  -> results/EX-GC34M/identities.json  (M1, M2, M3)
Reuses exgc34d (Kuhn lattice: ET, TETS, tet_edges) unchanged.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34d as D

OUT = ROOT / 'results/EX-GC34M'
ET = D.ET.astype(int)
DHAT = ET / np.linalg.norm(ET, axis=1)[:, None]
TETS = [np.array(o, int) for o in D.TETS]
EDGES = [D.tet_edges(o) for o in TETS]
MASS = 1.0
DUMP = dict(indent=1, default=lambda o: o.tolist() if isinstance(o, np.ndarray) else float(o) if isinstance(o, np.floating) else str(o))


def sh(a, o):
    """a(v + o) on the periodic grid (o integer 3-vector)."""
    return np.roll(a, tuple(-int(c) for c in o), axis=(0, 1, 2))


def put(a, o):
    """array b with b(v + o) = a(v)."""
    return np.roll(a, tuple(int(c) for c in o), axis=(0, 1, 2))


# ---------------- tetrahedron geometry from edge lengths ----------------
def gram(lmat):
    """lmat (...,4,4) symmetric edge lengths -> Gram (...,3,3) of edge vectors from vertex 0."""
    l2 = lmat ** 2
    G = np.empty(lmat.shape[:-2] + (3, 3))
    for a in range(3):
        for b in range(3):
            G[..., a, b] = 0.5 * (l2[..., 0, a + 1] + l2[..., 0, b + 1] - l2[..., a + 1, b + 1])
    return G


def flat_lmat(o, h):
    P = o.astype(float) * h
    return np.linalg.norm(P[:, None, :] - P[None, :, :], axis=-1)


def tet_index(o):
    return [k for k, x in enumerate(TETS) if (x == o).all()][0]


def geom(h):
    out = []
    for k, o in enumerate(TETS):
        L = flat_lmat(o, h); G = gram(L); Gi = np.linalg.inv(G); V = np.sqrt(np.linalg.det(G)) / 6
        Mb = np.zeros((3, 4)); Mb[:, 0] = -1; Mb[:, 1:] = np.eye(3)
        K = Mb.T @ (V * Gi) @ Mb
        w, dA, dV = [], [], []
        for (i, j, t) in EDGES[k]:
            w.append(-K[i, j])
            dL = np.zeros((4, 4)); dL[i, j] = dL[j, i] = 1.0
            dl2 = 2 * L * dL
            dG = np.array([[0.5 * (dl2[0, a + 1] + dl2[0, b + 1] - dl2[a + 1, b + 1]) for b in range(3)] for a in range(3)])
            dVe = 0.5 * V * np.trace(Gi @ dG)
            dA.append(dVe * Gi - V * Gi @ dG @ Gi); dV.append(dVe)
        out.append(dict(o=o, w=np.array(w), dA=np.array(dA), dV=np.array(dV), V=V))
    return out


# ---------------- analytic sources s = dS/dx (per unit time) ----------------
def sources(phi, phid, h, GE):
    """s_n (N^3), s_b (3,N^3), s_e (7,N^3) from the action, flat background, fields phi, phid."""
    m2 = MASS ** 2
    Vv = h ** 3
    s_n = -Vv * (0.5 * phid ** 2 + 0.5 * m2 * phi ** 2)
    s_b = np.zeros((3,) + phi.shape); s_e = np.zeros((7,) + phi.shape)
    Kd = 0.5 * phid ** 2 - 0.5 * m2 * phi ** 2           # volume coefficient
    for g in GE:
        o = g['o']
        ph = [sh(phi, o[a]) for a in range(4)]; pd = [sh(phid, o[a]) for a in range(4)]
        dphi = [ph[a + 1] - ph[0] for a in range(3)]
        Ksum = sum(sh(Kd, o[a]) for a in range(4))
        for e, (i, j, t) in enumerate(EDGES[tet_index(o)]):
            dA = g['dA'][e]
            q = sum(dA[a, b] * dphi[a] * dphi[b] for a in range(3) for b in range(3))
            s_e[t] += put(-0.5 * q + 0.25 * g['dV'][e] * Ksum, o[i])
            w = g['w'][e]
            if abs(w) > 1e-14 * h:
                assert t < 3, 'non-axis edge with weight'
                dph = ph[j] - ph[i]; pbar = 0.5 * (pd[i] + pd[j])
                s_n += put(-0.25 * w * dph ** 2, o[i]) + put(-0.25 * w * dph ** 2, o[j])
                s_b[t] += put(-w * h * pbar * dph, o[i])
    return s_n, s_b, s_e


def lap(phi, h):
    return sum(np.roll(phi, 1, a) + np.roll(phi, -1, a) - 2 * phi for a in range(3)) / h ** 2


def sources_dot(phi, phid, h, GE):
    """Exact time derivative by polarisation (sources quadratic in (phi, phid)); phidd from the field equation."""
    phidd = -MASS ** 2 * phi + lap(phi, h)
    lam = 1.0
    sp = sources(phi + lam * phid, phid + lam * phidd, h, GE)
    sm = sources(phi - lam * phid, phid - lam * phidd, h, GE)
    return [(a - b) / (2 * lam) for a, b in zip(sp, sm)]


def identities(s, sd, h):
    s_n, s_b, s_e = s; sn_d, sb_d, se_d = sd
    IE = -sn_d + sum(s_b[a] - np.roll(s_b[a], 1, a) for a in range(3)) / h
    div = np.zeros_like(s_b)
    for t in range(7):
        dd = sh(s_e[t], -ET[t]) - s_e[t]
        for a in range(3):
            div[a] += DHAT[t, a] * dd
    IP = -sb_d + div
    return IE, IP, div


# ---------------- direct Lagrangian (separate path for M0) ----------------
def lagrangian(phi, phid, n, b, e, h):
    """L(x) evaluated directly from perturbed tet geometry (lengths l0 + e), vertex volumes, lapse, shift."""
    m2 = MASS ** 2
    Vv = np.zeros(phi.shape); Egrad = 0.0; Lns = 0.0
    for k, o in enumerate(TETS):
        L0 = flat_lmat(o, h)
        lmat = np.broadcast_to(L0, phi.shape + (4, 4)).copy()
        for (i, j, t) in EDGES[k]:
            lij = L0[i, j] + sh(e[t], o[i])
            lmat[..., i, j] = lmat[..., j, i] = lij
        G = gram(lmat); V = np.sqrt(np.linalg.det(G)) / 6; Gi = np.linalg.inv(G)
        dphi = np.stack([sh(phi, o[a + 1]) - sh(phi, o[0]) for a in range(3)], -1)
        Egrad += 0.5 * np.sum(np.einsum('...a,...ab,...b->...', dphi, V[..., None, None] * Gi, dphi))
        for a in range(4):
            Vv += put(V / 4, o[a])
        # flat weights for lapse and shift couplings
        Gf = gram(L0[None])[0]; Vf = np.sqrt(np.linalg.det(Gf)) / 6
        Mb = np.zeros((3, 4)); Mb[:, 0] = -1; Mb[:, 1:] = np.eye(3); Kf = Mb.T @ (Vf * np.linalg.inv(Gf)) @ Mb
        for (i, j, t) in EDGES[k]:
            w = -Kf[i, j]
            if abs(w) > 1e-14 * h:
                pi_, pj = sh(phi, o[i]), sh(phi, o[j]); dph = pj - pi_
                nbar = 0.5 * (sh(n, o[i]) + sh(n, o[j])); pbar = 0.5 * (sh(phid, o[i]) + sh(phid, o[j]))
                be = h * sh(b[t], o[i])
                Lns += np.sum(w * (0.5 * nbar * dph ** 2 + pbar * be * dph))
    Lkin = np.sum(Vv * (0.5 * (1 - n) * phid ** 2 - 0.5 * (1 + n) * m2 * phi ** 2))
    return Lkin - Egrad - Lns


def m0():
    t0 = time.process_time(); rng = np.random.default_rng(3434 + 13); res = {}
    for h in (1.0, 0.5):
        N = 8; GE = geom(h)
        phi = rng.standard_normal((N,) * 3); phid = rng.standard_normal((N,) * 3)
        n = rng.standard_normal((N,) * 3); b = rng.standard_normal((3,) + (N,) * 3); e = rng.standard_normal((7,) + (N,) * 3)
        s_n, s_b, s_e = sources(phi, phid, h, GE)
        rows = []
        for eps in (1e-4, 1e-5):
            Lp = lagrangian(phi, phid, eps * n, eps * b, eps * h * e, h)
            Lm = lagrangian(phi, phid, -eps * n, -eps * b, -eps * h * e, h)
            pred_e = np.sum(s_n * n) + np.sum(s_b * b) + np.sum(s_e * h * e)
            fd = (Lp - Lm) / (2 * eps)
            rows.append(dict(eps=eps, fd=float(fd), pred=float(pred_e), rel=float(abs(fd - pred_e) / abs(pred_e))))
        # separate parts: lapse only, shift only, edges only
        parts = {}
        for nm, (nn, bb, ee) in dict(lapse=(n, 0 * b, 0 * e), shift=(0 * n, b, 0 * e), edges=(0 * n, 0 * b, h * e)).items():
            eps = 1e-5
            fd = (lagrangian(phi, phid, eps * nn, eps * bb, eps * ee, h) - lagrangian(phi, phid, -eps * nn, -eps * bb, -eps * ee, h)) / (2 * eps)
            pr = np.sum(s_n * nn) + np.sum(s_b * bb) + np.sum(s_e * ee)
            parts[nm] = float(abs(fd - pr) / abs(pr))
        # Richardson (removes the O(eps^2) truncation of the central difference)
        def cd(eps):
            return (lagrangian(phi, phid, eps * n, eps * b, eps * h * e, h) - lagrangian(phi, phid, -eps * n, -eps * b, -eps * h * e, h)) / (2 * eps)
        rich = (4 * cd(1e-4) - cd(2e-4)) / 3
        rel_rich = float(abs(rich - pred_e) / abs(pred_e))
        res['h%g' % h] = dict(rows=rows, parts=parts, richardson_rel=rel_rich,
                              M0=bool(rel_rich <= 1e-8 and max(parts.values()) <= 1e-8),
                              note='first code used plain central differences (O(eps^2) truncation 1.2e-8 at eps 1e-5); Richardson added')
        print('richardson', rel_rich)
        print(h, rows, parts, flush=True)
    res['M0'] = all(v['M0'] for k, v in res.items() if k.startswith('h'))
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'm0.json').write_text(json.dumps(res, **DUMP))
    print('M0', res['M0'])


# ---------------- packet ----------------
LBOX = 32.0; K0 = 0.6; SIG = 3.0
DIRS = dict(axis=(1.0, 0, 0), generic=(0.8, 0.5, 0.33))


def packet(N, h, dname, t):
    kk = 2 * np.pi * np.fft.fftfreq(N, d=h)
    KX, KY, KZ = np.meshgrid(kk, kk, kk, indexing='ij')
    d = np.array(DIRS[dname]); d = d / np.linalg.norm(d); k0 = K0 * d
    x0 = np.full(3, LBOX / 2) - 6 * d             # start behind centre
    A = np.exp(-((KX - k0[0]) ** 2 + (KY - k0[1]) ** 2 + (KZ - k0[2]) ** 2) * SIG ** 2 / 2) * np.exp(-1j * (KX * x0[0] + KY * x0[1] + KZ * x0[2]))
    w = np.sqrt(MASS ** 2 + (4 * np.sin(KX * h / 2) ** 2 + 4 * np.sin(KY * h / 2) ** 2 + 4 * np.sin(KZ * h / 2) ** 2) / h ** 2)
    Z = np.fft.ifftn(A * np.exp(-1j * w * t)) * N ** 3 * (h / LBOX) ** 3 * (SIG ** 3)   # amplitude independent of N
    Zd = np.fft.ifftn(-1j * w * A * np.exp(-1j * w * t)) * N ** 3 * (h / LBOX) ** 3 * (SIG ** 3)
    return Z.real, Zd.real, (KX, KY, KZ)


def continuum(phi, phid, h, KK):
    """Integrated energy, momentum, stress of the band-limited interpolant (spectral)."""
    N = phi.shape[0]; c = h ** 3 / N ** 3
    P = np.fft.fftn(phi); Pd = np.fft.fftn(phid)
    grads = [np.fft.ifftn(1j * K * P).real for K in KK]
    Ekin = 0.5 * np.sum(phid ** 2) * h ** 3; Em = 0.5 * MASS ** 2 * np.sum(phi ** 2) * h ** 3
    Egr = 0.5 * sum(np.sum(g ** 2) for g in grads) * h ** 3
    mom = np.array([-np.sum(phid * g) * h ** 3 for g in grads])
    T = np.array([[np.sum(grads[i] * grads[j]) * h ** 3 for j in range(3)] for i in range(3)])
    T -= np.eye(3) * (Egr + Em - Ekin)
    return Ekin + Em + Egr, mom, T


def run():
    t0 = time.process_time(); res = dict(cases=[])
    for dname in DIRS:
        for h in (1.0, 0.5, 0.25):
            N = int(round(LBOX / h)); GE = geom(h)
            for t in (0.0, 6.0, 12.0):
                phi, phid, KK = packet(N, h, dname, t)
                s = sources(phi, phid, h, GE); sd = sources_dot(phi, phid, h, GE)
                IE, IP, div = identities(s, sd, h)
                rE = float(np.abs(IE).max() / np.abs(sd[0]).max())
                rP = float(np.abs(IP).max() / max(np.abs(sd[1]).max(), np.abs(div).max()))
                # M3
                Elat = -float(np.sum(s[0])); Plat = s[1].reshape(3, -1).sum(1)
                Tlat = sum(s[2][t_].sum() * np.linalg.norm(ET[t_]) * h * np.outer(DHAT[t_], DHAT[t_]) for t_ in range(7))
                Ec, Pc, Tc = continuum(phi, phid, h, KK)
                row = dict(dir=dname, h=h, N=N, t=t, rE=rE, rP=rP, maxIP=float(np.abs(IP).max()),
                           scale_sbdot=float(np.abs(sd[1]).max()), scale_div=float(np.abs(div).max()),
                           E_lat=Elat, E_cont=Ec, dE=abs(Elat - Ec) / Ec,
                           P_lat=Plat, P_cont=Pc, dP=float(np.linalg.norm(Plat - Pc) / np.linalg.norm(Pc)),
                           T_lat=Tlat, T_cont=Tc, dT=float(np.abs(Tlat - Tc).max() / np.abs(Tc).max()))
                res['cases'].append(row)
                print(dname, h, t, 'rE %.2e rP %.3e dE %.2e dP %.2e dT %.2e cpu %.0f' % (rE, rP, row['dE'], row['dP'], row['dT'], time.process_time() - t0), flush=True)
    C = res['cases']
    res['M1'] = bool(max(c['rE'] for c in C) <= 1e-10)
    exact = max(c['rP'] for c in C) <= 1e-10
    ratios = []
    for dname in DIRS:
        for t in (0.0, 6.0, 12.0):
            r = [next(c['rP'] for c in C if c['dir'] == dname and c['h'] == h and c['t'] == t) for h in (1.0, 0.5, 0.25)]
            ratios.append(dict(dir=dname, t=t, r1=r[0] / r[1], r2=r[1] / r[2]))
    res['M2_ratios'] = ratios
    res['M2'] = 'M2-EXACT' if exact else ('M2-SECOND' if all(min(x['r1'], x['r2']) >= 3.5 for x in ratios) else 'M2-FAIL')
    res['cpu_s'] = time.process_time() - t0
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'identities.json').write_text(json.dumps(res, **DUMP))
    print('M1', res['M1'], 'M2', res['M2'], [(x['dir'], x['t'], round(x['r1'], 2), round(x['r2'], 2)) for x in ratios], 'cpu', res['cpu_s'])


def posthoc():
    """POST HOC (labelled): (i) global momentum conservation; (ii) momentum density re-centred from the half-cell
    (forward-edge) placement to the vertex: momentum identity vs energy identity trade-off."""
    t0 = time.process_time(); rows = []
    for dname in DIRS:
        for h in (1.0, 0.5, 0.25):
            N = int(round(LBOX / h)); GE = geom(h); phi, phid, KK = packet(N, h, dname, 6.0)
            s = sources(phi, phid, h, GE); sd = sources_dot(phi, phid, h, GE)
            IE, IP, div = identities(s, sd, h)
            sc = max(np.abs(sd[1]).max(), np.abs(div).max())
            glob = float(np.abs(IP.reshape(3, -1).sum(1)).max() / np.abs(sd[1]).max())
            sb2 = np.array([0.5 * (s[1][a] + np.roll(s[1][a], 1, a)) for a in range(3)])
            sbd2 = np.array([0.5 * (sd[1][a] + np.roll(sd[1][a], 1, a)) for a in range(3)])
            IP2 = -sbd2 + div
            IE2 = -sd[0] + sum(sb2[a] - np.roll(sb2[a], 1, a) for a in range(3)) / h
            rows.append(dict(dir=dname, h=h, t=6.0, global_rel=glob, rP=float(np.abs(IP).max() / sc),
                             rP_recentred=float(np.abs(IP2).max() / sc), rE=float(np.abs(IE).max() / np.abs(sd[0]).max()),
                             rE_recentred=float(np.abs(IE2).max() / np.abs(sd[0]).max())))
            print(rows[-1], flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'posthoc.json').write_text(json.dumps(dict(note='POST HOC', rows=rows, cpu_s=time.process_time() - t0), **DUMP))


if __name__ == '__main__':
    dict(m0=m0, run=run, posthoc=posthoc)[sys.argv[1]]()
