"""EX-GC34B (exploratory): sources of the shared structure (lapse N, spatial metric gamma_ij) for knots and packets.

  python src/exgc34b.py states      # regenerate B1/B2 (h=0.2) and B1 (h=0.15) rest knots with the 26E1 optimiser
  python src/exgc34b.py sources     # numerical gamma/N derivatives + closed-form recomputation -> results/EX-GC34B
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[k] = '1'
import sys, json, hashlib, itertools, platform, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
OUT = ROOT / 'results/EX-GC34B'; RAW = ROOT / 'local_raw/EX-GC34B'
EPS = 1e-5


def e1mod(h):
    import exgc26e1 as e1
    e1.configure(h)
    e1.OUT = OUT / f'states_h{h:g}'; e1.RAW = RAW / f'h{h:g}'
    return e1


def states():
    for h, B in ((0.2, 1), (0.2, 2), (0.15, 1)):
        e1 = e1mod(h); e1.load()
        if (e1.OUT / f'static_B{B}.json').exists():
            continue
        e1.static(B, 300)


# ---------- knot energy with a uniform spatial metric (26E1 orientation stencil) ----------
def knot_parts(G, p):
    """Per-orientation masked D tensors summed into the needed contractions."""
    import exgc26b2 as b2
    out = []
    for s in itertools.product((-1, 1), repeat=3):
        m, d, D = b2.orient(G, p, s)
        out.append((m, np.array([[D[i][j] for j in range(3)] for i in range(3)])))
    return out


def knot_E(parts, gam, fac):
    gi = np.linalg.inv(gam); sg = np.sqrt(np.linalg.det(gam))
    E2 = E4 = 0.0
    for m, D in parts:
        GD = np.einsum('ik,kj...->ij...', gi, D)          # g^{ik} D_kj
        tr = np.einsum('ii...->...', GD)
        GDGD = np.einsum('ij...,ji...->...', GD, GD)
        E2 += (tr * m).sum() / 8; E4 += (0.5 * (tr * tr - GDGD) * m).sum() / 8
    return sg * E2 * fac, sg * E4 * fac


def knot_T_numeric(parts, fac):
    T = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            dg = np.zeros((3, 3)); dg[i, j] += EPS; dg[j, i] += EPS
            if i == j:
                dg[i, i] = EPS
            Ep = sum(knot_E(parts, np.eye(3) + dg, fac)); Em = sum(knot_E(parts, np.eye(3) - dg, fac))
            dE = (Ep - Em) / (2 * EPS)
            T[i, j] = T[j, i] = -2 * dE if i == j else -dE      # off-diagonal: both gamma_ij and gamma_ji moved
    return T


def knot_T_closed(parts, fac):
    T = np.zeros((3, 3)); E2 = E4 = 0.0
    for m, D in parts:
        tr = np.einsum('ii...->...', D); D2 = np.einsum('ik...,kj...->ij...', D, D)
        e2 = tr; e4 = 0.5 * (tr * tr - np.einsum('ii...->...', D2))
        Tij = 2 * D + 2 * (tr[None, None] * D - D2) - np.eye(3)[:, :, None, None, None] * (e2 + e4)[None, None]
        T += (Tij * m).sum((2, 3, 4)) / 8; E2 += (e2 * m).sum() / 8; E4 += (e4 * m).sum() / 8
    return T * fac, E2 * fac, E4 * fac


def knot_case(h, B):
    e1 = e1mod(h); import exgc26a_skyrmion as ref
    G = ref.Grid(h, 7.0); st = json.loads((e1.OUT / f'static_B{B}.json').read_text())
    p = np.load(e1.restpath(B)); fac = G.h ** 3 / e1.NORM
    parts = knot_parts(G, p)
    E2, E4 = knot_E(parts, np.eye(3), fac); E = E2 + E4
    Tn = knot_T_numeric(parts, fac); Tc, E2c, E4c = knot_T_closed(parts, fac)
    # lapse: L = K/N - N U, K = 0 at rest -> s_N = -dL/dN = U = E
    sN = E
    ref26 = json.loads((ROOT / ('results/EX-GC26E1' if h == 0.2 else 'results/EX-GC26E2') / f'static_B{B}.json').read_text())
    return dict(h=h, B=B, E=E, E2=E2, E4=E4, sN=sN, T_numeric=Tn.tolist(), T_closed=Tc.tolist(),
                trace_numeric=float(np.trace(Tn)), E4_minus_E2=E4 - E2, offset_rel=(E4 - E2) / E,
                closed_vs_numeric_max=float(abs(Tc - Tn).max()), E2c=E2c, E4c=E4c,
                regen_E=st['E'], committed_E=ref26['E'], sha_regen=st['sha256'], sha_committed=ref26['sha256'],
                sha_match=st['sha256'] == ref26['sha256'])


# ---------- packets on the 9-link lattice ----------
def packet_case(mu, tt):
    import exgc34a as S
    G = S.Grid((96, 64, 64), (-9.6, -6.4, -6.4), 0.2); bg = S.Background(G, 'flat')
    x0 = np.array([-3.0, 0, 0]); k0 = np.array([3.0, 0, 0])
    i0 = tuple(int(round((x0[a] - G.lo[a]) / G.h)) for a in range(3))
    z, zt = S.forward_packet(G, x0, k0, 2.0, 1e-3, bg.w[(slice(None),) + i0], 1.0, 1.0, mu)
    f = z.real.copy(); v = zt.real.copy(); Fs = S.Scalar(bg, mu); dt = 0.05
    a = Fs.accel(f)
    for _ in range(int(round(tt / dt))):
        v += .5 * dt * a; f += dt * v; a = Fs.accel(f); v += .5 * dt * a
    h3 = G.h ** 3

    def Lag(gam, N=1.0):
        sg = np.sqrt(np.linalg.det(gam)); A = N * sg * np.linalg.inv(gam)
        w = S.split_weights(A, 12.0)
        Ms = sg / N; Vs = N * sg
        Lk = 0.5 * Ms * (v ** 2).sum() - 0.5 * Vs * mu ** 2 * (f ** 2).sum()
        Lg = sum(0.5 * w[d] * ((S.shift(f, d, 1) - f) ** 2).sum() for d in range(9)) / G.h ** 2
        return (Lk - Lg) * h3
    E = float((Fs.edens(f, v)).sum() * h3)
    T = np.zeros((3, 3))
    for i in range(3):
        for j in range(i, 3):
            dg = np.zeros((3, 3)); dg[i, j] = dg[j, i] = EPS
            dL = (Lag(np.eye(3) + dg) - Lag(np.eye(3) - dg)) / (2 * EPS)
            T[i, j] = T[j, i] = 2 * dL if i == j else dL
    sN = -(Lag(np.eye(3), 1 + EPS) - Lag(np.eye(3), 1 - EPS)) / (2 * EPS)
    # closed-form lattice expectation: v_g * v_p at k0 from the dispersion
    w0 = S.split_weights(np.eye(3), 12.0)
    om = lambda k: float(S.omega_flat(np.array([k, 0, 0]), w0, 1.0, 1.0, mu, G.h))
    vp = om(3.0) / 3.0; vg = (om(3.0 + 1e-5) - om(3.0 - 1e-5)) / 2e-5
    return dict(mu=mu, E=E, sN=float(sN), T=T.tolist(), Txx_over_E=float(T[0, 0] / E),
                vg_vp=vg * vp, vg=vg, vp=vp, Tyy_over_E=float(T[1, 1] / E), Tzz_over_E=float(T[2, 2] / E))


def sources():
    OUT.mkdir(parents=True, exist_ok=True); t0 = time.process_time()
    knots = [knot_case(0.2, 1), knot_case(0.2, 2), knot_case(0.15, 1)]
    packs = []
    for mu in (0.0, 1.0):
        for T in (0.0, 5.0):
            r = packet_case(mu, T); r['t'] = T; packs.append(r)
    import scipy
    res = dict(knots=knots, packets=packs, cpu_s=time.process_time() - t0,
               platform=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__))
    (OUT / 'sources.json').write_text(json.dumps(res, indent=1) + '\n')
    for k in knots:
        print('knot h=%g B=%d E=%.5f sN=%.5f E4-E2=%.5f trT=%.5f Tdiag=%s closed-num=%.1e sha_match=%s' % (
            k['h'], k['B'], k['E'], k['sN'], k['E4_minus_E2'], k['trace_numeric'],
            np.round(np.diag(k['T_numeric']), 5), k['closed_vs_numeric_max'], k['sha_match']))
    for r in packs:
        print('packet mu=%g t=%g E=%.4e sN/E=%.6f Txx/E=%.4f vg*vp=%.4f Tyy/E=%.4f Tzz/E=%.4f' % (
            r['mu'], r['t'], r['E'], r['sN'] / r['E'], r['Txx_over_E'], r['vg_vp'], r['Tyy_over_E'], r['Tzz_over_E']))


if __name__ == '__main__':
    {'states': states, 'sources': sources}[sys.argv[1]]()
