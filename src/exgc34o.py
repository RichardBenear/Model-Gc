"""EX-GC34O (exploratory): a Y1 scalar stress that conserves momentum exactly. Holds 34N's Y1 layout, flat scalar
evolution and lapse/shift couplings; the parallel normal-metric coupling becomes 1/4 h_aa(v) G_a(v) G_a(v-e_a) at
sites. Plan: spec/EX-GC34O_plan.md.

  python src/exgc34o.py run  -> results/EX-GC34O/results.json  (O0-O3, null control, reported items)
Reuses exgc34n (fw, grads, msources, mlagrangian, msources_dot, midentities, PAIRS) and exgc34m.packet / continuum data.
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34n as N
import exgc34m as Mm

OUT = ROOT / 'results/EX-GC34O'
fw = N.fw; PAIRS = N.PAIRS; MASS = N.MASS; DUMP = N.DUMP


# ---------- Y1-S2 sources and action ----------
def msources2(phi, phid, h):
    m2 = MASS ** 2; V = h ** 3; G = N.grads(phi, h)
    s_n = V * (-0.5 * phid ** 2 - 0.5 * m2 * phi ** 2)
    s_d = [V * (0.25 * phid ** 2 - 0.25 * m2 * phi ** 2) for _ in range(3)]
    s_b = []
    for a in range(3):
        g2 = G[a] ** 2; g2s = g2 + fw(g2, a, -1)          # links v and v - e_a
        s_n = s_n - V * 0.25 * g2s                         # lapse coupling unchanged (34N)
        for c in range(3):
            if c != a:
                s_d[c] = s_d[c] - V * 0.125 * g2s          # -1/4 h-bar_cc G_a^2, c != a
        s_d[a] = s_d[a] + V * 0.25 * G[a] * fw(G[a], a, -1)   # +1/4 h_aa(v) G_a(v) G_a(v-e_a)
        s_b.append(-V * 0.5 * (phid + fw(phid, a)) * G[a])
    s_o = []
    for (a, b) in PAIRS[3:]:
        Ga = 0.5 * (G[a] + fw(G[a], b)); Gb = 0.5 * (G[b] + fw(G[b], a))
        s_o.append(V * Ga * Gb)
    return s_n, np.array(s_b), np.array(s_d), np.array(s_o)


def mlagrangian2(phi, phid, n, b, hd, ho, h):
    """Coded directly from the plan's action, independently of msources2."""
    m2 = MASS ** 2; V = h ** 3; G = N.grads(phi, h); T = hd.sum(0)
    L = np.sum(0.5 * (1 - n + 0.5 * T) * phid ** 2 - 0.5 * (1 + n + 0.5 * T) * m2 * phi ** 2)
    for a in range(3):
        nb = 0.5 * (n + fw(n, a))
        Tc = sum(0.5 * (hd[c] + fw(hd[c], a)) for c in range(3) if c != a)
        L += np.sum(-0.5 * (1 + nb + 0.5 * Tc) * G[a] ** 2 - 0.5 * (phid + fw(phid, a)) * b[a] * G[a])
        L += np.sum(0.25 * hd[a] * G[a] * fw(G[a], a, -1))
    for c, (a, bb) in enumerate(PAIRS[3:]):
        Ga = 0.5 * (G[a] + fw(G[a], bb)); Gb = 0.5 * (G[bb] + fw(G[bb], a))
        L += np.sum(ho[c] * Ga * Gb)
    return V * L


def sources_dot(fn, phi, phid, h):
    phidd = -MASS ** 2 * phi + N.lap(phi, h)
    sp = fn(phi + phid, phid + phidd, h); sm = fn(phi - phid, phid - phidd, h)
    return [(x - y) / 2 for x, y in zip(sp, sm)]


def rel(fn, phi, phid, h):
    s = fn(phi, phid, h); sd = sources_dot(fn, phi, phid, h)
    IE, IP, div = N.midentities(s, sd, h)
    rE = float(np.abs(IE).max() / np.abs(sd[0]).max())
    rP = float(np.abs(IP).max() / max(np.abs(sd[1]).max(), np.abs(div).max()))
    return rE, rP, s


# ---------- O0: 1D chain ----------
def chain(rng, Nn, h, form):
    m2 = MASS ** 2
    phi = rng.standard_normal(Nn); phid = rng.standard_normal(Nn)
    r = lambda a, s=1: np.roll(a, -s)
    G = (r(phi) - phi) / h
    phidd = (r(phi) + r(phi, -1) - 2 * phi) / h ** 2 - m2 * phi
    Gd = (r(phid) - phid) / h
    pdot = -0.5 * (phidd + r(phidd)) * G - 0.5 * (phid + r(phid)) * Gd    # d/dt of p = -phidot-bar G on links
    Gm = r(G, -1)
    T = 0.5 * (phid ** 2 - m2 * phi ** 2) + (0.5 * G * Gm if form == 'product' else 0.25 * (G ** 2 + Gm ** 2))
    res = pdot + (r(T) - T) / h
    return float(np.abs(res).max() / np.abs(pdot).max())


def run():
    t0 = time.process_time(); rng = np.random.default_rng(3434 + 16); res = {}
    o0 = [dict(h=h, product=chain(rng, 64, h, 'product'), avgsq=chain(rng, 64, h, 'avgsq')) for h in (1.0, 0.5)]
    res['O0_rows'] = o0
    res['O0'] = bool(max(r['product'] for r in o0) <= 1e-12 and min(r['avgsq'] for r in o0) >= 1e-3)
    print('O0', res['O0'], o0, flush=True)

    o1 = []; rnd = []
    for h in (1.0, 0.5):
        n8 = 8; phi = rng.standard_normal((n8,) * 3); phid = rng.standard_normal((n8,) * 3)
        n = rng.standard_normal((n8,) * 3); b = rng.standard_normal((3,) + (n8,) * 3)
        hd = rng.standard_normal((3,) + (n8,) * 3); ho = rng.standard_normal((3,) + (n8,) * 3)
        s_n, s_b, s_d, s_o = msources2(phi, phid, h)
        pred = np.sum(s_n * n) + np.sum(s_b * b) + np.sum(s_d * hd) + np.sum(s_o * ho)

        def cd(eps):
            return (mlagrangian2(phi, phid, eps * n, eps * b, eps * hd, eps * ho, h)
                    - mlagrangian2(phi, phid, -eps * n, -eps * b, -eps * hd, -eps * ho, h)) / (2 * eps)
        o1.append(dict(h=h, rel=float(abs((4 * cd(1e-4) - cd(2e-4)) / 3 - pred) / abs(pred))))
        rE, rP, _ = rel(msources2, phi, phid, h); rE0, rP0, _ = rel(N.msources, phi, phid, h)
        rnd.append(dict(h=h, rE=rE, rP=rP, rE_34N=rE0, rP_34N=rP0))
    res['O1_rows'] = o1; res['O1'] = bool(max(r['rel'] for r in o1) <= 1e-8)
    res['random_rows'] = rnd
    print('O1', res['O1'], o1, '\nrandom', rnd, flush=True)

    rows = []; integ = []
    old = json.loads((ROOT / 'results/EX-GC34N/matter.json').read_text())['rows']
    for dname in Mm.DIRS:
        for h in (1.0, 0.5, 0.25):
            Nn = int(round(Mm.LBOX / h))
            for t in (0.0, 6.0, 12.0):
                phi, phid, KK = Mm.packet(Nn, h, dname, t)
                rE, rP, s = rel(msources2, phi, phid, h); rE0, rP0, s0 = rel(N.msources, phi, phid, h)
                reg = next(x['rP'] for x in old if x['dir'] == dname and x['h'] == h and x['t'] == t)
                rows.append(dict(dir=dname, h=h, t=t, rE=rE, rP=rP, rP_34N_here=rP0, rP_34N_registered=reg))
                if t == 6.0:
                    integ.append(dict(dir=dname, h=h, S2=s[2].reshape(3, -1).sum(1), N34=s0[2].reshape(3, -1).sum(1),
                                      shear=s[3].reshape(3, -1).sum(1)))
                print(dname, h, t, 'rE %.1e rP %.1e | 34N rP %.3e (reg %.3e) cpu %.0f'
                      % (rE, rP, rP0, reg, time.process_time() - t0), flush=True)
    res['rows'] = rows
    allE = [r['rE'] for r in rows] + [r['rE'] for r in rnd]; allP = [r['rP'] for r in rows] + [r['rP'] for r in rnd]
    res['O2'] = bool(max(allE) <= 1e-10)
    res['O3'] = 'EXACT' if max(allP) <= 1e-10 else 'FAIL'
    res['null'] = bool(all(abs(r['rP_34N_here'] / r['rP_34N_registered'] - 1) < 5e-3 for r in rows))
    res['overall'] = bool(res['O0'] and res['O1'] and res['O2'] and res['O3'] == 'EXACT' and res['null'])
    # reported: integrated normal stress, product vs 34N form, and the continuum (from 34M's continuum helper if present)
    res['integrated_stress_t6'] = integ
    # reported: sign of the product-form parallel gradient energy at the zone boundary (1D chain)
    x = np.arange(16); phi = (-1.0) ** x; G = np.roll(phi, -1) - phi
    res['zone_boundary_1D'] = dict(P=float((G * np.roll(G, 1)).mean()), X=float((0.5 * (G ** 2 + np.roll(G, 1) ** 2)).mean()))
    res['cpu_s'] = time.process_time() - t0
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'results.json').write_text(json.dumps(res, **DUMP))
    print({k: res[k] for k in ('O0', 'O1', 'O2', 'O3', 'null', 'overall', 'zone_boundary_1D', 'cpu_s')})


if __name__ == '__main__':
    if sys.argv[1] == 'run':
        run()


def recompute():
    """Separate path: d/dt of the link momentum from phi-ddot directly (no polarisation, no midentities), and the
    stress divergence written out from the plan's T_aa and shear formulas, on the h = 0.5 packets and random fields."""
    rng = np.random.default_rng(99); out = []
    cases = [(d, *Mm.packet(int(round(Mm.LBOX / 0.5)), 0.5, d, 6.0)[:2]) for d in Mm.DIRS]
    cases.append(('random', rng.standard_normal((8,) * 3), rng.standard_normal((8,) * 3)))
    h = 0.5; m2 = MASS ** 2
    for name, phi, phid in cases:
        phidd = N.lap(phi, h) - m2 * phi
        G = [(np.roll(phi, -1, a) - phi) / h for a in range(3)]; Gd = [(np.roll(phid, -1, a) - phid) / h for a in range(3)]
        worst = 0.0
        for a in range(3):
            pdot = -0.5 * (phidd + np.roll(phidd, -1, a)) * G[a] - 0.5 * (phid + np.roll(phid, -1, a)) * Gd[a]
            Taa = 0.5 * (phid ** 2 - m2 * phi ** 2) + 0.5 * G[a] * np.roll(G[a], 1, a)
            for b in range(3):
                if b != a:
                    Taa = Taa - 0.25 * (G[b] ** 2 + np.roll(G[b], 1, b) ** 2)
            div = (np.roll(Taa, -1, a) - Taa) / h
            for b in range(3):
                if b != a:
                    Y = 0.25 * (G[a] + np.roll(G[a], -1, b)) * (G[b] + np.roll(G[b], -1, a))   # plaquette v+(a+b)/2
                    div = div + (Y - np.roll(Y, 1, b)) / h
            worst = max(worst, float(np.abs(pdot + div).max() / np.abs(pdot).max()))
        out.append(dict(case=name, rel=worst))
    print(out)
    r = json.loads((OUT / 'results.json').read_text()); r['recompute'] = out
    (OUT / 'results.json').write_text(json.dumps(r, **DUMP))


if __name__ == '__main__' and sys.argv[1] == 'recompute':
    recompute()
