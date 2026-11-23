"""EX-GC34C numerical check: build n and h_ij on a periodic grid from the solved spectra (coefficients of rho/K^2
from results/EX-GC34C/analysis.json), and compute straight-ray deflection integrals D = int d_x(h_zz/2 - n) dz vs the
lapse-only D0 = int d_x(-n) dz. Expect D/D0 = 1 + gamma.  -> results/EX-GC34C/check.json"""
import json
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / 'results/EX-GC34C'
N, Lb = 128, 64.0; h = Lb / N
k1 = 2 * np.pi * np.fft.fftfreq(N, d=h)
KX, KY, KZ = np.meshgrid(k1, k1, k1, indexing='ij'); K2 = KX ** 2 + KY ** 2 + KZ ** 2; K2[0, 0, 0] = 1
x = (np.arange(N) - N // 2) * h
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
rho = np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / (2 * 0.5 ** 2)); rho /= rho.sum() * h ** 3
rhok = np.fft.fftn(np.fft.ifftshift(rho))


def fields(cn, cA, cB):
    """n^ = cn rho^/K^2 etc.; h_zz = A + B khat_z^2."""
    f = lambda c, extra=1.0: np.fft.fftshift(np.real(np.fft.ifftn(np.where(K2 == 1, 0, c * rhok / K2 * extra) * (K2 != 0))))
    n = f(cn); hzz = f(cA) + f(cB, KZ ** 2 / K2)
    return n, hzz


def deflection(F, b):
    i = N // 2 + int(round(b / h)); j = N // 2
    dFdx = (np.roll(F, -1, 0) - np.roll(F, 1, 0)) / (2 * h)
    return dFdx[i, j, :].sum() * h


def coeffs(expr_n, expr_A, expr_B, subs):
    vals = {str(k): v for k, v in subs.items()}; vals.update(K=1, rho=1)
    def g(e):
        if e is None:
            return 0.0
        ex = sp.sympify(e)
        return float(ex.subs({sym: vals[sym.name] for sym in ex.free_symbols}))
    return g(expr_n), g(expr_A), g(expr_B)


def main():
    an = json.loads((OUT / 'analysis.json').read_text()); s = sp.Symbol('s', real=True)
    a, b, c1, c2, c3, c4 = sp.symbols('a b c1 c2 c3 c4', real=True)
    g2 = an['A2_general'][0]
    cases = {
        'A0': (coeffs(an['A0']['n'], an['A0']['A'], an['A0']['B'], {}), 0.0),
        'A1': (coeffs(an['A1']['n'], an['A1']['A'], None, {}), -1.0),
        'A3': (coeffs(an['A3'][0]['n'], an['A3'][0]['A'], None, {s: sp.Rational(1, 4)}), 1.0),
    }
    pts = {'A2_generic': {a: 1, b: 0.3, c1: 0.7, c2: 0.2, c3: -0.1, c4: 0.05}}
    gs = an['A2_gauge_solution']; sub = {b: 1, c1: 1, c4: 0}
    ex = sp.sympify(an['A2_gauge_gamma']); vv = {str(k): v for k, v in sub.items()}
    cases['A2_gaugeinv_other_ratio'] = (coeffs(gs['n'], gs['A'], None, sub),
                                        float(ex.subs({sym: vv[sym.name] for sym in ex.free_symbols})))
    for name, sub in pts.items():
        ex = sp.sympify(an['A2_gamma']); vv = {str(k): v for k, v in sub.items()}
        gam = float(ex.subs({sym: vv[sym.name] for sym in ex.free_symbols}))
        cases[name] = (coeffs(g2['n'], g2['A'], g2['B'], sub), gam)
    res = {}
    for name, ((cn, cA, cB), gam) in cases.items():
        n, hzz = fields(cn, cA, cB)
        rows = []
        for bb in (8, 12, 16):
            D = deflection(0.5 * hzz - n, bb); D0 = deflection(-n, bb)
            rows.append(dict(b=bb, ratio=D / D0, expected=1 + gam, rel_err=(D / D0) / (1 + gam) - 1 if gam != -1 else D / D0))
        res[name] = dict(coeffs=[cn, cA, cB], gamma=gam, rows=rows)
        print(name, 'gamma', round(gam, 6), [(r['b'], round(r['ratio'], 6)) for r in rows])
    (OUT / 'check.json').write_text(json.dumps(res, indent=1) + '\n')


if __name__ == '__main__':
    main()
