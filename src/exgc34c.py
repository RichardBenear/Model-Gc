"""EX-GC34C (exploratory): static linear gravity actions on the shared structure -> PPN gamma.

Scalar sector in Fourier space: h_ij = A delta_ij + B khat_i khat_j, lapse n, source rho (T^ij = 0).
Substitution rules for a real Fourier mode (common factor 1/2 dropped): d_a X d_b Y -> k_a k_b X Y,
n d_i d_j h -> -k_i k_j n h. gamma = -A/(2n) at k -> 0 (B drops out of straight-ray deflection).

  python src/exgc34c.py      -> results/EX-GC34C/analysis.json (+ printed summary)
"""
import json, platform, time
from pathlib import Path
import numpy as np
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34C'

K, n, A, B, rho = sp.symbols('K n A B rho', real=True)
a, b, c1, c2, c3, c4, s = sp.symbols('a b c1 c2 c3 c4 s', real=True)
hh = 3 * A ** 2 + 2 * A * B + B ** 2            # h_ij h_ij
kh = A + B                                        # khat_i khat_j h_ij
tr = 3 * A + B                                    # h_kk
khkh = (A + B) ** 2                               # (khat_i h_ij)(khat_k h_kj)


def solve(L, fields):
    eqs = [sp.diff(L, f) for f in fields]
    sol = sp.solve(eqs, fields, dict=True)
    return sol


def pure_clock_L():
    Lam = -K ** 2 * (a * kh + b * tr)
    Q = K ** 2 * (c1 * hh + c2 * khkh + c3 * kh * tr + c4 * tr ** 2)
    return -(n * Lam + Q) - n * rho


def fp_static_L():
    """Fierz-Pauli with h00 = -2n, h0i = 0, static; h = 2n + tr."""
    h = 2 * n + tr
    L = (-sp.Rational(1, 2) * K ** 2 * (4 * n ** 2 + hh) + K ** 2 * khkh
         - K ** 2 * kh * h + sp.Rational(1, 2) * K ** 2 * h ** 2)
    return sp.expand(s * L) - n * rho


def gamma_of(sol):
    return sp.simplify(-sol[A] / (2 * sol[n]))


def main():
    t0 = time.process_time(); res = {}
    # A0: independent fields
    L0 = -K ** 2 * n ** 2 - sp.Rational(1, 4) * K ** 2 * hh - n * rho
    s0 = solve(L0, [n, A, B])[0]
    res['A0'] = dict(n=str(s0[n]), A=str(s0[A]), B=str(s0[B]), gamma=str(gamma_of(s0)))
    # A1: Nordstrom (n = phi, A = 2 phi, B = 0), L = -K^2 phi^2 - phi rho
    phi = sp.symbols('phi'); s1 = sp.solve(sp.diff(-K ** 2 * phi ** 2 - phi * rho, phi), phi)[0]
    res['A1'] = dict(n=str(s1), A=str(2 * s1), gamma=str(sp.simplify(-2 * s1 / (2 * s1))))
    # A3: static Fierz-Pauli
    L3 = fp_static_L()
    nn = sp.Poly(sp.expand(L3), n).coeff_monomial(n ** 2)
    res['A3_n2_coefficient'] = str(sp.simplify(nn))                # should vanish: lapse is a pure multiplier
    s3 = solve(L3, [n, A, B])
    res['A3_solutions'] = [{str(k): str(v) for k, v in d.items()} for d in s3]
    # FP has a residual gauge freedom in B (spatial diffeo): solve with B free
    s3b = sp.solve([sp.diff(L3, n), sp.diff(L3, A)], [n, A], dict=True)
    res['A3'] = [{str(k): str(sp.simplify(v)) for k, v in d.items()} for d in s3b]
    if s3b:
        g3 = sp.simplify(-s3b[0][A] / (2 * s3b[0][n]))
        res['A3_gamma'] = str(g3)
        res['A3_dL_dB_residual'] = str(sp.simplify(sp.diff(L3, B).subs(s3b[0])))
        res['A3_n_over_A0_n'] = str(sp.simplify(s3b[0][n] / s0[n]))
    # A2: general pure clock
    L2 = pure_clock_L()
    s2 = sp.solve([sp.diff(L2, f) for f in (n, A, B)], [n, A, B], dict=True)
    res['A2_general'] = [{str(k): str(sp.factor(v)) for k, v in d.items()} for d in s2]
    if s2:
        g2 = sp.factor(-s2[0][A] / (2 * s2[0][n]))
        res['A2_gamma'] = str(g2)
        # identify FP point inside the A2 family: Lambda and Q read off from L3
        res['A2_gamma_at_FPlike'] = None
    # FP coefficients expressed in A2 language: expand L3 + n rho and match
    L3m = sp.expand(L3 + n * rho)
    lin = sp.expand(sp.diff(L3m, n))                    # = -Lam  (since n appears linearly)
    res['A3_Lambda_fourier'] = str(sp.factor(-lin))
    Q3 = sp.expand(-(L3m.subs(n, 0)))
    res['A3_Q_fourier'] = str(sp.factor(Q3))
    # A2 restricted to linear gauge invariance (no dependence on B): a + b = 0, Q independent of B
    Lg = pure_clock_L().subs(a, -b)
    Qc = sp.Poly(sp.expand(K ** 2 * (c1 * hh + c2 * khkh + c3 * kh * tr + c4 * tr ** 2) / K ** 2), A, B)
    cons = [Qc.coeff_monomial(A * B), Qc.coeff_monomial(B ** 2)]
    csol = sp.solve(cons, [c2, c3], dict=True)[0]
    Lg = sp.expand(Lg.subs(csol))
    res['A2_gauge_constraints'] = {str(k): str(v) for k, v in csol.items()}
    res['A2_gauge_dL_dB'] = str(sp.simplify(sp.diff(Lg, B)))
    sg = sp.solve([sp.diff(Lg, n), sp.diff(Lg, A)], [n, A], dict=True)[0]
    res['A2_gauge_solution'] = {str(k): str(sp.factor(v)) for k, v in sg.items()}
    res['A2_gauge_gamma'] = str(sp.factor(-sg[A] / (2 * sg[n])))
    res['cpu_s'] = time.process_time() - t0
    res['platform'] = dict(python=platform.python_version(), sympy=sp.__version__)
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'analysis.json').write_text(json.dumps(res, indent=1) + '\n')
    for k, v in res.items():
        print(k, ':', v)


if __name__ == '__main__':
    main()
