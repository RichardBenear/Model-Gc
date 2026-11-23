"""EX-GC34D part D5: is a shift-like variable forced for a moving source? (linearised GR, sympy, Fourier space)

Plane wave h_{mu nu} exp(i(k.x - w t)), k = (k, 0, 0), dust moving with v = (vx, vy, 0), w = k vx (conservation).
Full theory: all ten G_{mu nu} = kappa T_{mu nu}. Reduced ("no shift"): h_{0i} = 0, only the lapse (00) and the
six ij equations. Comparison through the gauge-invariant linearised Riemann tensor.

  python src/exgc34d_shift.py      -> results/EX-GC34D/shift.json
"""
import json, time
from pathlib import Path
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results/EX-GC34D'

k, kap, rho = sp.symbols('k kappa rho', positive=True)
vx, vy = sp.symbols('v_x v_y', real=True)
eta = sp.diag(-1, 1, 1, 1)


def setup(w):
    kl = sp.Matrix([-w, k, 0, 0])                 # k_mu (lower): d_mu -> i k_mu, d_t -> -i w
    ku = eta * kl
    k2 = (kl.T * ku)[0]
    v = [vx, vy, 0]
    T = sp.zeros(4, 4); T[0, 0] = rho
    for i in range(3):
        T[0, i + 1] = T[i + 1, 0] = -rho * v[i]
        for j in range(3):
            T[i + 1, j + 1] = rho * v[i] * v[j]
    return kl, ku, k2, T


def einstein(hm, kl, ku, k2):
    tr = sum(eta[m, m] * hm[m, m] for m in range(4))
    kh = hm * ku                                   # k^l h_{l nu}
    khk = (ku.T * hm * ku)[0]
    G = sp.zeros(4, 4)
    for m in range(4):
        for n in range(4):
            G[m, n] = sp.Rational(1, 2) * (k2 * hm[m, n] - kl[m] * kh[n] - kl[n] * kh[m] + kl[m] * kl[n] * tr
                                           + eta[m, n] * khk - eta[m, n] * k2 * tr)
    return G


def riemann(hm, kl):
    """Linearised R_{abcd} in Fourier space (common factor dropped); gauge invariant."""
    R = {}
    for a in range(4):
        for b in range(a + 1, 4):
            for c in range(4):
                for d in range(c + 1, 4):
                    R[(a, b, c, d)] = sp.Rational(1, 2) * (kl[b] * kl[c] * hm[a, d] + kl[a] * kl[d] * hm[b, c]
                                                           - kl[a] * kl[c] * hm[b, d] - kl[b] * kl[d] * hm[a, c])
    return R


def solve_case(w, reduced):
    kl, ku, k2, T = setup(w)
    syms = {}
    hm = sp.zeros(4, 4)
    for m in range(4):
        for n in range(m, 4):
            if reduced and m == 0 and n > 0:
                continue
            s = sp.Symbol(f'h{m}{n}'); syms[(m, n)] = s
            hm[m, n] = hm[n, m] = s
    G = einstein(hm, kl, ku, k2)
    idx = [(m, n) for m in range(4) for n in range(m, 4)]
    if reduced:
        idx = [p for p in idx if not (p[0] == 0 and p[1] > 0)]
    eqs = [sp.expand(G[m, n] - kap * T[m, n]) for m, n in idx]
    sol = sp.solve(eqs, list(syms.values()), dict=True)
    return hm, G, T, kl, sol, syms


def main():
    t0 = time.process_time(); res = {}
    w = k * vx
    # full and reduced, moving source with w != 0
    hf, Gf, Tf, kl, solf, _ = solve_case(w, False)
    hr, Gr, Tr, _, solr, symr = solve_case(w, True)
    res['full_solvable'] = bool(solf); res['reduced_solvable'] = bool(solr)
    hf_s = hf.subs(solf[0]); hr_s = hr.subs(solr[0])
    res['reduced_free_params'] = sorted(str(s) for s in hr_s.free_symbols - {k, kap, rho, vx, vy})
    Rf = riemann(hf_s, kl); Rr = riemann(hr_s, kl)
    diff = {str(key): sp.simplify(Rr[key] - Rf[key]) for key in Rf}
    res['riemann_reduced_minus_full_nonzero'] = {kk: str(v) for kk, v in diff.items() if v != 0}
    mom = [sp.simplify((Gr[0, i] - kap * Tr[0, i]).subs(solr[0])) for i in range(1, 4)]
    res['reduced_momentum_constraint_residual'] = [str(m) for m in mom]
    # vector part of h_ij in the reduced solution: h_xy (k along x, v_y transverse)
    res['reduced_h_xy'] = str(sp.factor(hr_s[1, 2]))
    res['reduced_h_xy_times_w'] = str(sp.factor(sp.simplify(hr_s[1, 2] * w)))
    res['reduced_n'] = str(sp.factor(-hr_s[0, 0] / 2))
    # static-like mode: k perpendicular to v (w = 0, v = (0, vy, 0))
    hs, Gs, Ts, kls, sols, _ = solve_case(sp.Integer(0), True)
    sub0 = {vx: 0}
    Gs0 = Gs.subs(sub0); Ts0 = Ts.subs(sub0)
    eqs0 = [sp.expand(Gs0[m, n] - kap * Ts0[m, n]) for m in range(4) for n in range(m, 4) if not (m == 0 and n > 0)]
    syms0 = sorted(hs.free_symbols - {k, kap, rho, vx, vy}, key=str)
    sol0 = sp.solve(eqs0, syms0, dict=True)
    res['k_perp_v_reduced_solvable'] = bool(sol0)
    if sol0:
        res['k_perp_v_momentum_constraint_residual'] = [str(sp.simplify((Gs0[0, i] - kap * Ts0[0, i]).subs(sol0[0])))
                                                        for i in range(1, 4)]
    # check of the operator: static Newtonian limit of the full theory (v = 0): n = -kappa rho / (2 k^2)
    hN, GN, TN, klN, solN, _ = solve_case(sp.Integer(0), False)
    hN_s = hN.subs(solN[0]).subs({vx: 0, vy: 0})
    RN = riemann(hN_s, klN)
    res['static_R0101'] = str(sp.factor(RN[(0, 1, 0, 1)]))   # = k^2 n for h00 = -2n -> expect -kappa rho / 2
    res['cpu_s'] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / 'shift.json').write_text(json.dumps(res, indent=1) + '\n')
    for kk, v in res.items():
        print(kk, ':', v)


if __name__ == '__main__':
    main()
