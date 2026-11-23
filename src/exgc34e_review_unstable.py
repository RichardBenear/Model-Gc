"""EX-GC34E post hoc review check (2026-11-23): the vacuum kernel has a growing (imaginary-frequency) mode.

Prompted by ChatGPT's review of white paper v0.5. Not a re-scoring of the registered E3 criteria (P1 stays FAIL as
registered). Uses the committed 34E kernel (src/exgc34e.py); the analysis path here is separate from wave_modes():
  1. roots of the quadratic pencil (axis edges dropped = displacement gauge), each substituted back into ALL 11
     equations of M(w) = C + w B + w^2 A (A, B, C from M at w = 0, +1, -1; M is exactly quadratic in w);
  2. independently, the smallest singular value of the gauge-fixed M(w) minimised along the imaginary-w axis (one more
     singular value is skipped where the pencil has an exact null direction);
  3. how much of the growing mode's metric perturbation is pure gauge (continuum h -> h + d xi + d xi), and its
     field content (lapse, shift, edges).
Convention e^{-i w t}: +Im w grows.

  python src/exgc34e_review_unstable.py   -> results/EX-GC34E/review_unstable.json
"""
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(k, '1')
import json, platform, sys, time
from pathlib import Path
import numpy as np
import mpmath as mp
import scipy.linalg as sl
import scipy.optimize as so

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
import exgc34e as E

H = E.H
KEEP = [0, 1, 2, 3] + [4 + t for t in range(3, 7)]
DIRS = {'100': (1, 0, 0), '110': (1, 1, 0), '111': (1, 1, 1), 'generic': (0.8, 0.5, 0.33), 'near100': (1, 0.05, 0)}
IU = np.triu_indices(4)


def pencil(st, q):
    Ms = [E.full_kernel(st, q[None], w)[0] for w in (0.0, 1.0, -1.0)]
    C = Ms[0]
    return C, (Ms[1] - Ms[2]) / 2, (Ms[1] + Ms[2] - 2 * C) / 2


def pencil_roots(C, B, A, k):
    Ck, Bk, Ak = [X[np.ix_(KEEP, KEEP)] for X in (C, B, A)]
    n = len(KEEP); I = np.eye(n); Z = np.zeros((n, n))
    ev, V = sl.eig(np.block([[Z, I], [-Ck, -Bk]]), np.block([[I, Z], [Z, Ak]]))
    out = []
    for e, v in zip(ev, V.T):
        if not np.isfinite(e) or abs(e) / k > 20:
            continue
        x = np.zeros(11, complex); x[KEEP] = v[:n]; x /= np.linalg.norm(x)
        M = C + e * B + e * e * A
        out.append(dict(w_over_k=[float(e.real / k), float(e.imag / k)],
                        resid_all11=float(np.linalg.norm(M @ x) / np.linalg.norm(M, 2))))
    return out


def sv(C, B, A, w, idx):
    s = np.linalg.svd((C + w * B + w * w * A)[np.ix_(KEEP, KEEP)], compute_uv=False)
    return s[-idx] / s[0]


def null_vec(C, B, A, w, idx):
    U, s, Vh = np.linalg.svd((C + w * B + w * w * A)[np.ix_(KEEP, KEEP)])
    x = np.zeros(11, complex); x[KEEP] = Vh[-idx].conj()
    return x


def main():
    mp.mp.dps = 50
    t0 = time.process_time()
    st = E.tet_static()
    res = dict(convention='e^{-i w t}; +Im w grows', pencil_roots=[], imag_axis=[], real_axis_floor=[])
    for kh in (0.025, 0.05, 0.1, 0.2):
        k = kh / H
        for dn, d in DIRS.items():
            q = kh * np.array(d, float) / np.linalg.norm(d)
            C, B, A = pencil(st, q)
            res['pencil_roots'].append(dict(kh=kh, dir=dn, roots=pencil_roots(C, B, A, k)))
            # exact null directions: singular values that stay tiny at generic real w
            floor = [[float(sv(C, B, A, g * k, i)) for i in (1, 2, 3)] for g in (0.37, 0.81, 1.3, 2.2)]
            nexact = int(sum(max(f[i] for f in floor) < 1e-15 for i in range(3)))
            res['real_axis_floor'].append(dict(kh=kh, dir=dn, sv123_at_generic_real_w=floor, n_exact_null=nexact))
            if dn == '100':
                continue                       # three exact null directions: this method cannot resolve 100
            idx = nexact + 1
            br = (2.0, 2.3) if dn == '110' else ((1.45, 1.75) if dn == 'near100' else (1.3, 1.55))
            r = so.minimize_scalar(lambda g: np.log(sv(C, B, A, 1j * g * k, idx)), bounds=br, method='bounded',
                                   options=dict(xatol=1e-10))
            g = float(r.x); w = 1j * g * k
            x = null_vec(C, B, A, w, idx)
            M = C + w * B + w * w * A
            h = E.metric_from_x(x, q); kl = np.r_[-w, q / H]
            Gb = []
            for m in range(4):
                hg = np.zeros((4, 4), complex); hg[m, :] += 1j * kl; hg[:, m] += 1j * kl; Gb.append(hg[IU])
            Gb = np.array(Gb).T; hv = h[IU]
            c, *_ = np.linalg.lstsq(Gb, hv, rcond=None)
            qh = q / np.linalg.norm(q)
            res['imag_axis'].append(dict(
                kh=kh, dir=dn, sv_index=idx, growth_Im_w_over_k=g, sv_at_root=float(np.exp(r.fun)),
                sv_at_bracket_ends=[float(sv(C, B, A, 1j * b * k, idx)) for b in br],
                resid_all11=float(np.linalg.norm(M @ x) / np.linalg.norm(M, 2)),
                abs_lapse=float(abs(x[0])), abs_shift_long=float(abs(qh @ x[1:4])),
                abs_shift_trans=float(np.linalg.norm(x[1:4] - qh * (qh @ x[1:4]))),
                abs_diag_edges=float(np.linalg.norm(x[7:])),
                non_gauge_fraction_of_h=float(np.linalg.norm(hv - Gb @ c) / np.linalg.norm(hv))))
            a = res['imag_axis'][-1]
            print(kh, dn, 'Im w/k = %.5f  sv %.1e (ends %.1e %.1e)  resid %.1e  non-gauge %.2e' % (
                g, a['sv_at_root'], *a['sv_at_bracket_ends'], a['resid_all11'], a['non_gauge_fraction_of_h']), flush=True)
    res['light_speed_kh005'] = {dn: float(E.light_speed(0.05 * np.array(d, float) / np.linalg.norm(d)))
                                for dn, d in DIRS.items()}
    res['platform'] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    res['cpu_s'] = time.process_time() - t0
    (ROOT / 'results/EX-GC34E/review_unstable.json').write_text(json.dumps(res, indent=1))
    print('cpu', res['cpu_s'])


if __name__ == '__main__':
    main()
