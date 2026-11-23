"""EX-GC34E E1 separate path: kinetic energy summed directly in real space (per-tetrahedron strain-rate solve) on a
periodic 8^3 lattice for random fields, vs the Fourier quadratic form of exgc34e.kinetic_kernel.
  python src/exgc34e_check.py -> results/EX-GC34E/check.json"""
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))
import exgc34d as D
import exgc34e as E

N = 8; rng = np.random.default_rng(3405)
beta = rng.standard_normal((3, N, N, N)); ldot = rng.standard_normal((7, N, N, N)) * E.H
sh = lambda f, off: np.roll(f, shift=tuple(-np.array(off)), axis=(-3, -2, -1))
total = 0.0
for o in D.TETS:
    es = D.tet_edges(o)
    Wm = np.zeros((6, 6)); s = []
    for e, (i, j, t) in enumerate(es):
        Dv = E.ET[t] * E.H; l0 = np.linalg.norm(Dv); Dh = Dv / l0
        for c, (a, b) in enumerate(E.SYM):
            Wm[e, c] = Dv[a] * Dv[b] * (1 if a == b else 2)
        db = sh(beta, o[j]) - sh(beta, o[i])
        s.append(l0 * (sh(ldot[t], o[i]) - np.einsum('a,a...->...', Dh, db)))
    K = np.einsum('cd,d...->c...', np.linalg.inv(Wm), np.stack(s))
    total += (E.VTET / 2 * np.einsum('c...,cd,d...->...', K, E.GHAT, K)).sum()
x = np.zeros((11, N, N, N), complex); x[1:4] = beta; x[4:] = ldot / (-1j)       # w = 1: -i w dl = ldot
xh = np.fft.fftn(x, axes=(1, 2, 3)); kk = 2 * np.pi * np.fft.fftfreq(N)
q = np.stack(np.meshgrid(kk, kk, kk, indexing='ij'), -1)
T = E.kinetic_kernel(q, 1.0)
four = 0.5 * np.einsum('axyz,xyzab,bxyz->', np.conj(xh), T, xh).real / N ** 3
out = dict(direct=float(total), fourier=float(four), rel=float(abs(total / four - 1)))
(ROOT / 'results/EX-GC34E/check.json').write_text(json.dumps(out, indent=1) + '\n'); print(out)
