"""EX-GC26F pre-plan check (Claude): static Skyrme interaction of two B = 1 knots (product ansatz of the closed 26C
rest knot) vs separation d and relative orientation, against the predicted gravitational force at G_N = 0.01."""
import sys, json, hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
import exgc26e1 as e1, exgc26d as d26
e1.configure(0.2); e1.load()
man = json.loads((d26.OUT / 'inputs.json').read_text()); p1 = np.load(d26.RAW / 'c26_p.npy')
assert hashlib.sha256((d26.RAW / 'c26_p.npy').read_bytes()).hexdigest() == man['c26_p.npy']
G1 = e1.ref.Grid(0.2, 7.0)
def qmul(a, b):
    a0, av = a[0], a[1:]; b0, bv = b[0], b[1:]
    return np.concatenate([[a0 * b0 - np.sum(av * bv, 0)], a0 * bv + b0 * av + np.cross(av, bv, axis=0)])
def place(G, c, R=None):
    X = np.array(G.X, float); X[0] -= c; idx = X / G1.h + (G1.n - 1) / 2
    q = np.stack([map_coordinates(p1[k], idx, order=1, mode='constant', cval=(1. if k == 0 else 0.)) for k in range(4)])
    if R is not None: q[1:] = np.einsum('ij,j...->i...', R, q[1:])      # isorotation
    return q / np.sqrt(np.sum(q * q, 0))
Rz = np.diag([-1., -1., 1.])   # pi about isospin z (perpendicular to separation x)
Rx = np.diag([1., -1., -1.])   # pi about isospin x (along separation)
out = {}
for L in (14.0,):
    G = e1.ref.Grid(0.2, L); E1 = float(e1.evaluate(G, G.fix(place(G, 0.)))[0].sum()); print("E1 parts", e1.evaluate(G, G.fix(place(G, 0.)))[0])
    out['E_single'] = E1
    for dd in (5., 6., 7., 8., 10., 12.):
        row = {}
        for name, R in (('same', None), ('pi_perp', Rz), ('pi_par', Rx)):
            p = G.fix(qmul(place(G, -dd / 2), place(G, dd / 2, R))); row[name] = float(e1.evaluate(G, p)[0].sum()) - 2 * E1
        out[str(dd)] = row; print(dd, row, flush=True)
Path(__file__).with_name('pair_static.json').write_text(json.dumps(out, indent=1) + '\n')
