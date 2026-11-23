"""EX-GC34D POST HOC (not registered): D4 rerun with the gauge part of the solved spatial metric removed.
The post hoc conformal fit shows the Regge solution is exactly l_e = l0_e (1 - mean lapse) plus a lattice gauge
(vertex displacement) mode, so the same geometry in isotropic coordinates is gamma_ij = (1 - 2n) delta_ij with the
solved n. This checks whether the registered D4 excess over 2 is a coordinate (gauge) effect of the finite
early/late readout windows.   -> results/EX-GC34D/posthoc_iso.json"""
import json, sys
from pathlib import Path
from types import SimpleNamespace
from multiprocessing import Pool
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT / 'src'))
import exgc34d_bend as Bd
S = Bd.S; R = Bd.R


def iso(n):
    return -2 * n[..., None, None] * np.eye(3)


class IsoBg(Bd.InterpBg):
    def coeffs(self, x):
        shp = x.shape[:-1]
        idx = [((x[..., a] - Bd.LO[a]) / Bd.H).ravel() for a in range(3)]
        n = Bd.ndimage.map_coordinates(self.c[0], idx, order=3, mode='grid-wrap', prefilter=False).reshape(shp)
        N, sg, A = Bd.A_of(n, iso(n), True)
        return sg / N, N * sg, R.weights_from_A(A, 12.0)


def ray_iso(_):
    rec = json.loads((Bd.OUT / 'D4_regge_iso_F1.json').read_text())
    return R.rays(IsoBg(True), Bd.X0, Bd.K0, Bd.SIG, 0.0, Bd.H, Bd.T, Bd.DT, nq=4,
                  times=np.array(rec['hist']['t']), box=Bd.BOX).tolist()


def main():
    z = np.load(Bd.RAW / 'gauss_fields.npz'); nh, hh = z['nh'], z['hh']
    G = S.Grid(Bd.NM2, Bd.LO, Bd.H)
    n0, _ = Bd.fields_at(nh, hh, (0, 0, 0))
    nl = [Bd.fields_at(nh, hh, 0.5 * S.DIRS[d])[0] for d in range(9)]
    N, sg, _ = Bd.A_of(n0, iso(n0), True)
    w = np.stack([S.split_weights(Bd.A_of(nl[d], iso(nl[d]), True)[2], 12.0)[d] for d in range(9)])
    bg = SimpleNamespace(G=G, Ms=sg / N, Vs=N * sg, w=w)
    tag = 'D4_regge_iso_F1'
    if not (Bd.OUT / (tag + '.json')).exists():
        hist = S.run('scalar', bg, Bd.X0, Bd.K0, Bd.SIG, Bd.T, Bd.DT, 10, save_final=Bd.RAW / (tag + '.npz'))
        S.save(tag, dict(tag=tag, bg='regge_iso', field='F1', h=Bd.H, dt=Bd.DT, sigma=Bd.SIG, T=Bd.T, hist=hist,
                         platform=S.platform_info()))
    rec = json.loads((Bd.OUT / (tag + '.json')).read_text()); t = np.array(rec['hist']['t'])
    ray = np.array(ray_iso(None))
    b = json.loads((Bd.OUT / 'bend.json').read_text())
    a = Bd.alpha(t, Bd.unwrap(rec['hist']['cc'])); ar = Bd.alpha(t, Bd.unwrap(ray))
    out = dict(alpha=a, alpha_ray=ar, ratio_field=a / b['lapse']['alpha'], ratio_ray=ar / b['lapse']['alpha_ray'],
               field_minus_ray_over_lapse_ray=(a - ar) / abs(b['lapse']['alpha_ray']))
    (Bd.OUT / 'posthoc_iso.json').write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps(out, indent=1))


if __name__ == '__main__':
    main()
