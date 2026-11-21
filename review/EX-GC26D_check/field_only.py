"""EX-GC26D post-run diagnostic (Claude): is the fore/aft asymmetry at d >= 3 a taper artifact?

Field-only evolution Phi_tt = lap Phi - 4 pi G_N rho with a rigidly moving, Lorentz-contracted source
rho(x, t) = rho_rest(gamma (x - x0 - u t), y, z) (no back-reaction; the knot's own response is ~1e-4).
  case 'taper': the PC initial data (exgc26d.rest_field_fn + boost_field), Phi = 0 on the wall  -> should reproduce the PC co-motion
  case 'exact': no far taper; wall held at the exactly boosted -A/r at every step                -> clean co-motion
Inputs: local_raw/EX-GC26D (hash-checked as in exgc26d).  python review/EX-GC26D_check/field_only.py 0.3
"""
import sys, json, hashlib, time, platform
from pathlib import Path
import numpy as np, scipy
from scipy.ndimage import map_coordinates, spline_filter
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'src'))
import exgc26c as c, exgc26e1 as e1, exgc26d as d

u = float(sys.argv[1]); GN = 0.01; T = d.T_RUN; x0 = -u * T / 2; g = 1 / np.sqrt(1 - u * u)
e1.configure(0.2); e1.load(); G = e1.ref.Grid(e1.H, e1.L); Fd = c.Field(G, c.L_PHI)
man = json.loads((d.OUT / 'inputs.json').read_text())
def inp(n):
    b = (d.RAW / n).read_bytes(); assert hashlib.sha256(b).hexdigest() == man[n]; return np.load(d.RAW / n)
p = inp('c26_p.npy'); Fc = inp('c26_Fcrop.npy'); k = man['crop_half_cells']; cc = Fd.n // 2
Frest = np.zeros(Fd.r.shape); Frest[cc - k:cc + k + 1, cc - k:cc + k + 1, cc - k:cc + k + 1] = Fc
A = d.fit_A(Fd, Frest)
rho_rest = c.rho(G, p, np.zeros_like(p), np.exp(2 * Frest[Fd.sub]))
src_coef = spline_filter(rho_rest, order=3, mode='constant')
Xm = np.array(np.meshgrid(G.ax, G.ax, G.ax, indexing='ij'))
def source(t):   # on the matter sub-grid
    idx = Xm.copy(); idx[0] = g * (Xm[0] - x0 - u * t)
    return map_coordinates(src_coef, idx / G.h + (G.n - 1) / 2, order=3, mode='constant', prefilter=False)
Xf = np.array(np.meshgrid(Fd.ax, Fd.ax, Fd.ax, indexing='ij'))
wall = np.zeros(Fd.r.shape, bool); wall[[0, -1]] = wall[:, [0, -1]] = wall[:, :, [0, -1]] = True
Xw = Xf[:, wall]
def wall_val(t):
    xr = g * (Xw[0] - x0 - u * t); return -A / np.sqrt(xr ** 2 + Xw[1] ** 2 + Xw[2] ** 2)

def run(case):
    t0 = time.process_time()
    if case == 'taper':
        fn = d.rest_field_fn(Fd, Frest, A)
    else:
        coeff = spline_filter(Frest, order=3, mode='constant')
        def fn(X):
            r = np.sqrt(np.sum(X ** 2, axis=0)); inner = map_coordinates(coeff, X / Fd.h + (Fd.n - 1) / 2, order=3, mode='nearest', prefilter=False)
            w = np.clip((10.5 - r) / 1.0, 0, 1); w = w * w * (3 - 2 * w); return w * inner + (1 - w) * (-A / np.maximum(r, 1e-9))
    F, Pi = d.boost_field(Fd, fn, u, x0)
    if case == 'exact':
        F[wall] = wall_val(0.); Pi[wall] = 0.
    def deriv(F, Pi, t):
        PiD = Fd.lap(F); PiD[Fd.sub] -= 4 * np.pi * GN * source(t); PiD[wall] = 0.
        FD = Pi.copy()
        if case == 'exact':
            FD[wall] = (wall_val(t + 1e-4) - wall_val(t - 1e-4)) / 2e-4   # wall follows the exact moving field
        else:
            FD[wall] = 0.
        return FD, PiD
    dt = .05 * G.h; steps = round(T / dt)
    for i in range(steps):
        t = i * dt
        k1 = deriv(F, Pi, t); k2 = deriv(F + dt / 2 * k1[0], Pi + dt / 2 * k1[1], t + dt / 2)
        k3 = deriv(F + dt / 2 * k2[0], Pi + dt / 2 * k2[1], t + dt / 2); k4 = deriv(F + dt * k3[0], Pi + dt * k3[1], t + dt)
        F = F + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]); Pi = Pi + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
    X = x0 + u * T; coeff = spline_filter(F, order=3, mode='constant')
    at = lambda x, y, z: float(map_coordinates(coeff, np.array([[x], [y], [z]]) / Fd.h + (Fd.n - 1) / 2, order=3, prefilter=False)[0])
    out = []
    for dd in (2., 3., 4., 5.):
        fo, af, si = at(X + dd, 0, 0), at(X - dd, 0, 0), at(X, dd, 0)
        out.append(dict(d=dd, fore=fo, aft=af, side=si, fore_over_aft=fo / af, axial_over_side=(fo + af) / 2 / si))
    return dict(case=case, comotion=out, cpu_s=time.process_time() - t0)

res = dict(u=u, GN=GN, A=A, gamma_inv=1 / g, runs=[run('taper'), run('exact')],
           platform=f'Python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} {platform.platform()}')
Path(__file__).with_name(f'field_only_u{u:g}.json').write_text(json.dumps(res, indent=1) + '\n'); print(json.dumps(res, indent=1))
