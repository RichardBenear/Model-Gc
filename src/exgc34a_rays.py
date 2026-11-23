"""EX-GC34A predictions by ray tracing (separate path from the field solver).

Hamiltonian from the exact lattice dispersion with local coefficients:
    H(x,k) = sqrt( [ sum_d w_d(x) 4 sin^2(h k.d/2)/h^2 + Vs(x) mu^2 ] / Ms(x) )
Rays sample the initial packet's Wigner function, a Gaussian in x (width sigma/sqrt2 per axis for |env|^2) and k
(width 1/(sigma sqrt2)); Gauss-Hermite nodes, weighted by energy ~ omega^2. The predicted centroid is the weighted
mean ray position. Only the background *formulas* are shared with the solver (re-implemented here, not imported).
"""
import numpy as np

DIRS = np.array([(1, 0, 0), (0, 1, 0), (0, 0, 1),
                 (1, 1, 0), (1, -1, 0), (1, 0, 1), (1, 0, -1), (0, 1, 1), (0, 1, -1)], float)


def weights_from_A(A, beta_div):
    """A: (...,3,3) -> (...,9) link weights; same null-space rule as registered."""
    w = np.zeros(A.shape[:-2] + (9,))
    pairs = {(0, 1): (3, 4), (0, 2): (5, 6), (1, 2): (7, 8)}
    beta = {}
    for (i, j), (p, m) in pairs.items():
        b = (A[..., i, i] + A[..., j, j]) / beta_div; beta[(i, j)] = b
        w[..., p] = b + A[..., i, j] / 2; w[..., m] = b - A[..., i, j] / 2
    for i in range(3):
        w[..., i] = A[..., i, i] - sum(2 * beta[q] for q in pairs if i in q)
    return w


class Bg:
    def __init__(self, kind, beta_div=12.0, gam=None, A=0.05, a=2.0, gppn=1.0, nordstrom=False, scale=1.0,
                 dt_lf=None):
        self.dt_lf = dt_lf
        self.kind = kind; self.bd = beta_div; self.gam = None if gam is None else np.array(gam, float)
        self.A = A; self.a = a; self.gppn = gppn; self.nord = nordstrom; self.scale = scale

    def coeffs(self, x):
        """x (...,3) -> Ms, Vs, w (...,9)."""
        if self.kind == 'bump':
            r = np.linalg.norm(x, axis=-1); Phi = -self.A / np.sqrt(r ** 2 + self.a ** 2)
            if self.nord:
                N = np.exp(Phi); B = Phi
            else:
                N = 1 + Phi; B = -self.gppn * Phi
            Ms = np.exp(3 * B) / N; Vs = N * np.exp(3 * B)
            Am = (N * np.exp(B))[..., None, None] * np.eye(3)
        elif self.kind == 'aniso':
            sg = np.sqrt(np.linalg.det(self.gam)); shp = x.shape[:-1]
            Ms = np.full(shp, sg); Vs = np.full(shp, sg)
            Am = np.broadcast_to(np.linalg.inv(self.gam) * sg, shp + (3, 3))
        else:
            shp = x.shape[:-1]; Ms = np.ones(shp); Vs = np.ones(shp); Am = np.broadcast_to(np.eye(3), shp + (3, 3))
        return Ms, Vs, weights_from_A(Am * self.scale, self.bd)

    def H(self, x, k, mu, h):
        Ms, Vs, w = self.coeffs(x)
        q = (w * 4 * np.sin(h * (k @ DIRS.T) / 2) ** 2).sum(-1) / h ** 2
        om = np.sqrt((q + Vs * mu ** 2) / Ms)
        if self.dt_lf:   # exact dispersion of the leapfrog time step: sin(W dt/2) = om dt/2
            om = 2 / self.dt_lf * np.arcsin(om * self.dt_lf / 2)
        return om


def gh_nodes(n):
    x, wt = np.polynomial.hermite_e.hermegauss(n)   # weight exp(-x^2/2)
    return x, wt / wt.sum()


def rays(bg, x0, k0, sigma, mu, h, T, dt, nq=5, times=None, eps=1e-6, box=None):
    """Energy-weighted mean ray position at the requested times."""
    sx = sigma / np.sqrt(2); sk = 1 / (sigma * np.sqrt(2))
    u, wu = gh_nodes(nq)
    grids = np.meshgrid(*([u] * 6), indexing='ij'); wg = np.meshgrid(*([wu] * 6), indexing='ij')
    U = np.stack([g.ravel() for g in grids], -1); W = np.prod(np.stack([g.ravel() for g in wg], -1), -1)
    X = np.asarray(x0, float) + sx * U[:, :3]; K = np.asarray(k0, float) + sk * U[:, 3:]
    om = bg.H(X, K, mu, h)
    wt = W * om ** 2; wt /= wt.sum()

    def rhs(X, K):
        dX = np.zeros_like(X); dK = np.zeros_like(K)
        for a in range(3):
            e = np.zeros(3); e[a] = eps
            dX[:, a] = (bg.H(X, K + e, mu, h) - bg.H(X, K - e, mu, h)) / (2 * eps)
            dK[:, a] = -(bg.H(X + e, K, mu, h) - bg.H(X - e, K, mu, h)) / (2 * eps)
        return dX, dK

    times = np.asarray(times); out = []
    t = 0.0; nstep = int(round(T / dt)); ti = 0
    for it in range(nstep + 1):
        while ti < len(times) and abs(times[ti] - t) < 0.5 * dt:
            if box is None:
                out.append((wt[:, None] * X).sum(0))
            else:   # 34A-r: periodic (circular) mean, same definition as the solver's readout
                lo, Lb = np.asarray(box[0], float), np.asarray(box[1], float)
                z = (wt[:, None] * np.exp(2j * np.pi * (X - lo) / Lb)).sum(0)
                out.append(lo + Lb / (2 * np.pi) * (np.angle(z) % (2 * np.pi)))
            ti += 1
        if it == nstep:
            break
        a1, b1 = rhs(X, K); a2, b2 = rhs(X + dt / 2 * a1, K + dt / 2 * b1)
        a3, b3 = rhs(X + dt / 2 * a2, K + dt / 2 * b2); a4, b4 = rhs(X + dt * a3, K + dt * b3)
        X = X + dt / 6 * (a1 + 2 * a2 + 2 * a3 + a4); K = K + dt / 6 * (b1 + 2 * b2 + 2 * b3 + b4)
        t = (it + 1) * dt
    return np.array(out)


def fit_v(t, c, t_from):
    t = np.asarray(t); m = t >= t_from
    return [float(np.polyfit(t[m], c[m, a], 1)[0]) for a in range(3)]
