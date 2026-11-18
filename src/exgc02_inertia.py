"""EX-GC02 (exploratory): inertia in Model Gc via species response times (drift bifurcation).

tau_s d_t u_s = F_s(u). Steady particle unchanged. l = 1 dipole problem: J1 psi = lambda T psi.
Drag gamma(tau) = sum_s tau_s <psi+_s, U'_s>, mobility mu = -<psi+, P>/gamma, velocity relaxation rate lambda_v,
effective mass m = gamma / |lambda_v|.
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc01_unequal as g1  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC02")
NAMES = ("X", "Y", "W", "G")


def load(name="ref"):
    p = g1.TYPES[name]
    z = np.load(os.path.join(e.ROOT, "results", "EX-GC01", f"radial_{name}.npz"))
    M = g1.RadialP(120.0, 0.05, p["k"], p["c"], p["dw"])
    u, res, ok = M.newton(z["u"])
    assert ok
    return M, u


def Jl(M, u, l):
    """Jacobian for angular index l (l = 0 Neumann-type as model; l >= 1 Dirichlet-zero outer plus -l(l+1)/r^2)."""
    if l <= 1:
        return M.J(u, l)
    J1 = M.J(u, 1)
    extra = sp.diags(np.tile(-(l * (l + 1) - 2) / M.r ** 2, 4) * np.repeat([1.0, e.DY, M.dw, 1.0], M.n))
    return (J1 + extra).tocsc()


def null_vectors(M, u):
    """Left/right null vectors of J1 (translation), right scaled to U'."""
    J1 = Jl(M, u, 1)
    Up = np.array([np.gradient(u[s], M.r) for s in range(4)]).ravel()
    lu, luT = spla.splu(J1), spla.splu(J1.T.tocsc())
    rv, lv = Up.copy(), np.random.default_rng(0).standard_normal(J1.shape[0])
    for _ in range(10):
        rv = lu.solve(rv); rv /= np.linalg.norm(rv)
        lv = luT.solve(lv); lv /= np.linalg.norm(lv)
    rv *= (Up @ Up) / (rv @ Up)
    return J1, lv, rv, Up


def forcing(M, u):
    X, Y, W, G = u
    P = np.zeros((4, M.n))
    P[0] = M.cx * M.k * W * M.r
    P[1] = M.cy * M.k * W * M.r
    P[3] = -M.k * W * M.r
    return P.ravel()


def weights(M, u):
    J1, lv, rv, Up = null_vectors(M, u)
    n = M.n
    c = np.array([lv[s * n:(s + 1) * n] @ rv[s * n:(s + 1) * n] for s in range(4)])
    p = float(lv @ forcing(M, u))
    s = np.sign(c.sum())                   # normalise so gamma(tau = 1) = +1
    c, p = c * s / abs(c.sum()), p * s / abs(c.sum())
    return c, p                            # gamma(tau) = c . tau ; mu = -p / gamma


def gen_eigs(J, Tdiag, sigma, k=8):
    """Eigenvalues of J psi = lambda T psi near sigma (T diagonal positive)."""
    A = (sp.diags(1 / Tdiag) @ J).tocsc()
    vals, vecs = spla.eigs(A, k=k, sigma=sigma, which="LM")
    return vals, vecs


def tau_vec(M, tau):
    return np.repeat([tau.get(s, 1.0) for s in NAMES], M.n)


def velocity_mode(M, u, tau, k=10):
    """Slow l = 1 modes; the velocity mode is the core-localised one (largest W share) other than translation."""
    J1 = Jl(M, u, 1)
    vals, vecs = gen_eigs(J1, tau_vec(M, tau), sigma=1e-6, k=k)
    n = M.n
    out = []
    Up = np.array([np.gradient(u[s], M.r) for s in range(4)]).ravel()
    for lam, v in zip(vals, vecs.T):
        v = v / np.linalg.norm(v)
        core = float(np.sum(np.abs(v[2 * n:3 * n]) ** 2))         # W share (core-only species)
        ovl = float(abs(np.vdot(v, Up)) / np.linalg.norm(Up))
        out.append(dict(re=float(lam.real), im=float(lam.imag), W_share=core, overlap_Up=ovl))
    return sorted(out, key=lambda d: -d["re"])


def max_growth(M, u, tau, ls=(0, 1, 2), sigmas=(0.0, 0.3j, 1j, 2j, 4j)):
    """Largest real part found near several shifts (excluding the l = 1 translation mode)."""
    best = {}
    T = tau_vec(M, tau)
    for l in ls:
        J = Jl(M, u, l)
        m = -np.inf
        for sg in sigmas:
            try:
                vals, _ = gen_eigs(J, T, sigma=sg + 1e-6, k=6)
            except Exception:
                continue
            for lam in vals:
                if l == 1 and abs(lam) < 1e-3 and abs(lam.imag) < 1e-6 and lam.real > -5e-4:
                    continue            # translation (finite-domain shifted) mode
                m = max(m, float(lam.real))
        best[l] = m
    return best


if __name__ == "__main__" and sys.argv[1:2] == ["weights"]:
    os.makedirs(OUT, exist_ok=True)
    for name in ("ref", "dw0.2", "k0.4"):
        M, u = load(name)
        c, p = weights(M, u)
        print(name, "weights X,Y,W,G:", np.round(c, 5), "P:", p, "mu(tau=1):", -p, "drift():", M.drift(u)["mu"])


def mass(M, u, tau):
    """Jordan-chain effective mass: J1 w = -(T U' - alpha U'), <psi+, w> = 0, m = -<psi+, T w> (normalised like gamma)."""
    J1, lv, rv, Up = null_vectors(M, u)
    n = M.n
    T = tau_vec(M, tau)
    c = np.array([lv[s * n:(s + 1) * n] @ rv[s * n:(s + 1) * n] for s in range(4)])
    norm = np.sign(c.sum()) / abs(c.sum())
    gam = float(lv @ (T * rv))
    alpha = gam / float(lv @ rv)
    rhs = -(T * rv - alpha * rv)
    B = sp.bmat([[J1, sp.csr_matrix(rv[:, None])], [sp.csr_matrix(lv[None, :]), None]]).tocsc()
    sol = spla.spsolve(B, np.concatenate([rhs, [0.0]]))
    w = sol[:-1]
    m = -float(lv @ (T * w))
    parts = [-float(lv[s * n:(s + 1) * n] @ (T * w)[s * n:(s + 1) * n]) * norm for s in range(4)]
    return dict(gamma=gam * norm, m=m * norm, m_parts=parts, mu=-float(lv @ forcing(M, u)) / gam, relax_time=m / gam)


def time_response(M, u, tau, E=1e-3, dt=1.0, T_end=1500.0):
    """Linear l = 1 response: T dx/dt = J1 x + E P, x(0) = 0. Position from the W dipole, normalised to a shift."""
    J1 = Jl(M, u, 1)
    Tm = tau_vec(M, tau)
    P = forcing(M, u)
    n = M.n
    W0p = np.gradient(u[2], M.r)
    norm = float(np.sum(M.V * M.r * (-W0p)))
    A = (sp.diags(Tm / dt) - J1).tocsc()
    lu = spla.splu(A)
    x = np.zeros(4 * n)
    ts, Xs = [0.0], [0.0]
    for i in range(int(T_end / dt)):
        x = lu.solve(Tm / dt * x + E * P)
        ts.append((i + 1) * dt)
        Xs.append(float(np.sum(M.V * M.r * x[2 * n:3 * n])) / norm)
    return np.array(ts), np.array(Xs)


def ode_prediction(ts, gamma, m, ks, f):
    """m X'' + gamma X' + ks X = f, X(0) = X'(0) = 0 (all in the normalised units of gamma)."""
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, y: [y[1], (f - gamma * y[1] - ks * y[0]) / m] if m > 0 else [(f - ks * y[0]) / gamma, 0],
                    (0, ts[-1]), [0.0, 0.0], t_eval=ts, rtol=1e-9, atol=1e-14)
    return sol.y[0]


def run():
    os.makedirs(OUT, exist_ok=True)
    out = dict(weights={}, scan={}, stability={}, response={})
    variants = {}
    for name in ("ref", "dw0.2", "k0.4"):
        M, u = load(name)
        variants[name] = (M, u)
        c, p = weights(M, u)
        out["weights"][name] = dict(c=c.tolist(), p=p, tauW_threshold=1 + c.sum() / -c[2] if c[2] < 0 else None)
        print("weights", name, out["weights"][name], flush=True)
    taus = [1, 5, 10, 15, 20, 22.5, 25, 26, 27]
    for name, (M, u) in variants.items():
        out["scan"][name] = []
        for tw in taus:
            d = mass(M, u, {"W": tw}); d["tauW"] = tw
            out["scan"][name].append(d)
        print("scan", name, [(r["tauW"], round(r["gamma"], 4), round(r["m"], 3), round(r["mu"], 5)) for r in out["scan"][name]], flush=True)
    M, u = variants["ref"]
    for tw in (1, 10, 20, 25, 27):
        out["stability"][str(tw)] = max_growth(M, u, {"W": tw})
        print("stability tauW", tw, out["stability"][str(tw)], flush=True)
    # direct response on a finer grid (weaker grid pinning)
    Mf = g1.RadialP(120.0, 0.025, 0.2, 1.0, 0.1)
    uf, res, ok = Mf.newton(np.array([np.interp(Mf.r, M.r, u[s]) for s in range(4)]))
    ks = -velocity_mode(Mf, uf, {}, k=4)[0]["re"]       # pinning rate at tau = 1 (gamma = 1 units)
    for tw in (1, 20, 25):
        ts, Xs = time_response(Mf, uf, {"W": tw})
        d = mass(Mf, uf, {"W": tw})
        f = d["mu"] * 1e-3 * d["gamma"]
        pred_m = ode_prediction(ts, d["gamma"], d["m"], ks, f)
        pred_0 = ode_prediction(ts, d["gamma"], 0.0, ks, f)
        out["response"][str(tw)] = dict(t=ts[::10].tolist(), X=Xs[::10].tolist(), pred_with_mass=pred_m[::10].tolist(),
                                       pred_no_mass=pred_0[::10].tolist(), gamma=d["gamma"], m=d["m"], ks=ks, mu=d["mu"])
        err_m = float(np.max(np.abs(Xs - pred_m)) / np.max(np.abs(Xs)))
        err_0 = float(np.max(np.abs(Xs - pred_0)) / np.max(np.abs(Xs)))
        out["response"][str(tw)].update(rel_err_with_mass=err_m, rel_err_no_mass=err_0)
        print("response tauW", tw, "X(end)", Xs[-1], "err with mass", err_m, "without", err_0, flush=True)
    json.dump(out, open(os.path.join(OUT, "results.json"), "w"), indent=1, default=float)


if __name__ == "__main__" and sys.argv[1:2] == ["run"]:
    run()
