"""EX-GC38B: lattice Bisognano-Wichmann and Tolman for the half-chain entanglement Hamiltonian.
Plan: spec/EX-GC38B_plan.md (A1 = implementation spec). Our construction, not LaViolette's Model G.

Open chain of 2L sites with Dirichlet walls, region = right half (sites 1..L), one cut between sites 0 and 1.
H = sum_j N_j (pi_j^2/2 + m^2 phi_j^2/2) + sum_links N_link (phi_{j+1} - phi_j)^2 / 2.
K = pi^T Nk pi / 2 + phi^T Mk phi / 2, computed in arbitrary precision (python-flint arb).
"""
import json, sys, time
from pathlib import Path
import numpy as np
import flint
from flint import arb, arb_mat, acb_mat

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC38B"


def funm(S, fs):
    """Matrix functions of a symmetric (diagonalisable, real-spectrum) arb_mat. fs: list of arb->arb maps."""
    n = S.nrows()
    lam, R = acb_mat(S).eig(right=True, algorithm="approx")
    V = arb_mat([[R[i, j].real for j in range(n)] for i in range(n)])
    lam = [l.real for l in lam]
    Vi = V.inv()
    outs = []
    for f in fs:
        D = arb_mat(n, n)
        for i in range(n):
            D[i, i] = f(lam[i])
        outs.append(V * D * Vi)
    return outs, lam


def ground_XP(L, m, lapse, mass_mode="physical"):
    """lapse: numpy array of N_j for the 2L sites (sites -L+1..L, index 0..2L-1)."""
    n = 2 * L
    m2 = arb(m) ** 2
    Ns = [arb(float(x)) for x in lapse]
    V = arb_mat(n, n)
    for j in range(n):
        # physical: N m^2 (H = sum N h). "conformal": m^2 / N, so that in optical length dy = dx/N the mass is uniform
        V[j, j] += (Ns[j] if mass_mode == "physical" else 1 / Ns[j]) * m2
    for j in range(-1, n):          # links incl. the walls (phi = 0 beyond the ends)
        Nl = Ns[j] if j == n - 1 else (Ns[j + 1] if j < 0 else (Ns[j] + Ns[j + 1]) / 2)
        for a in (j, j + 1):
            if 0 <= a < n:
                V[a, a] += Nl
        if 0 <= j and j + 1 < n:
            V[j, j + 1] -= Nl; V[j + 1, j] -= Nl
    Dh = arb_mat(n, n); Dmh = arb_mat(n, n)
    for j in range(n):
        Dh[j, j] = Ns[j].sqrt(); Dmh[j, j] = 1 / Ns[j].sqrt()
    S = Dh * V * Dh
    (Smh, Sh), _ = funm(S, [lambda l: 1 / l.sqrt(), lambda l: l.sqrt()])
    X = Dh * Smh * Dh / 2
    P = Dmh * Sh * Dmh / 2
    return X, P


def ent_hamiltonian(X, P):
    (Xh, Xmh), _ = funm(X, [lambda l: l.sqrt(), lambda l: 1 / l.sqrt()])
    A = Xh * P * Xh
    n = A.nrows()
    lam, R = acb_mat(A).eig(right=True, algorithm="approx")
    U = arb_mat([[R[i, j].real for j in range(n)] for i in range(n)])
    Ui = U.inv()
    nu = [l.real.sqrt() for l in lam]
    # modes beyond the working precision (2 nu - 1 below 2^-(prec/2), i.e. eps > ~prec ln2 / 2) are localised far
    # from the cut; their eps is capped there. Convergence in prec is checked (A1).
    floor = arb(2) ** (-(flint.ctx.prec // 2))
    eps = [((2 * v + 1) / max(2 * v - 1, floor)).log() for v in nu]
    D1 = arb_mat(n, n); D2 = arb_mat(n, n)
    for i in range(n):
        D1[i, i] = eps[i] / nu[i]; D2[i, i] = eps[i] * nu[i]
    Nk = Xh * U * D1 * Ui * Xh
    Mk = Xmh * U * D2 * Ui * Xmh
    return Nk, Mk, eps, nu


def tofloat(M, rows):
    return np.array([[float(M[i, j].mid()) for j in range(M.ncols())] for i in rows])


def run_case(L, m, eps_lapse=0.0, prec=500, jmax=25, mass_mode="physical"):
    flint.ctx.prec = prec
    x_site = np.arange(-L + 1, L + 1) - 0.5          # distance of site j from the cut (sites j = -L+1..L)
    f = np.where(x_site > 0, np.exp(-(x_site - 10) ** 2 / (2 * 16.0)), 0.0)
    lapse = 1 + eps_lapse * f
    t = time.process_time()
    X, P = ground_XP(L, m, lapse, mass_mode)
    idx = list(range(L, 2 * L))                       # region sites j = 1..L
    Xr = X.submatrix(L, L, L, L) if hasattr(X, "submatrix") else arb_mat([[X[i, k] for k in idx] for i in idx])
    Pr = arb_mat([[P[i, k] for k in idx] for i in idx])
    Xr = arb_mat([[X[i, k] for k in idx] for i in idx])
    Nk, Mk, eps, nu = ent_hamiltonian(Xr, Pr)
    rows = list(range(max(jmax, 21)))
    Nf, Mf = tofloat(Nk, rows), tofloat(Mk, rows)
    j = np.arange(1, jmax + 1)
    beta_pi = np.diag(Nf[:, :jmax]).copy()
    beta_link = -np.array([Mf[i, i + 1] for i in range(jmax - 1)])
    epsf = np.sort([float(e.mid()) for e in eps])
    # P3 non-locality in rows j <= 20
    nl_M = max(abs(Mf[i, k]) for i in range(20) for k in range(L) if abs(i - k) >= 2)
    nn_M = max(abs(Mf[i, i + 1]) for i in range(20))
    nl_N = max(abs(Nf[i, k]) for i in range(20) for k in range(L) if k != i)
    S_K = float(np.sum([e / np.expm1(e) - np.log(-np.expm1(-e)) for e in epsf if e < 700]))
    Xd = np.array([[float(Xr[a, b].mid()) for b in range(L)] for a in range(L)])
    Pd = np.array([[float(Pr[a, b].mid()) for b in range(L)] for a in range(L)])
    return dict(L=L, m=m, eps_lapse=eps_lapse, prec=prec, j=j.tolist(), x=(j - 0.5).tolist(),
                beta_pi=beta_pi.tolist(), beta_link=beta_link.tolist(), lapse_region=lapse[L:L + jmax].tolist(),
                f_region=f[L:L + jmax].tolist(), eps_low=epsf[:12].tolist(), S_from_K=S_K,
                nonlocal_M_over_nn=float(nl_M / nn_M), nonlocal_N_over_diag=float(nl_N / np.max(np.abs(beta_pi[:20]))),
                M_diag=np.diag(Mf[:, :jmax]).tolist(), cpu_s=time.process_time() - t), (Xd, Pd)


def sym_entropy(X, P):
    nu = np.sqrt(np.clip(np.linalg.eigvals(X @ P).real, 0.25, None))
    a, b = nu + 0.5, nu - 0.5
    return float(np.sum(a * np.log(a) - np.where(b > 1e-14, b * np.log(np.where(b > 1e-14, b, 1)), 0)))


if __name__ == "__main__":
    part = sys.argv[1]
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}
    if part == "test":
        r, _ = run_case(int(sys.argv[2]), float(sys.argv[3]), prec=int(sys.argv[4]))
        print(np.round(np.array(r["beta_pi"]) / (2 * np.pi * np.array(r["x"])), 5), r["cpu_s"])
        sys.exit()
    if part == "main":
        flat, (Xd, Pd) = run_case(100, 0.05, 0.0, 1000)
        flat800, _ = run_case(100, 0.05, 0.0, 1400)
        lap, _ = run_case(100, 0.05, 0.1, 1000)
        res["flat"] = flat; res["lapse"] = lap
        res["prec_check_max_rel"] = float(np.max(np.abs(np.array(flat800["beta_pi"]) / np.array(flat["beta_pi"]) - 1)))
        res["G1"] = dict(S_from_K=flat["S_from_K"], S_symplectic_double=sym_entropy(Xd, Pd))
        x = np.array(flat["x"]); bf = np.array(flat["beta_pi"]); bl = np.array(lap["beta_pi"])
        ratio = bf / (2 * np.pi * x)
        sel = (np.array(flat["j"]) >= 3) & (np.array(flat["j"]) <= 20)
        res["P1"] = dict(ratio=ratio.tolist(), min=float(ratio[sel].min()), max=float(ratio[sel].max()),
                         pass_=bool(np.all(np.abs(ratio[sel] - 1) <= 0.10)))
        # P2: predicted delta beta = 2 pi eps [x f(x) - int_0^x f], f on a fine grid
        xx = np.linspace(0, x.max(), 20001)
        ff = np.exp(-(xx - 10) ** 2 / 32.0)
        If = np.concatenate([[0], np.cumsum((ff[1:] + ff[:-1]) / 2 * np.diff(xx))])
        fx = np.exp(-(x - 10) ** 2 / 32.0)
        pred = 2 * np.pi * 0.1 * (x * fx - np.interp(x, xx, If))
        naive = 2 * np.pi * 0.1 * x * fx
        d = bl - bf
        rms = lambda v: float(np.sqrt(np.mean(v[sel] ** 2)))
        res["P2"] = dict(delta_beta=d.tolist(), predicted=pred.tolist(), naive_beta_prop_N=naive.tolist(),
                         rel_rms_err=rms(d - pred) / rms(pred), rel_rms_err_naive=rms(d - naive) / rms(naive),
                         pass_=bool(rms(d - pred) <= 0.10 * rms(pred)))
        print(json.dumps({k: v for k, v in res.items() if k in ("P1", "P2", "G1", "prec_check_max_rel")}, default=str)[:3000])
    if part == "g2":
        out = {}
        for m in (0.1, 0.05, 0.025):
            r, _ = run_case(int(round(5 / m)) if m > 0.03 else 200, m, 0.0, 1000, jmax=int(0.5 / m) + 2)
            x = np.array(r["x"]); b = np.array(r["beta_pi"])
            out[m] = dict(x_at_xm_half=0.5 / m, ratio=float(np.interp(0.5 / m, x, b / (2 * np.pi * x))), cpu=r["cpu_s"])
            print(m, out[m], flush=True)
        res["G2"] = out
    if part == "p2diag":
        # descriptive, after the registered P2 failed: is the miss the m^2 N^2 term? (a) smaller m, (b) "conformal" mass
        x = np.arange(1, 26) - 0.5
        xx = np.linspace(0, x.max(), 20001); ff = np.exp(-(xx - 10) ** 2 / 32.0)
        If = np.concatenate([[0], np.cumsum((ff[1:] + ff[:-1]) / 2 * np.diff(xx))])
        pred = 2 * np.pi * 0.1 * (x * np.exp(-(x - 10) ** 2 / 32.0) - np.interp(x, xx, If))
        sel = (x >= 2.5) & (x <= 19.5)
        rms = lambda v: float(np.sqrt(np.mean(v[sel] ** 2)))
        out = {}
        for m, mode in ((0.05, "conformal"), (0.02, "physical"), (0.005, "physical"), (0.005, "conformal")):
            f0, _ = run_case(100, m, 0.0, 1000, mass_mode=mode)
            f1, _ = run_case(100, m, 0.1, 1000, mass_mode=mode)
            d = np.array(f1["beta_pi"]) - np.array(f0["beta_pi"])
            out[f"m{m}_{mode}"] = dict(delta_beta=d.tolist(), rel_rms_err=rms(d - pred) / rms(pred),
                                       flat_ratio_j20=f0["beta_pi"][19] / (2 * np.pi * 19.5))
            print(m, mode, out[f"m{m}_{mode}"]["rel_rms_err"], out[f"m{m}_{mode}"]["flat_ratio_j20"], flush=True)
        res["predicted"] = pred.tolist(); res["cases"] = out
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=2, default=float) + "\n")
