"""EX-GC38A: vacuum entanglement area law of the Y1 lattice scalar vs the 37A induced 1/(4G).
Plan: spec/EX-GC38A_plan.md. Our construction, not LaViolette's Model G.

A planar cut separates by transverse momentum into 1D chains V = 1/2 sum [a phi_u^2 - 2 b phi_u phi_{u+1}].
The ground-state entanglement depends only on r = a/b (rescaling V leaves X P invariant).
"""
import json, sys, time
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC38A"
Kf = lambda k: 2 * np.sin(k / 2)


def sym_entropy(X, P):
    nu = np.sqrt(np.clip(np.linalg.eigvals(X @ P).real, 0.25, None))
    a, b = nu + 0.5, nu - 0.5
    return float(np.sum(a * np.log(a) - np.where(b > 1e-14, b * np.log(np.where(b > 1e-14, b, 1)), 0)))


def chain_block_entropy(r, L, Np=1 << 16):
    """Infinite chain with E^2 = r - 2 cos p (b = 1). Entropy of an L-site block (two boundaries)."""
    p = 2 * np.pi * np.arange(Np) / Np
    E = np.sqrt(r - 2 * np.cos(p))
    x = np.real(np.fft.ifft(1 / (2 * E)))      # X(d) = (1/2pi) int cos(pd)/(2E)
    pp = np.real(np.fft.ifft(E / 2))
    d = np.abs(np.subtract.outer(np.arange(L), np.arange(L)))
    return sym_entropy(x[d], pp[d])


def chain_table(rs, L):
    return np.array([chain_block_entropy(r, L) / 2 for r in rs])   # per boundary


def area_law(m, orient, nt=1024, table=None):
    """Entanglement entropy per unit area of a single planar cut."""
    k = 2 * np.pi * (np.arange(nt) + 0.5) / nt - np.pi
    if orient == "z":
        KX, KY = np.meshgrid(k, k, indexing="ij")
        r = 2 + m * m + Kf(KX) ** 2 + Kf(KY) ** 2
        w, area = 1.0, 1.0
    else:   # (1,1,0): kappa in [-pi/2, pi/2), kz in [-pi, pi)
        kap = 0.5 * k
        KA, KZ = np.meshgrid(kap, k, indexing="ij")
        b = 2 * np.cos(KA)
        r = (4 + Kf(KZ) ** 2 + m * m) / b
        w, area = 1.0, np.sqrt(2)
    s = table(np.log(r - 2))
    return float(np.mean(s) * w / area)


def build_table(L, npts=400, lo=1e-4, hi=1e5):
    """Cubic spline of the per-boundary chain entropy vs ln(r - 2)."""
    from scipy.interpolate import CubicSpline
    lr = np.linspace(np.log(lo), np.log(hi), npts)
    return CubicSpline(lr, chain_table(2 + np.exp(lr), L))


def quarterG(n, m, pol):
    """1/(4G_ind) = 4 pi A, with A = dPi/dOmega^2 at (0, 0) done analytically: A = sum F^2(k, q=0) / (32 E^5).
    Same kernel as exgc37a.coeffs (which takes finite Omega steps); computed in kx slabs."""
    import exgc37a as A37
    k = 2 * np.pi * (np.arange(n) + 0.5) / n - np.pi
    ky, kz = np.meshgrid(k, k, indexing="ij")
    tot = 0.0
    for kx in k:
        E = A37.energy(kx, ky, kz, m)
        tot += np.sum(A37.F2(kx, ky, kz, 0.0, pol) / (32 * E ** 5))
    return 4 * np.pi * tot / n ** 3


def check_G2b(Lb=12, m=0.5):
    """(1,1,0) cut: brute force on a periodic Lb^3 box, region (x + y) mod Lb < Lb/2, vs the layer separation
    (layer n = x + y; transverse momenta q along the layer and kz; complex hopping 1 + e^{iq} between layers)."""
    k = 2 * np.pi * np.arange(Lb) / Lb
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    E = np.sqrt(m * m + Kf(KX) ** 2 + Kf(KY) ** 2 + Kf(KZ) ** 2)
    xc = np.real(np.fft.ifftn(1 / (2 * E))); pc = np.real(np.fft.ifftn(E / 2))
    S = np.array([(x, y, z) for x in range(Lb) for y in range(Lb) for z in range(Lb) if (x + y) % Lb < Lb // 2])
    D = (S[:, None, :] - S[None, :, :]) % Lb
    brute = sym_entropy(xc[D[..., 0], D[..., 1], D[..., 2]], pc[D[..., 0], D[..., 1], D[..., 2]])
    sep = 0.0
    nl = np.arange(Lb // 2)
    for q in k:
        t = 1 + np.exp(1j * q)
        for kz in k:
            a = 4 + Kf(kz) ** 2 + m * m
            # chain over layers n (period Lb): phi_n(x) ~ e^{i q x}; layer n couples to n+1 with amplitude t
            Ez = np.sqrt(a - 2 * np.abs(t) * np.cos(k + np.angle(t)))
            d = np.subtract.outer(nl, nl) % Lb
            X = (np.fft.ifft(1 / (2 * Ez)))[d]; P = (np.fft.ifft(Ez / 2))[d]
            sep += sym_entropy(X, P)
    return dict(Lb=Lb, m=m, brute=brute, separated=float(sep), rel=abs(brute - sep) / brute)


def check_G2(Lb=12, m=0.5):
    """Brute-force 3D periodic box Lb^3, region = slab z in [0, Lb/2): two planar cuts of area Lb^2 each."""
    k = 2 * np.pi * np.arange(Lb) / Lb
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    E = np.sqrt(m * m + Kf(KX) ** 2 + Kf(KY) ** 2 + Kf(KZ) ** 2)
    xc = np.real(np.fft.ifftn(1 / (2 * E)))
    pc = np.real(np.fft.ifftn(E / 2))
    sites = [(x, y, z) for x in range(Lb) for y in range(Lb) for z in range(Lb // 2)]
    S = np.array(sites)
    D = (S[:, None, :] - S[None, :, :]) % Lb
    X = xc[D[..., 0], D[..., 1], D[..., 2]]
    P = pc[D[..., 0], D[..., 1], D[..., 2]]
    brute = sym_entropy(X, P)
    # separated: transverse (kx, ky) on the Lb^2 grid, finite periodic chain of Lb, block Lb/2
    sep = 0.0
    for kx in k:
        for ky in k:
            r = 2 + m * m + Kf(kx) ** 2 + Kf(ky) ** 2
            Ez = np.sqrt(r - 2 * np.cos(k))
            x = np.real(np.fft.ifft(1 / (2 * Ez))); pp = np.real(np.fft.ifft(Ez / 2))
            d = np.abs(np.subtract.outer(np.arange(Lb // 2), np.arange(Lb // 2)))
            sep += sym_entropy(x[d], pp[d])
    return dict(Lb=Lb, m=m, brute=brute, separated=sep, rel=abs(brute - sep) / brute)


def check_G1():
    """Calabrese-Cardy: per-boundary entropy ~ (1/6) ln(xi), xi = 1/M, for M << 1."""
    Ms = np.array([0.004, 0.008, 0.016, 0.032])
    S = [chain_block_entropy(2 + M * M, 2000) / 2 for M in Ms]
    slope = np.polyfit(np.log(1 / Ms), S, 1)[0]
    return dict(M=Ms.tolist(), S=S, slope=float(slope), expected=1 / 6)


if __name__ == "__main__":
    t0 = time.process_time()
    part = sys.argv[1]
    res = {}
    if part == "checks":
        res["G1"] = check_G1()
        res["G2"] = check_G2()
        res["G2b_110"] = check_G2b()
        tab = build_table(800)
        rr = np.array([2.0013, 2.04, 2.7, 9.0, 300.0])
        direct = np.array([chain_block_entropy(r, 800) / 2 for r in rr])
        res["spline_check"] = dict(r=rr.tolist(), direct=direct.tolist(), spline=tab(np.log(rr - 2)).tolist())
        res["nt_convergence"] = {o: {nt: area_law(0.05, o, nt=nt, table=tab) for nt in (512, 1024, 2048)} for o in ("z", "110")}
        res["quarterG_n_convergence"] = {n: quarterG(n, 0.05, "cross") for n in (256, 512, 768)}
        # L convergence of the chain table at the smallest gap used
        res["L_convergence"] = {L: chain_block_entropy(2 + 0.05 ** 2, L) / 2 for L in (100, 200, 400, 800)}
    else:
        tab = build_table(800)
        ms = [0.02, 0.03, 0.05, 0.07, 0.1, 0.14, 0.2, 0.3, 0.5, 1.0]
        rows = []
        for m in ms:
            row = dict(m=m, S_z=area_law(m, "z", table=tab), S_110=area_law(m, "110", table=tab))
            for pol in ("cross", "plus34o"):
                n = 256 if m >= 0.1 else 768
                row[f"quarterG_{pol}"] = quarterG(n, m, pol)
            rows.append(row)
            print({k: round(v, 7) for k, v in row.items()}, flush=True)
        res["rows"] = rows
        # universal term: fit f(m) = a0 + a1 m^2 + a2 m^2 ln m + a3 m^4 on m <= 0.3
        sel = [r for r in rows if r["m"] <= 0.3]
        mm = np.array([r["m"] for r in sel])
        X = np.column_stack([np.ones_like(mm), mm ** 2, mm ** 2 * np.log(mm), mm ** 4])
        fits, variants = {}, {}
        sel2 = [r for r in rows if r["m"] <= 0.2]
        m2 = np.array([r["m"] for r in sel2])
        X2 = np.column_stack([np.ones_like(m2), m2 ** 2, m2 ** 2 * np.log(m2), m2 ** 4])
        X5 = np.column_stack([X, mm ** 4 * np.log(mm)])
        for key in ("S_z", "S_110", "quarterG_cross", "quarterG_plus34o"):
            y = np.array([r[key] for r in sel])
            c = np.linalg.lstsq(X, y, rcond=None)[0]
            fits[key] = dict(a0=float(c[0]), a1=float(c[1]), a2_m2lnm=float(c[2]), a3=float(c[3]),
                             max_resid=float(np.max(np.abs(X @ c - y))))
            y2 = np.array([r[key] for r in sel2])
            variants[key] = dict(m_le_0p2=float(np.linalg.lstsq(X2, y2, rcond=None)[0][2]),
                                 with_m4lnm=float(np.linalg.lstsq(X5, y, rcond=None)[0][2]))
        res["fits"] = fits          # primary (as written before the run): 4 terms, m <= 0.3
        res["fit_variants"] = variants
        # The plan stated -m^2 ln(m^2)/(24 pi), i.e. 1/(12 pi) on m^2 ln m: a factor-2 slip. The heat-kernel derivation
        # (1/(4G) = (Lambda^2 - m^2 ln(Lambda^2/m^2))/(48 pi)) and the area-law integral of -(1/12) ln(m^2 + k^2)
        # both give 1/(24 pi). Both values are recorded; the note scores against the registered one and says so.
        res["continuum_a2_m2lnm_registered"] = 1 / (12 * np.pi)
        res["continuum_a2_m2lnm_derived"] = 1 / (24 * np.pi)
        for r in rows:
            for o in ("S_z", "S_110"):
                for pol in ("cross", "plus34o"):
                    r[f"ratio_{o}_{pol}"] = r[o] / r[f"quarterG_{pol}"]
        print(json.dumps(fits, indent=1))
    res["cpu_seconds"] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=2) + "\n")
    print({k: v for k, v in res.items() if k not in ("rows",)})
