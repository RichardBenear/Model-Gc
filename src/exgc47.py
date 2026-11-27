"""EX-GC47: two bound heavy pairs (45/46 model) coupled through the prescribed 44 fields.
47A benchmark pairs; 47B static cross energy from microscopic charge profiles (pair form factors); 47C release
(COM effective dynamics, validated against the exact two-body problem at L = 6).
Plan: spec/EX-GC47_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
Conventions: k = 2 pi n / L; o(q) = sum_v e^{+iq.v} o_v maps k -> k+q; a pair state at total momentum K is
|Psi_K> = sum_k phi_K(k) |1:k+K>|2:-k>, phi_K ~ A_K/(D_K - E_K), A_K = conj(u6(k+K)) conj(u6(-k)) (natural gauge:
the band-projected resolvent of the on-site heavy pair Bloch sum; gauge invariant in the Bloch phases).
"""
import json, sys, time, itertools
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import exgc45 as X

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC47"
LAM0, C2CHI = 0.1, 1.001
BENCH = {"B1": (0.5, 0.025, 0.15), "B2": (0.7, 0.025, 0.15)}
T46 = {"B1": 5.8805199427528924e-05, "B2": 4.580653348367047e-05}   # 46B one-step t_pair (four bands)
neg = lambda a: np.roll(np.flip(a, axis=(0, 1, 2)), 1, axis=(0, 1, 2))
sh = lambda a, n: np.roll(a, [-x for x in n], axis=(0, 1, 2))          # f(k + q) for q = 2 pi n / L


def setup(L, m0, g, frac=None, U=None):
    e, u, K = X.eig(L, LAM0, m0=m0, g=g)
    eh, uh = e[..., 3], u[..., :, 3]
    W = float(eh.max() - eh.min())
    if U is None:
        U = X.U_for(eh, uh, frac * W)
    return dict(L=L, m0=m0, g=g, eh=eh, uh=uh, W=W, U=U, K=K)


def bound(S, n=(0, 0, 0)):
    L = S["L"]; N = L ** 3; e, u6 = S["eh"], S["uh"][..., 6]
    e1, u1, e2, u2 = sh(e, n), sh(u6, n), neg(e), neg(u6)
    D = e1 + e2; A = np.conj(u1) * np.conj(u2); A2 = np.abs(A) ** 2; Dmin = D.min()
    f = lambda E: S["U"] * np.sum(A2 / (D - E)) / N - 1
    if f(Dmin - 1e-14) < 0:
        return None
    E0 = brentq(f, Dmin - 50, Dmin - 1e-14, xtol=1e-15, rtol=1e-14)
    phi = A / (D - E0); phi /= np.linalg.norm(phi)
    a = np.sum(phi * u1 * u2)
    return dict(phi=phi, E0=float(E0), Dmin=float(Dmin), a=a, V=float(-S["U"] / N * abs(a) ** 2))


def Bu(S, uv, nq, strain=False):
    """B(k, q) u(k): lapse (site + half-weighted bond) or strain (bond only) operator applied to Bloch vectors."""
    L = S["L"]; K = S["K"]; q = 2 * np.pi * np.array(nq) / L; out = np.zeros_like(uv)
    if not strain:
        Hs = LAM0 * (np.ones((7, 7)) - np.eye(7)); Hs[6, :] = S["g"]; Hs[:, 6] = S["g"]; Hs[6, 6] = S["m0"]
        out += uv @ Hs.T
    for a in range(3):
        ph = np.exp(1j * K[..., a])
        out[..., 2 * a] += -ph * (1 + np.exp(1j * q[a])) / 2 * uv[..., 2 * a + 1]
        out[..., 2 * a + 1] += -np.conj(ph) * (1 + np.exp(-1j * q[a])) / 2 * uv[..., 2 * a]
    return out


def M(S, nq, strain=False):
    """single-carrier matrix element <k+q| o(q) |k> on the heavy band, as an array over k."""
    return np.sum(np.conj(sh(S["uh"], nq)) * Bu(S, S["uh"], nq, strain), axis=-1)


def formfactor(S, nq, b0, bK=None, nK=(0, 0, 0)):
    """F(q; K) = <Psi_{K+q}| o(q) |Psi_K>, lapse (N) and strain (a); bK is the bound state at K+q."""
    nKq = tuple(np.add(nK, nq) % S["L"])
    if bK is None:
        bK = bound(S, nKq)
    if bK is None:
        return None
    out = {}
    for lab, st in (("N", False), ("a", True)):
        Mq = M(S, nq, st)
        # species 1: |1:k+K> -> |1:k+K+q>, k unchanged;  species 2: |2:-k> -> |2:-k+q>, k -> k-q
        F1 = np.sum(np.conj(bK["phi"]) * b0["phi"] * sh(Mq, nK))
        F2 = np.sum(np.conj(np.roll(bK["phi"], list(nq), axis=(0, 1, 2))) * b0["phi"] * neg(Mq))
        out[lab] = F1 + F2
    out["V"] = -S["U"] / S["L"] ** 3 * np.conj(bK["a"]) * b0["a"]
    out["N"] = out["N"] + out["V"]
    return out


def wedge(L):
    h = L // 2
    return [(a, b, c) for a in range(h + 1) for b in range(a, h + 1) for c in range(b, h + 1)]


def fill(L, table):
    n = np.arange(L); m = np.minimum(n, L - n)
    G = np.stack(np.meshgrid(m, m, m, indexing="ij"), -1); G.sort(axis=-1)
    out = np.zeros((L, L, L), complex)
    for idx in itertools.product(range(L), repeat=3):
        out[idx] = table[tuple(G[idx])]
    return out


def pair_tables(S):
    L = S["L"]; b0 = bound(S); FN, Fa, Eb = {}, {}, {}
    for w in wedge(L):
        bq = bound(S, w); Eb[w] = (bq["Dmin"] - bq["E0"], bq["E0"]) if bq else None
        F = formfactor(S, w, b0, bq)
        FN[w], Fa[w] = (F["N"], F["a"]) if F else (np.nan, np.nan)
    return b0, fill(L, FN), fill(L, Fa), Eb


def kernel(L, G4pi):
    n = 2 * np.pi * np.fft.fftfreq(L) ; Q = np.meshgrid(n, n, n, indexing="ij")
    qt2 = sum(2 * (1 - np.cos(x)) for x in Q); Kq = np.zeros_like(qt2)
    Kq[qt2 > 0] = -G4pi / ((1 + G4pi * C2CHI) * qt2[qt2 > 0])
    return Kq


def Vtilde(FS, FP, Kq):
    """Ṽ'(q) for point source at 0 and point probe: probe energy (1/N) sum_q Ṽ'(q) e^{iqR}."""
    return Kq * ((FS["N"] + FS["a"]) * FP["N"] + FS["N"] * FP["a"])


def envelope(L, R, sigma, sign=+1):
    n = np.arange(L); d = [((n - r + L // 2) % L) - L // 2 for r in R]
    if sigma == 0:
        A = np.zeros((L, L, L)); A[tuple(np.array(R) % L)] = 1.0
    else:
        D = np.meshgrid(*d, indexing="ij"); A = np.exp(-sum(x ** 2 for x in D) / (2 * sigma ** 2)); A /= A.sum()
    return L ** 3 * np.fft.ifftn(A) if sign > 0 else np.fft.fftn(A)        # sum_R |c|^2 e^{+-iqR}


def V12(FS, FP, Kq, R, sS=0, sP=0):
    L = Kq.shape[0]
    AS = envelope(L, (0, 0, 0), sS, -1); AP = envelope(L, R, sP, +1)
    return np.sum(Vtilde(FS, FP, Kq) * AS * AP) / L ** 3


# ---------------------------------------------------------------- 47A
def part_A(L=32, h=1e-5):
    res = {}
    for name, (m0, g, fr) in BENCH.items():
        S = setup(L, m0, g, fr); b0 = bound(S)
        dirs = {"x": (1, 0, 0), "y": (0, 1, 0), "z": (0, 0, 1), "xyz": (1, 1, 1)}; disp = {}
        for dn, v in dirs.items():
            rows = []
            for n in range(0, 9):
                b = bound(S, tuple(n * np.array(v))); kk = (2 * np.pi * n / L) ** 2 * sum(v)
                rows.append((n, kk, b["E0"] - b0["E0"], b["Dmin"] - b["E0"]))
            r = np.array(rows); sel = r[1:5]
            c = np.linalg.lstsq(np.stack([sel[:, 1], sel[:, 1] ** 2], 1), sel[:, 2], rcond=None)[0]
            disp[dn] = dict(rows=r.tolist(), t_fit=float(c[0]), u_fit=float(c[1]), t_onestep_n2=float(r[2, 2] / r[2, 1]),
                            Eb_min_over_K=float(r[:, 3].min()))
        # charges
        rp = bound(dict(S, **{k: v for k, v in zip(("eh", "uh"), (lambda e, u: (e[..., 3], u[..., :, 3]))(*X.eig(L, LAM0, m0=m0, g=g, t=(1 + h,) * 3)[:2]))}))
        rm = bound(dict(S, **{k: v for k, v in zip(("eh", "uh"), (lambda e, u: (e[..., 3], u[..., :, 3]))(*X.eig(L, LAM0, m0=m0, g=g, t=(1 - h,) * 3)[:2]))}))
        Qa = (rp["E0"] - rm["E0"]) / (2 * h)
        # pair size: P(r) = N sum_ij |(1/N) sum_k phi(k) u_i(k) u_j(-k) e^{ikr}|^2
        uh = S["uh"]; un = neg(uh); P = np.zeros((L, L, L))
        for i in range(7):
            for j in range(7):
                P += np.abs(L ** 3 * np.fft.ifftn(b0["phi"] * uh[..., i] * un[..., j])) ** 2 / L ** 3
        n = np.arange(L); d = (n + L // 2) % L - L // 2; Dg = np.meshgrid(d, d, d, indexing="ij"); r2 = sum(x ** 2 for x in Dg)
        t_iso = disp["x"]["t_fit"]
        res[name] = dict(m0=m0, g=g, frac=fr, U=S["U"], W=S["W"], E=b0["E0"], Eb=b0["Dmin"] - b0["E0"], V=b0["V"], Qa=Qa,
                         Qa_over_E=Qa / b0["E0"], P_norm=float(P.sum()), rms_sep=float(np.sqrt((P * r2).sum() / P.sum())),
                         P_onsite=float(P[0, 0, 0] / P.sum()), disp=disp, t_pair=t_iso, t46=T46[name],
                         t_vs_46=t_iso / T46[name] - 1, t_spread=float(np.ptp([disp[k]["t_fit"] for k in disp]) / t_iso),
                         inertia_over_E=1 / (2 * t_iso * b0["E0"]))
        print("47A", name, {k: v for k, v in res[name].items() if k != "disp"}, flush=True)
    return res


# ---------------------------------------------------------------- 47B
def carrier_tables(S):
    """one heavy carrier at the band minimum, gauge u6 real positive (natural single-carrier gauge)."""
    L = S["L"]; uh = S["uh"] * (np.conj(S["uh"][..., 6:7]) / np.abs(S["uh"][..., 6:7])); S1 = dict(S, uh=uh)
    i = np.unravel_index(np.argmin(S["eh"]), S["eh"].shape); FN, Fa = {}, {}
    for w in wedge(L):
        for lab, st, T in (("N", False, FN), ("a", True, Fa)):
            T[w] = M(S1, w, st)[i]
    return dict(N=fill(L, FN), a=fill(L, Fa)), float(S["eh"][i])


def part_B(Ls=(16, 24, 32), G4pis=(0.005, 0.01, 0.02)):
    A = json.loads((OUT / "A.json").read_text())["result"]; res = {}
    for L in Ls:
        T = {}
        for name, (m0, g, fr) in BENCH.items():
            S = setup(L, m0, g, U=A[name]["U"]); b0, FN, Fa, Eb = pair_tables(S)
            T[name] = dict(N=FN, a=Fa, E=b0["E0"], Qa=float(FN.real.sum() * 0 + Fa[0, 0, 0].real), E_F0=float(FN[0, 0, 0].real),
                           Fimag=float(max(np.abs(FN.imag).max(), np.abs(Fa.imag).max())), unbound=int(sum(v is None for v in Eb.values())))
        S = setup(L, *BENCH["B1"][:2], U=A["B1"]["U"]); Fc, ec = carrier_tables(S)
        T["C1"] = dict(N=Fc["N"], a=Fc["a"], E=ec, Qa=float(Fc["a"][0, 0, 0].real), E_F0=float(Fc["N"][0, 0, 0].real))
        rows = {}
        for G4pi in G4pis:
            Kq = kernel(L, G4pi); Geff = G4pi / (4 * np.pi * (1 + G4pi * C2CHI))
            for (s, p) in (("B1", "B1"), ("B1", "B2"), ("B2", "B1"), ("C1", "C1")):
                FS, FP = T[s], T[p]
                pS = dict(N=np.full_like(FS["N"], FS["E_F0"]), a=np.full_like(FS["a"], FS["Qa"]))
                pP = dict(N=np.full_like(FP["N"], FP["E_F0"]), a=np.full_like(FP["a"], FP["Qa"]))
                dS = dict(N=pS["N"], a=0 * pS["a"]); dP = dict(N=pP["N"], a=0 * pP["a"])
                for R in [(r, 0, 0) for r in range(2, L // 2 + 1)] + [(r, r, 0) for r in range(2, L // 4 + 1)] + [(r, r, r) for r in range(1, L // 4 + 1)]:
                    for sg in ((0, 0), (2, 2)):
                        v = V12(FS, FP, Kq, R, *sg); vp = V12(pS, pP, Kq, R, *sg); vd = V12(dS, dP, Kq, R, *sg)
                        r = float(np.sqrt(sum(x * x for x in R)))
                        bare = -Geff * (FS["E_F0"] * FP["E_F0"] + FS["E_F0"] * FP["Qa"] + FS["Qa"] * FP["E_F0"]) / r
                        rows[f"{G4pi}|{s}>{p}|{R}|{sg}"] = dict(V=v.real, Vimag=v.imag, Vpoint=vp.real, Vdust=vd.real, Vbare=bare,
                                                               ratio_point=v.real / vp.real, ratio_dust=v.real / vd.real, r=r, VoverGeff=v.real / Geff)
            # G = 0 control
        rows["G0"] = float(V12(T["B1"], T["B1"], kernel(L, 0.0), (6, 0, 0)).real)
        res[str(L)] = dict(tables={k: {kk: vv for kk, vv in v.items() if kk not in ("N", "a")} for k, v in T.items()}, rows=rows)
        print("47B L", L, res[str(L)]["tables"], flush=True)
        for key in (f"0.01|B1>B1|(6, 0, 0)|(0, 0)", f"0.01|B1>B1|(8, 0, 0)|(2, 2)", f"0.01|B1>B2|(6, 0, 0)|(0, 0)", f"0.01|C1>C1|(6, 0, 0)|(0, 0)"):
            if key in rows:
                print("  ", key, rows[key], flush=True)
        if L == 24:
            np.savez_compressed(OUT / "partial_F24.npz", **{f"{k}_{c}": T[k][c] for k in ("B1", "B2") for c in ("N", "a")})
    return res


# ---------------------------------------------------------------- 47C
def pair_band(S):
    L = S["L"]; tab = {}
    for w in wedge(L):
        b = bound(S, w); tab[w] = b["E0"] if b else np.nan
    return fill(L, tab).real


def exact_validation(L=6, G4pi=0.01, nev=5, nsolve=12):
    """exact two-body (heavy band, full pair space) in a static source field vs H_eff on the COM states."""
    from scipy.sparse.linalg import LinearOperator, eigsh
    A = json.loads((OUT / "A.json").read_text())["result"]; out = {}
    S = setup(L, *BENCH["B1"][:2], U=A["B1"]["U"]); N = L ** 3
    b0, FN, Fa, _ = pair_tables(S); Kq = kernel(L, G4pi)
    # source (B1, point at origin): Phi'(q) = Kq (F_N + F_a), Psi'(q) = Kq F_N  (e^{-iqv} convention)
    Phq, Psq = Kq * (FN + Fa), Kq * FN
    Phi_v = np.fft.ifftn(Phq).real                     # Phi_v = (1/N) sum_q Phi'(q) e^{iqv}
    # single-carrier field operator on the heavy band, k index flattened
    idx = np.arange(N).reshape(L, L, L); dh = np.zeros((N, N), complex)
    for nq in itertools.product(range(L), repeat=3):
        MN, Ma = M(S, nq), M(S, nq, True); kp = sh(idx, nq).reshape(-1)      # k+q index for each k
        dh[kp, idx.reshape(-1)] += (Phq[nq] * MN + Psq[nq] * Ma).reshape(-1) / N
    e = S["eh"].reshape(-1); u6 = S["uh"][..., 6].reshape(-1)
    kv = S["K"].reshape(-1, 3); sites = np.array(list(itertools.product(range(L), repeat=3)))
    Xv = u6[None, :] * np.exp(1j * sites @ kv.T)                               # X_{v,k} = u6(k) e^{ik.v}
    wv = 1 + Phi_v[tuple(sites.T)]
    def mv(x, field=True):
        p = x.reshape(N, N); y = e[:, None] * p + p * e[None, :]
        chi = np.einsum("vk,vk->v", Xv @ p, Xv) / N
        c = (wv if field else 1.0) * chi
        y += -S["U"] / N * (Xv.conj().T @ (c[:, None] * Xv.conj()))
        if field:
            y += dh @ p + p @ dh.T
        return y.reshape(-1)
    ex = {}
    for fl in (False, True):
        op = LinearOperator((N * N, N * N), matvec=lambda x, f=fl: mv(x, f), dtype=complex)
        ex[fl] = np.sort(eigsh(op, k=nsolve, which="SA", tol=1e-13, ncv=96)[0])[:nev]
    # H_eff: COM momentum basis
    EK = pair_band(S).reshape(-1)
    Kgrid = np.array(list(itertools.product(range(L), repeat=3)))
    He = np.diag(EK).astype(complex); H1 = np.diag(EK).astype(complex)
    bK = {tuple(k): bound(S, tuple(k)) for k in Kgrid}
    for i, Kp in enumerate(Kgrid):
        for j, K in enumerate(Kgrid):
            q = tuple((Kp - K) % L)
            He[i, j] += (Phq[q] * FN[q] + Psq[q] * Fa[q]) / N
            F = formfactor(S, q, bK[tuple(K)], bK[tuple(Kp)], tuple(K))
            H1[i, j] += (Phq[q] * F["N"] + Psq[q] * F["a"]) / N
    ee = np.sort(np.linalg.eigvalsh(0.5 * (He + He.conj().T)))[:nev]; e1 = np.sort(np.linalg.eigvalsh(0.5 * (H1 + H1.conj().T)))[:nev]
    e0 = np.sort(EK)[:nev]
    out = dict(L=L, G4pi=G4pi, exact_nofield=ex[False].tolist(), pair_band_lowest=e0.tolist(), exact_field=ex[True].tolist(),
               Heff=ee.tolist(), H1st=e1.tolist(), shift_exact=(ex[True] - ex[False]).tolist(), shift_Heff=(ee - e0).tolist(),
               shift_H1st=(e1 - e0).tolist(), herm_Heff=float(np.abs(He - He.conj().T).max()), herm_H1=float(np.abs(H1 - H1.conj().T).max()))
    out["rel_Heff"] = (np.abs(np.array(out["shift_Heff"]) / np.array(out["shift_exact"]) - 1)).tolist()
    out["rel_H1st"] = (np.abs(np.array(out["shift_H1st"]) / np.array(out["shift_exact"]) - 1)).tolist()
    out["C_nofield"] = float(np.abs(ex[False] - e0).max())
    print("47C validation", out, flush=True)
    return out


def evolve(EK, V, c0, dt, nsteps, every):
    L = EK.shape[0]; ph = np.exp(-1j * (EK - EK.flat[0]) * dt); pv = np.exp(-0.5j * V * dt)
    n = np.arange(L); d = (n + L // 2) % L - L // 2; xg = np.meshgrid(d, d, d, indexing="ij")[0]
    c = c0.copy(); xs = [(0.0, float(np.sum(np.abs(c) ** 2 * xg)))]; norm0 = np.sum(np.abs(c) ** 2)
    for s in range(1, nsteps + 1):
        c = pv * c; c = np.fft.ifftn(ph * np.fft.fftn(c)); c = pv * c
        if s % every == 0:
            xs.append((s * dt, float(np.sum(np.abs(c) ** 2 * xg))))
    return np.array(xs), float(abs(np.sum(np.abs(c) ** 2) / norm0 - 1)), c


def part_C(L=24, G4pi=0.01, sigma=3.0, seps=(6, 8)):
    A = json.loads((OUT / "A.json").read_text())["result"]; res = dict(validation=exact_validation())
    Fz = np.load(OUT / "partial_F24.npz"); F = {k: dict(N=Fz[f"{k}_N"], a=Fz[f"{k}_a"]) for k in ("B1", "B2")}
    Kq = kernel(L, G4pi); EK = {}
    for name, (m0, g, fr) in BENCH.items():
        EK[name] = pair_band(setup(L, m0, g, U=A[name]["U"]))
    n = np.arange(L); d = (n + L // 2) % L - L // 2; Dg = np.meshgrid(d, d, d, indexing="ij")
    Phi_S = np.fft.ifftn(Kq * (F["B1"]["N"] + F["B1"]["a"])).real      # source B1 lapse field (point, origin)
    out = {}
    for probe in ("B1", "B2"):
        V = np.fft.ifftn(Vtilde(F["B1"], F[probe], Kq)).real            # V(R) for source B1 at 0, point probe at R
        t = A[probe]["t_pair"]
        for d0 in seps:
            c0 = np.exp(-((Dg[0] - d0) ** 2 + Dg[1] ** 2 + Dg[2] ** 2) / (4 * sigma ** 2)).astype(complex); c0 /= np.linalg.norm(c0)
            w = np.abs(c0) ** 2
            gV = (np.roll(V, -1, 0) - np.roll(V, 1, 0)) / 2; gP = (np.roll(Phi_S, -1, 0) - np.roll(Phi_S, 1, 0)) / 2
            a_pred = float(-2 * t * np.sum(w * gV)); a_N = float(-np.sum(w * gP))
            T = np.sqrt(2 * 0.05 / abs(a_pred)); runs = {}
            for dtf in (1, 2):
                nst = 400 * dtf; dt = T / nst
                xs, nerr, _ = evolve(EK[probe], V, c0, dt, nst, 4 * dtf)
                xo, _, _ = evolve(EK[probe], 0 * V, c0, dt, nst, 4 * dtf)
                dx = xs[:, 1] - xo[:, 1]; tt = xs[:, 0]
                cf = np.linalg.lstsq(np.stack([tt ** 2 / 2, tt ** 3, tt ** 4], 1), dx, rcond=None)[0]
                runs[dtf] = dict(a_meas=float(cf[0]), dx_end=float(dx[-1]), norm_err=nerr, x_off_drift=float(xo[-1, 1] - xo[0, 1]))
            out[f"{probe}@{d0}"] = dict(t_pair=t, a_pred=a_pred, a_N=a_N, a_meas=runs[2]["a_meas"], rel_pred=runs[2]["a_meas"] / a_pred - 1,
                                        a_over_aN=runs[2]["a_meas"] / a_N, dt_change=abs(runs[1]["dx_end"] / runs[2]["dx_end"] - 1),
                                        runs=runs, T=T, Vvar_across_pair=float(abs(gV[d0, 0, 0]) * 2 * A[probe]["rms_sep"]), Eb=A[probe]["Eb"])
            print("47C probe", probe, d0, {k: v for k, v in out[f"{probe}@{d0}"].items() if k != "runs"}, flush=True)
    res["probe"] = out
    # both released: relative coordinate, P = 0
    rel = {}
    for (s, p) in (("B1", "B1"), ("B1", "B2")):
        Erel = EK[s] + neg(EK[p]); V = np.fft.ifftn(Vtilde(F[s], F[p], Kq)).real; d0 = 8
        c0 = np.exp(-((Dg[0] - d0) ** 2 + Dg[1] ** 2 + Dg[2] ** 2) / (4 * sigma ** 2)).astype(complex); c0 /= np.linalg.norm(c0)
        gV = (np.roll(V, -1, 0) - np.roll(V, 1, 0)) / 2
        tr = A[s]["t_pair"] + A[p]["t_pair"]; a_pred = float(-2 * tr * np.sum(np.abs(c0) ** 2 * gV)); T = np.sqrt(2 * 0.05 / abs(a_pred))
        rr = {}
        for dtf in (1, 2):
            nst = 400 * dtf; xs, nerr, _ = evolve(Erel, V, c0, T / nst, nst, 4 * dtf); xo, _, _ = evolve(Erel, 0 * V, c0, T / nst, nst, 4 * dtf)
            dx = xs[:, 1] - xo[:, 1]; tt = xs[:, 0]; cf = np.linalg.lstsq(np.stack([tt ** 2 / 2, tt ** 3, tt ** 4], 1), dx, rcond=None)[0]
            rr[dtf] = dict(a_meas=float(cf[0]), dx_end=float(dx[-1]), norm_err=nerr)
        rel[f"{s}+{p}"] = dict(a_pred=a_pred, a_meas=rr[2]["a_meas"], rel_pred=rr[2]["a_meas"] / a_pred - 1, dx_end=rr[2]["dx_end"],
                               dt_change=abs(rr[1]["dx_end"] / rr[2]["dx_end"] - 1), runs=rr)
        print("47C both", s, p, {k: v for k, v in rel[f"{s}+{p}"].items() if k != "runs"}, flush=True)
    res["both"] = rel
    return res


def part_checks(L=16):
    A = json.loads((OUT / "A.json").read_text())["result"]; S = setup(L, *BENCH["B1"][:2], U=A["B1"]["U"]); b0 = bound(S)
    rng = np.random.default_rng(47); out = []
    _, FN, Fa, _ = pair_tables(S)
    for _ in range(6):
        nq = tuple(int(x) for x in rng.integers(0, L, 3)); F = formfactor(S, nq, b0)
        out.append(dict(q=nq, dN=float(abs(F["N"] - FN[nq])), da=float(abs(F["a"] - Fa[nq]))))
    print("C2", out, flush=True)
    return dict(C2=out)


if __name__ == "__main__":
    part = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = {"A": part_A, "B": part_B, "C": part_C, "checks": part_checks}[part]()
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else (str(o) if isinstance(o, complex) else float(o))) + "\n")
    print("done", part, res["wall_seconds"])
