"""EX-GC47R: corrections to EX-GC47 after the Codex review. (1) projector-derivative (Hellmann-Feynman consistent)
form factors; (2) release packet periodic about its own centre; (3) force vs distance and box sweep.
Plan: spec/EX-GC47R_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time
from pathlib import Path
import numpy as np
import exgc45 as X
import exgc47 as D

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC47R"
A47 = Path(__file__).resolve().parents[1] / "results" / "EX-GC47" / "A.json"
sh, neg = D.sh, D.neg


def setup_full(L, name):
    A = json.loads(A47.read_text())["result"][name]; m0, g, _ = D.BENCH[name]
    S = D.setup(L, m0, g, U=A["U"]); e, u, _ = X.eig(L, D.LAM0, m0=m0, g=g); S["e_all"], S["u_all"] = e, u
    return S, A


def Mbb(S, nq, bto, bfrom, strain):
    """<k+q, bto| o(q) |k, bfrom> over k."""
    u = S["u_all"]
    return np.sum(np.conj(sh(u[..., :, bto], nq)) * D.Bu(S, u[..., :, bfrom], nq, strain), axis=-1)


def extra(S, nq, b0, bq, strain, bands):
    """projector-derivative terms of <Psi_q| d(P H P) |Psi_0> for a field coupling o(q), per carrier."""
    L = S["L"]; N = L ** 3; U = S["U"]; e = S["e_all"]; u6 = S["u_all"][..., 6, :]
    e3, u3 = e[..., 3], u6[..., 3]; p0, pq = b0["phi"], bq["phi"]; tot = 0j
    for b in bands:
        eb, ub = e[..., b], u6[..., b]
        Mb3, M3b = Mbb(S, nq, b, 3, strain), Mbb(S, nq, 3, b, strain)
        # carrier 1
        A1 = -(U / N) * np.conj(bq["a"]) * np.sum(sh(ub, nq) * neg(u3) * Mb3 * p0 / (e3 - sh(eb, nq)))
        B1 = -(U / N) * b0["a"] * np.sum(np.conj(pq) * M3b * np.conj(ub * neg(u3)) / (sh(e3, nq) - eb))
        # carrier 2
        A2 = -(U / N) * np.conj(bq["a"]) * np.sum(u3 * neg(sh(ub, nq)) * neg(Mb3) * p0 / (neg(e3) - neg(sh(eb, nq))))
        B2 = -(U / N) * b0["a"] * np.sum(np.conj(np.roll(pq, list(nq), axis=(0, 1, 2))) * neg(M3b) * np.conj(u3 * neg(ub)) / (neg(sh(e3, nq)) - neg(eb)))
        tot += A1 + B1 + A2 + B2
    return tot


def tables_corr(S, bands):
    L = S["L"]; b0 = D.bound(S); FN, Fa, FN0, Fa0 = {}, {}, {}, {}
    for w in D.wedge(L):
        bq = D.bound(S, w); F = D.formfactor(S, w, b0, bq)
        FN0[w], Fa0[w] = F["N"], F["a"]
        FN[w] = F["N"] + extra(S, w, b0, bq, False, bands); Fa[w] = F["a"] + extra(S, w, b0, bq, True, bands)
    return b0, D.fill(L, FN), D.fill(L, Fa), D.fill(L, FN0), D.fill(L, Fa0)


def part_F(h=1e-5):
    """Q1 at L = 16 and 32 (q = 0), and corrected tables at L = 24 and 32."""
    res = {}; ALL = (0, 1, 2, 4, 5, 6); UP = (4, 5, 6)
    for name in ("B1", "B2"):
        for L in (16, 32):
            S, A = setup_full(L, name); b0 = D.bound(S)
            F0 = D.formfactor(S, (0, 0, 0), b0, b0)
            ea, eu = extra(S, (0, 0, 0), b0, b0, True, ALL), extra(S, (0, 0, 0), b0, b0, True, UP)
            eN = extra(S, (0, 0, 0), b0, b0, False, ALL)
            Ep = D.bound(dict(S, eh=X.eig(L, D.LAM0, m0=S["m0"], g=S["g"], t=(1 + h,) * 3)[0][..., 3], uh=X.eig(L, D.LAM0, m0=S["m0"], g=S["g"], t=(1 + h,) * 3)[1][..., :, 3]))["E0"]
            Em = D.bound(dict(S, eh=X.eig(L, D.LAM0, m0=S["m0"], g=S["g"], t=(1 - h,) * 3)[0][..., 3], uh=X.eig(L, D.LAM0, m0=S["m0"], g=S["g"], t=(1 - h,) * 3)[1][..., :, 3]))["E0"]
            tdE = (Ep - Em) / (2 * h)
            res[f"Q1_{name}_L{L}"] = dict(Fa0_proj=F0["a"].real, extra_all=ea.real, extra_up=eu.real, Fa0_corr_all=(F0["a"] + ea).real,
                                          Fa0_corr_up=(F0["a"] + eu).real, t_dEdt=tdE, rel_all=abs((F0["a"] + ea).real / tdE - 1), extra_N_q0=abs(eN))
            print("Q1", name, L, res[f"Q1_{name}_L{L}"], flush=True)
    for L in (24, 32):
        T = {}
        for name in ("B1", "B2"):
            S, A = setup_full(L, name); b0, FN, Fa, FN0, Fa0 = tables_corr(S, ALL)
            T[name] = dict(N=FN, a=Fa, N0=FN0, a0=Fa0)
        np.savez_compressed(OUT / f"partial_Fcorr{L}.npz", **{f"{k}_{c}": T[k][c] for k in T for c in ("N", "a", "N0", "a0")})
        Kq = D.kernel(L, 0.01); rows = {}
        for (s, p) in (("B1", "B1"), ("B1", "B2")):
            for r in range(2, L // 2 + 1):
                vc = D.V12(dict(N=T[s]["N"], a=T[s]["a"]), dict(N=T[p]["N"], a=T[p]["a"]), Kq, (r, 0, 0)).real
                v0 = D.V12(dict(N=T[s]["N0"], a=T[s]["a0"]), dict(N=T[p]["N0"], a=T[p]["a0"]), Kq, (r, 0, 0)).real
                rows[f"{s}>{p}|{r}"] = dict(V_corr=vc, V_47=v0, change=vc / v0 - 1)
        res[f"V_L{L}"] = rows
        ch = [abs(v["change"]) for k, v in rows.items() if 3 <= int(k.split("|")[1]) <= 12]
        print("Q2 L", L, "max |change| r 3-12:", max(ch), flush=True)
    return res


def evolve(EK, V, c0, xg, dt, nsteps, every):
    L = EK.shape[0]; E0 = EK.flat[0]; ph = np.exp(-1j * (EK - E0) * dt); pv = np.exp(-0.5j * V * dt)
    def energy(c):
        ck = np.fft.fftn(c); return float((np.sum(np.abs(ck) ** 2 * (EK - E0)) / L ** 3 + np.sum(np.abs(c) ** 2 * V)) / np.sum(np.abs(c) ** 2))
    c = c0.copy(); xs = [(0.0, float(np.sum(np.abs(c) ** 2 * xg)))]; en0 = energy(c)
    for s in range(1, nsteps + 1):
        c = pv * c; c = np.fft.ifftn(ph * np.fft.fftn(c)); c = pv * c
        if s % every == 0:
            xs.append((s * dt, float(np.sum(np.abs(c) ** 2 * xg))))
    return np.array(xs), energy(c) - en0, float(abs(np.sum(np.abs(c) ** 2) - np.sum(np.abs(c0) ** 2)))


def part_R(G4pi=0.01):
    F = json.loads((OUT / "F.json").read_text()) if (OUT / "F.json").exists() else None
    A = {k: json.loads(A47.read_text())["result"][k] for k in ("B1", "B2")}; res = {}
    for L in (24, 32):
        Z = np.load(OUT / f"partial_Fcorr{L}.npz"); T = {k: dict(N=Z[f"{k}_N"], a=Z[f"{k}_a"]) for k in ("B1", "B2")}
        Kq = D.kernel(L, G4pi); EK = {k: D.pair_band(setup_full(L, k)[0]) for k in ("B1", "B2")}
        Phi_S = np.fft.ifftn(Kq * (T["B1"]["N"] + T["B1"]["a"])).real
        n = np.arange(L)
        for mode in ("probe", "both"):
            for probe in ("B1", "B2"):
                V = np.fft.ifftn(D.Vtilde(T["B1"], T[probe], Kq)).real
                Ek = EK[probe] if mode == "probe" else EK["B1"] + neg(EK[probe])
                tp = A[probe]["t_pair"] if mode == "probe" else A["B1"]["t_pair"] + A[probe]["t_pair"]
                for d0 in (6, 8):
                    for sigma in (2.0, 3.0):
                        dx = (n - d0 + L // 2) % L - L // 2                  # periodic about the packet's own centre
                        Dg = np.meshgrid(dx, (n + L // 2) % L - L // 2, (n + L // 2) % L - L // 2, indexing="ij")
                        c0 = np.exp(-(Dg[0] ** 2 + Dg[1] ** 2 + Dg[2] ** 2) / (4 * sigma ** 2)).astype(complex); c0 /= np.linalg.norm(c0)
                        xg = d0 + Dg[0]                                       # centroid coordinate, unwrapped about d0
                        w = np.abs(c0) ** 2; gV = (np.roll(V, -1, 0) - np.roll(V, 1, 0)) / 2; gP = (np.roll(Phi_S, -1, 0) - np.roll(Phi_S, 1, 0)) / 2
                        a_pred = float(-2 * tp * np.sum(w * gV)); a_N = float(-np.sum(w * gP)); T0 = np.sqrt(2 * 0.05 / abs(a_pred))
                        Vscale = float(abs(np.sum(w * V)))
                        for win in (1.0, 0.5):
                            Tw = T0 * win; out = {}
                            for dtf in (1, 2):
                                nst = 400 * dtf; xs, dE, dn = evolve(Ek, V, c0, xg, Tw / nst, nst, 4 * dtf)
                                xo, _, _ = evolve(Ek, 0 * V, c0, xg, Tw / nst, nst, 4 * dtf)
                                dxs = xs[:, 1] - xo[:, 1]; tt = xs[:, 0]
                                cf = np.linalg.lstsq(np.stack([tt ** 2 / 2, tt ** 3, tt ** 4], 1), dxs, rcond=None)[0]
                                out[dtf] = dict(a=float(cf[0]), dx_end=float(dxs[-1]), dE=dE, dE_rel=abs(dE) / Vscale, dnorm=dn,
                                                off_drift=float(xo[-1, 1] - xo[0, 1]))
                            key = f"L{L}|{mode}|{probe}|d{d0}|s{sigma}|w{win}"
                            res[key] = dict(a_meas=out[2]["a"], a_pred=a_pred, rel_pred=out[2]["a"] / a_pred - 1, a_N=a_N,
                                            a_over_aN=out[2]["a"] / a_N if mode == "probe" else None, dt_change=abs(out[1]["dx_end"] / out[2]["dx_end"] - 1),
                                            dE_rel=out[2]["dE_rel"], off_drift=out[2]["off_drift"], dx_end=out[2]["dx_end"], seam_weight=float(w[(d0 + L // 2) % L, 0, 0] / w.max()))
                            print(key, {k: (round(v, 6) if isinstance(v, float) and abs(v) > 1e-3 else v) for k, v in res[key].items()}, flush=True)
    return res


def part_force(G4pi=0.01):
    res = {}; Geff = G4pi / (4 * np.pi * (1 + G4pi * D.C2CHI))
    for L in (24, 32, 48):
        if L == 48:
            T = {}
            for name in ("B1",):
                S, A = setup_full(L, name); b0, FN, Fa, _, _ = tables_corr(S, (0, 1, 2, 4, 5, 6)); T[name] = dict(N=FN, a=Fa)
        else:
            Z = np.load(OUT / f"partial_Fcorr{L}.npz"); T = {"B1": dict(N=Z["B1_N"], a=Z["B1_a"])}
        F = T["B1"]; E, Q = F["N"][0, 0, 0].real, F["a"][0, 0, 0].real; coef = E * E + 2 * E * Q
        Kq = D.kernel(L, G4pi)
        V = np.fft.ifftn(D.Vtilde(F, F, Kq)).real
        P = dict(N=np.full_like(F["N"], E), a=np.full_like(F["a"], Q)); Vp = np.fft.ifftn(D.Vtilde(P, P, Kq)).real
        rows = []
        for r in range(2, L // 2):
            f = -(V[r + 1, 0, 0] - V[r - 1, 0, 0]) / 2; fp = -(Vp[r + 1, 0, 0] - Vp[r - 1, 0, 0]) / 2
            f_inf = -Geff * coef * (1 / (r + 1) - 1 / (r - 1)) / 2          # same central difference of -G coef / r
            rows.append(dict(r=r, F=f, F_point=fp, ratio_point=f / fp, r2F_over_Gcoef=-(r ** 2) * f / (Geff * coef) * -1, ratio_inf=f / f_inf))
        res[f"L{L}"] = dict(E=E, Qa=Q, rows=rows)
        print("force L", L, [(x["r"], round(x["ratio_point"], 5), round(x["ratio_inf"], 4)) for x in rows], flush=True)
    return res


if __name__ == "__main__":
    part = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = {"F": part_F, "R": part_R, "force": part_force}[part]()
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)) + "\n")
    print("done", part, res["wall_seconds"])
