"""EX-GC46: heavy pairs with all four upper bands (rank-1 contact in a larger basis), and the small-g limit.
46A: 45's nine points with four bands per carrier (fixed U and matched E_b/W). 46B: small-g sweep at two m0.
Plan: spec/EX-GC46_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import exgc45 as X

OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC46"
C45 = Path(__file__).resolve().parents[1] / "results" / "EX-GC45"
LAM0 = 0.1
neg = lambda a: np.roll(np.flip(a, axis=(0, 1, 2)), 1, axis=(0, 1, 2))


def basis(e, u, Ks=(0, 0, 0), bands=(3, 4, 5, 6)):
    D = []; P = []
    for b1 in bands:
        e1 = np.roll(e[..., b1], [-s for s in Ks], axis=(0, 1, 2)); h1 = np.roll(u[..., 6, b1], [-s for s in Ks], axis=(0, 1, 2))
        for b2 in bands:
            D.append((e1 + neg(e[..., b2])).reshape(-1)); P.append((np.conj(h1) * np.conj(neg(u[..., 6, b2]))).reshape(-1))
    return np.concatenate(D), np.concatenate(P)


def solve(e, u, U, Ks=(0, 0, 0), bands=(3, 4, 5, 6)):
    N = e.shape[0] ** 3; D, P = basis(e, u, Ks, bands); A = np.abs(P) ** 2; Dmin = D.min()
    f = lambda E: U * np.sum(A / (D - E)) / N - 1
    if f(Dmin - 1e-14) < 0:
        return None
    E0 = brentq(f, Dmin - 50, Dmin - 1e-14, xtol=1e-15, rtol=1e-14)
    w = A / (D - E0) ** 2; Dexp = float(np.sum(w * D) / np.sum(w))
    return dict(E0=float(E0), Dmin=float(Dmin), Eb=float(Dmin - E0), V=float(E0 - Dexp), T=float(Dexp - Dmin))


def U_for(e, u, Eb, bands=(3, 4, 5, 6)):
    N = e.shape[0] ** 3; D, P = basis(e, u, bands=bands)
    return float(N / np.sum(np.abs(P) ** 2 / (D - (D.min() - Eb))))


def point(L, m0, g, U=None, frac=None, h=1e-5, bands=(3, 4, 5, 6), mob=True):
    e, u, _ = X.eig(L, LAM0, m0=m0, g=g); W = float(e[..., 3].max() - e[..., 3].min())
    if U is None:
        U = U_for(e, u, frac * W, bands)
    r0 = solve(e, u, U, bands=bands)
    ep, up, _ = X.eig(L, LAM0, m0=m0, g=g, t=(1 + h,) * 3); em, um, _ = X.eig(L, LAM0, m0=m0, g=g, t=(1 - h,) * 3)
    rp = solve(ep, up, U, bands=bands); rm = solve(em, um, U, bands=bands)
    Qp = (rp["E0"] - rm["E0"]) / (2 * h); Qe = (rp["Dmin"] - rm["Dmin"]) / (2 * h); dQ = Qp - Qe
    row = dict(m0=m0, g=g, W=W, U=U, E0=r0["E0"], Eb=r0["Eb"], Eb_over_W=r0["Eb"] / W, Eb_over_dE=r0["Eb"] / r0["E0"], V=r0["V"], T=r0["T"],
               Qa_pair=Qp, Qa_thr=Qe, dQa=dQ, eta=dQ / r0["Eb"], Qa_over_dE=Qp / r0["E0"], dQa_over_dE=dQ / r0["E0"])
    if mob:
        Kx = 2 * np.pi * 2 / L; rK = solve(e, u, U, (2, 0, 0), bands)
        row["t_pair"] = (rK["E0"] - r0["E0"]) / Kx ** 2 if rK else None
    return row


def part_A(L=32):
    C = json.loads((C45 / "C.json").read_text())["result"]; res = {}
    for key, c in C.items():
        m0, g = c["m0"], c["g"]; fr = float(key.split("_")[2])
        fx = point(L, m0, g, U=c["U"], mob=False); mt = point(L, m0, g, frac=fr)
        res[key] = dict(proj_eta=c["eta"], proj_Eb=c["Eb"], proj_dQa=c["dQa"], fixedU=fx, matched=mt,
                        Eb_change=fx["Eb"] / c["Eb"] - 1, dQa_change=fx["dQa"] / c["dQa"] - 1 if c["dQa"] else None)
        print("46A", key, "eta proj %.4f fixedU %.4f matched %.4f | dEb %.3g" % (c["eta"], fx["eta"], mt["eta"], res[key]["Eb_change"]), flush=True)
    return res


def part_B(L=32, m0s=(0.5, 0.7), gs=(0.0125, 0.025, 0.05, 0.1), fracs=(0.15, 0.30)):
    res = {}
    for m0 in m0s:
        for g in gs:
            for fr in fracs:
                r = point(L, m0, g, frac=fr); res[f"{m0}_{g}_{fr}"] = r
                print("46B", m0, g, fr, "eta %.4f QaOverdE %.5f EbOverdE %.2e tpair %s" % (r["eta"], r["Qa_over_dE"], r["Eb_over_dE"], r["t_pair"]), flush=True)
    return res


def part_checks(L=32, h=1e-5):
    C = json.loads((C45 / "C.json").read_text())["result"]["0.5_0.1_0.15"]; U = C["U"]; out = {}
    # C2: band-3-only equals 45's projected solver; HF
    e, u, _ = X.eig(L, LAM0, m0=0.5, g=0.1)
    a = solve(e, u, U, bands=(3,)); b = X.pair_solve(e[..., 3], u[..., :, 3], U)
    rU = solve(e, u, U * (1 + 1e-6)); rUm = solve(e, u, U * (1 - 1e-6)); r4 = solve(e, u, U)
    out["C2"] = dict(diff_E0=abs(a["E0"] - b["E0"]), HF_rel=abs((rU["E0"] - rUm["E0"]) / 2e-6 / r4["V"] - 1))
    # C1: homogeneity
    E = {}
    for s in (1, -1):
        H, _ = X.bloch(L, LAM0, m0=0.5, g=0.1, Nsc=1 + s * h); w, v = np.linalg.eigh(H); E[s] = solve(w, v, U * (1 + s * h))["E0"]
    out["C1"] = dict(QN=(E[1] - E[-1]) / (2 * h), dE=r4["E0"], rel=abs((E[1] - E[-1]) / (2 * h) / r4["E0"] - 1))
    # C3: dense four-band L = 6 against the solver
    e6, u6, _ = X.eig(6, LAM0, m0=0.5, g=0.1)
    out["C3"] = dict(dense=float(X.pair_dense_full(6, LAM0, 0.5, 0.1, U, (3, 4, 5, 6))), solver=solve(e6, u6, U)["E0"])
    out["C3"]["diff"] = abs(out["C3"]["dense"] - out["C3"]["solver"])
    for k, v in out.items():
        print(k, v, flush=True)
    return out


if __name__ == "__main__":
    part = sys.argv[1]; OUT.mkdir(parents=True, exist_ok=True); t0 = time.time()
    res = {"A": part_A, "B": part_B, "checks": part_checks}[part]()
    res = dict(result=res, wall_seconds=time.time() - t0, platform=f"python {sys.version.split()[0]}, numpy {np.__version__}")
    (OUT / f"{part}.json").write_text(json.dumps(res, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else float(o)) + "\n")
    print("done", part, res["wall_seconds"])
