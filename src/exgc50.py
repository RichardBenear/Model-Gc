"""EX-GC50: vacuum kernel chi_ij(q), i, j in {N, a}, of the 40A pair-network vacuum and the 48 Wilson-Dirac sea under 44's
homogeneous clock weights (bonds Nbar, sites N; bonds ~ 1/abar), vs 37B weights (bonds N_v N_w; anchor to EX-GC43).
Second-order perturbation theory for the filled sea: paramagnetic -2 Re S_ij / N_s plus diamagnetic weight terms.
Plan: spec/EX-GC50_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc40a as A40
import exgc48 as D48

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC50"
J6 = np.ones((6, 6)) - np.eye(6)
sh = lambda a, n: np.roll(a, [-x for x in n], axis=(0, 1, 2))


# ------------------------------------------------------------------ network
def net_setup(L, lam0):
    H, K = A40.bloch_h(L, lam0); e, u = np.linalg.eigh(H)
    return dict(L=L, lam0=lam0, K=K, e=e, u=u, H=H)


def net_bond(K, a, g):
    """Bloch transition operator <k+q| bond_a weighted |k>, as (L,L,L,6,6); g = weight factor for this axis."""
    L = K.shape[0]; O = np.zeros(K.shape[:3] + (6, 6), complex); ph = np.exp(1j * K[..., a])
    O[..., 2 * a, 2 * a + 1] = -ph * g; O[..., 2 * a + 1, 2 * a] = -np.conj(ph) * np.conj(g)
    return O


def net_kernel(S, nq, weights="H"):
    L = S["L"]; K = S["K"]; q = 2 * np.pi * np.array(nq) / L; N = L ** 3
    gs = [(1 + np.exp(1j * q[a])) / 2 for a in range(3)]
    bonds = [net_bond(K, a, gs[a]) for a in range(3)]
    fac = 2.0 if weights == "37B" else 1.0
    ON = S["lam0"] * J6 + fac * sum(bonds); Oa = -sum(bonds)
    e, u = S["e"], S["u"]; eq, uq = sh(e, nq), sh(u, nq)
    def Mx(O):
        return np.einsum("...ib,...ij,...jc->...bc", uq.conj(), np.broadcast_to(O, u.shape), u)[..., 3:, :3]   # emp(k+q) x occ(k)
    Dl = eq[..., 3:, None] - e[..., None, :3]
    MN, Ma = Mx(ON), Mx(Oa)
    Snn = np.sum(np.abs(MN) ** 2 / Dl); Sna = np.sum(np.conj(MN) * Ma / Dl); Saa = np.sum(np.abs(Ma) ** 2 / Dl)
    ea = bond_energies_net(S)
    c = np.cos(q)
    out = dict(para_NN=float(-2 * Snn.real / N), para_Na=float(-2 * Sna.real / N), para_aa=float(-2 * Saa.real / N))
    if weights == "37B":
        out["dia_NN"] = float(2 * np.sum(ea * c)); out["dia_Na"] = out["dia_aa"] = 0.0
    else:
        out["dia_NN"] = 0.0; out["dia_Na"] = float(-0.5 * np.sum(ea * (1 + c))); out["dia_aa"] = float(np.sum(ea * (1 + c)))
    for k in ("NN", "Na", "aa"):
        out[k] = out["para_" + k] + out["dia_" + k]
    out["qt2"] = float(np.sum(2 * (1 - c)))
    return out


def bond_energies_net(S):
    if "ea" not in S:
        K = S["K"]; u = S["u"][..., :, :3]; N = S["L"] ** 3; ea = []
        for a in range(3):
            O = net_bond(K, a, 1.0); ea.append(float(np.real(np.sum(np.einsum("...ib,...ij,...jb->...", u.conj(), O, u))) / N))
        S["ea"] = np.array(ea)
    return S["ea"]


def net_uniform(L, lam0, h=1e-4):
    def Evac(aa):
        H, K = A40.bloch_h(L, lam0); Hs = lam0 * J6; Hl = H - Hs
        e = np.linalg.eigvalsh(Hl / aa + Hs); return float(np.sum(e[..., :3]) / L ** 3)
    E0, Ep, Em = Evac(1.0), Evac(1 + h), Evac(1 - h)
    return dict(E_vac=E0, dE_da=(Ep - Em) / (2 * h), d2E_da2=(Ep + Em - 2 * E0) / h ** 2)


# ------------------------------------------------------------------ Dirac sea (W1, Wilson as bond operators)
G = D48.GAM


def dirac_bond(K, a, g, r, q_a):
    ka = K[..., a]; al = G[1 + a]; be = G[4]
    t1 = np.exp(1j * ka)[..., None, None] * al / 2j - np.exp(-1j * (ka + q_a))[..., None, None] * al / 2j
    t2 = (r / 2) * (1 + np.exp(-1j * q_a)) * be - (r / 2) * (np.exp(1j * ka) + np.exp(-1j * (ka + q_a)))[..., None, None] * be
    return g * (t1 + t2)


def dirac_proj(K, m, r):
    sk, M, E = D48.parts(K, m, r)
    h = np.einsum("...j,jab->...ab", sk, np.array(G[1:4])) + M[..., None, None] * G[4]
    I = np.eye(4); return 0.5 * (I + h / E[..., None, None]), 0.5 * (I - h / E[..., None, None]), E


def dirac_kernel(L, m, nq, r=1.0, cache={}):
    key = (L, m, r)
    if key not in cache:
        k = 2 * np.pi * np.fft.fftfreq(L); K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
        Pp, Pm, E = dirac_proj(K, m, r); ea = []
        for a in range(3):
            B0 = dirac_bond(K, a, 1.0, r, 0.0); ea.append(float(np.real(np.sum(np.einsum("...ij,...ji->...", Pm, B0))) / L ** 3))
        cache[key] = (K, Pp, Pm, E, np.array(ea))
    K, Pp, Pm, E, ea = cache[key]; q = 2 * np.pi * np.array(nq) / L; N = L ** 3
    gs = [(1 + np.exp(1j * q[a])) / 2 for a in range(3)]
    bonds = [dirac_bond(K, a, gs[a], r, q[a]) for a in range(3)]
    ON = m * G[4] + sum(bonds); Oa = -sum(bonds)
    Ppq = sh(Pp, nq); Eq = sh(E, nq); Dl = (Eq + E)[..., None, None]
    def S(Oi, Oj):
        return np.sum(np.einsum("...ab,...bc,...cd,...ed->...ae", Ppq, Oi, Pm, Oj.conj()).trace(axis1=-2, axis2=-1) / Dl[..., 0, 0])
    Snn, Sna, Saa = S(ON, ON), S(ON, Oa), S(Oa, Oa); c = np.cos(q)
    out = dict(para_NN=float(-2 * Snn.real / N), para_Na=float(-2 * Sna.real / N), para_aa=float(-2 * Saa.real / N),
               dia_NN=0.0, dia_Na=float(-0.5 * np.sum(ea * (1 + c))), dia_aa=float(np.sum(ea * (1 + c))), qt2=float(np.sum(2 * (1 - c))))
    for kk in ("NN", "Na", "aa"):
        out[kk] = out["para_" + kk] + out["dia_" + kk]
    out["e_a"] = ea.tolist()
    return out


def dirac_uniform(L, m, r=1.0, h=1e-4):
    k = 2 * np.pi * np.fft.fftfreq(L); K = np.stack(np.meshgrid(k, k, k, indexing="ij"), -1)
    def Evac(aa):
        S = np.sin(K) / aa; W = np.sum(1 - np.cos(K), -1); E = np.sqrt(np.sum(S ** 2, -1) + (m + r * W / aa) ** 2)
        return float(-2 * np.sum(E) / L ** 3)
    E0, Ep, Em = Evac(1.0), Evac(1 + h), Evac(1 - h)
    return dict(E_vac=E0, dE_da=(Ep - Em) / (2 * h), d2E_da2=(Ep + Em - 2 * E0) / h ** 2)


# ------------------------------------------------------------------ fits and parts
QS_AX = [(mm, 0, 0) for mm in range(1, 9)]; QS_DG = [(1, 1, 1), (2, 2, 2)]


def fitq(rows, key, nfit=4):
    x = np.array([r["qt2"] for r in rows[:nfit]]); y = np.array([r[key] for r in rows[:nfit]])
    c = np.polyfit(x, y, 2)                     # y = c0 x^2 + c1 x + c2  ->  chi0 = c2, c2coef = c1
    return dict(chi0=float(c[2]), c2=float(c[1]), c4=float(c[0]))


def part_kernels(L=32):
    res = {}; C43 = json.loads((ROOT / "results" / "EX-GC43" / "vacuum.json").read_text())["result"]
    for lam0 in (0.1, 0.2):
        S = net_setup(L, lam0); uni = net_uniform(L, lam0); ea = bond_energies_net(S)
        rowsH = [dict(nq=list(n), **net_kernel(S, n, "H")) for n in QS_AX + QS_DG]
        rows37 = [dict(nq=list(n), **net_kernel(S, n, "37B")) for n in QS_AX]
        anchor = [rows37[i]["NN"] / C43[str(lam0)]["axis"][str(i + 1)]["chi"] - 1 for i in range(8)]
        fN = fitq(rowsH, "NN"); fNa = fitq(rowsH, "Na"); faa = fitq(rowsH, "aa"); f37 = fitq(rows37, "NN")
        res[f"net_{lam0}"] = dict(uniform=uni, eps_link=float(ea.sum()), e_a=ea.tolist(), H=rowsH, B37=rows37, anchor_rel=anchor,
                                  fit_NN=fN, fit_Na=fNa, fit_aa=faa, fit_37B=f37,
                                  Na0_rel=fNa["chi0"] / (-ea.sum()) - 1, aa0_rel=faa["chi0"] / uni["d2E_da2"] - 1,
                                  NN_m1_over_37B=rowsH[0]["NN"] / rows37[0]["NN"], NN_max=max(r["NN"] for r in rowsH))
        print("net", lam0, "anchor max %.2e" % max(map(abs, anchor)), "| NN chi0 %.2e c2 %.4f max %.2e | m1 ratio %.2e | Na0 %.6f (-eps_link %.6f) | aa0 %.6f (uni %.6f) | 37B chi0 %.4f c2 %.4f" % (
            fN["chi0"], fN["c2"], res[f"net_{lam0}"]["NN_max"], res[f"net_{lam0}"]["NN_m1_over_37B"], fNa["chi0"], -ea.sum(), faa["chi0"], uni["d2E_da2"], f37["chi0"], f37["c2"]), flush=True)
    for m in (0.2, 0.5):
        uni = dirac_uniform(L, m); rows = [dict(nq=list(n), **dirac_kernel(L, m, n)) for n in QS_AX + QS_DG]
        ea = np.array(rows[0]["e_a"]); fN = fitq(rows, "NN"); fNa = fitq(rows, "Na"); faa = fitq(rows, "aa")
        res[f"dirac_{m}"] = dict(uniform=uni, eps_link=float(ea.sum()), H=rows, fit_NN=fN, fit_Na=fNa, fit_aa=faa,
                                 Na0_rel=fNa["chi0"] / (-ea.sum()) - 1, aa0_rel=faa["chi0"] / uni["d2E_da2"] - 1, NN_max=max(r["NN"] for r in rows),
                                 dEda_vs_minus_eps=uni["dE_da"] / (-ea.sum()) - 1)
        print("dirac", m, "| NN chi0 %.2e c2 %.4f max %.2e | Na0 %.6f (-eps_link %.6f) | aa0 %.6f (uni %.6f) | dE/da %.6f" % (
            fN["chi0"], fN["c2"], res[f"dirac_{m}"]["NN_max"], fNa["chi0"], -ea.sum(), faa["chi0"], uni["d2E_da2"], uni["dE_da"]), flush=True)
    return res


def save(name, obj, t0):
    OUT.mkdir(parents=True, exist_ok=True)
    rec = dict(study="EX-GC50", part=name, label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(), numpy=np.__version__,
               wall_s=time.time() - t0, result=obj)
    (OUT / f"{name}.json").write_text(json.dumps(rec, indent=1, default=float))


if __name__ == "__main__":
    t0 = time.time(); save("kernels", part_kernels(), t0); print("done %.0fs" % (time.time() - t0))
