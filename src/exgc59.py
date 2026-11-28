"""EX-GC59: spin source to spin probe, end to end. A resting polarised packet sources N and the link shift (58's PS rule:
j = P + 1/4 C^T s; links-only rule P: j = P); static shift with ZERO boundary values (no J handed to the solver); a transverse probe
spin evolves in the frozen fields (48^3 box around the probe); precession angle phi = atan2(<Sigma_y>, <Sigma_x>), odd in the source spin.
Prediction: free-space linearised GR from the source's own count, averaged over the probe's initial density. Plan: spec/EX-GC59_plan.md.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np
import numba

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc55_posthoc2 as SY
import exgc57 as E57
import exgc58 as X

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC59"
M = 0.6; SIGMA = 3.0; LAM = 4.0; NB = 48; DT = 0.1; T = 10.0
UP, DOWN = T53.UP, T53.DOWN; CHI_X = (1 / np.sqrt(2), 1 / np.sqrt(2), 0, 0)


# ------------------------------------------------------------------ source and fields
def source(n, chi):
    c = n // 2; psi, _, _ = S.packet_k(n, M, 0.0, (c, c, c), SIGMA, (0, 0, 0), chi=chi); return psi


def source_terms(n):
    one = np.ones((n, n, n)); out = {}
    for nm, chi in (("up", UP), ("down", DOWN)):
        p = source(n, chi); P = T53.link_P(p); s = X.sdens(p)
        out[nm] = dict(q=G.charge(p, one, M, 0.0), P=P, PS=P + 0.25 * X.curl_T(s), L=float(T53.angmom(P)[2]), Sz=0.5 * float(s[2].sum()))
    q = out["up"]["q"]; kap = 4 * np.pi * 12 * 0.01 / float(q.sum()); c = n // 2
    N, info = G.static([], [], 0.0, kap, n, (c, c, c), qfix=q); Nd, _ = G.static([], [], 0.0, kap, n, (c, c, c), qfix=out["down"]["q"])
    return out, kap, N, Nd


def shift(j, kap, bc):
    if bc == "dip":
        return T53.shift_static(j, kap, T53.angmom(j))
    return np.stack([T53.poisson(LAM * kap * j[a], np.zeros_like(j[a])) for a in range(3)])


def omega_gr(j, kap, pts, n):
    """Omega_GR,z at points pts (k, 3) (coordinates relative to the box centre), free space."""
    LP = T53.link_pos(n); out = np.zeros(len(pts)); thr = 1e-10 * np.abs(j).max()
    for a in range(3):
        msk = np.abs(j[a]) > thr; r = np.stack([LP[a][i][msk] for i in range(3)], 1); v = j[a][msk]
        for i0 in range(0, len(pts), 256):
            d = pts[i0:i0 + 256, None, :] - r[None]; dd = np.linalg.norm(d, axis=2) ** 3
            cz = d[..., 0] * (v if a == 1 else 0) - d[..., 1] * (v if a == 0 else 0)  # z-component of d x (v e_a)
            out[i0:i0 + 256] += np.sum(cz / dd, 1)
    return -(LAM * kap / (8 * np.pi)) * out


def beta_gr_loop(j, kap, steps, n):
    LP = T53.link_pos(n); c = n // 2; tot = 0.0; thr = 1e-10 * np.abs(j).max()
    srcs = []
    for a in range(3):
        msk = np.abs(j[a]) > thr; srcs.append((np.stack([LP[a][i][msk] for i in range(3)], 1), j[a][msk]))
    for (x, ax, sg) in steps:
        p = np.array(x, float) - c; p[ax] += 0.5 * sg; r, v = srcs[ax]; dd = np.linalg.norm(p - r, axis=1); ok = dd > 1e-9  # drop the singular self-term (loop link = source link)
        tot += sg * (-(LAM * kap / (4 * np.pi)) * np.sum(v[ok] / dd[ok]))
    return tot


# ------------------------------------------------------------------ probe evolution
def probe():
    c = NB // 2; psi, _, _ = S.packet_k(NB, M, 0.0, (c, c, c), SIGMA, (0, 0, 0), chi=CHI_X); return psi


def spin_xy(psi):
    sx = float(np.real(np.vdot(psi, psi @ X.SIG[0].T))); sy = float(np.real(np.vdot(psi, psi @ X.SIG[1].T))); return sx, sy


def cut_site(N, x0):
    h = NB // 2; return N[tuple(slice(x0[i] - h, x0[i] + h) for i in range(3))]


def run_ham(Ns, Bs, rule):
    psi = probe(); C = X.curl_op(Bs) if (rule == "PS" and Bs is not None) else None
    def H(p):
        out = G.Hm(p, Ns, M, 0.0)
        if Bs is not None:
            out = out + T53.Hb(p, Bs)
            if rule == "PS":
                out = out + X.Hpl(p, C)
        return out
    phis, norms = [], []; t = 0.0
    def rec(p):
        sx, sy = spin_xy(p); phis.append(np.arctan2(sy, sx)); norms.append(abs(np.linalg.norm(p) - 1))
    rec(psi)
    while t < T - 1e-9:
        k1 = -1j * H(psi); k2 = -1j * H(psi + DT / 2 * k1); k3 = -1j * H(psi + DT / 2 * k2); k4 = -1j * H(psi + DT * k3)
        psi = psi + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += DT; rec(psi)
    return dict(phi=phis, norm_dev=float(max(norms)))


@numba.njit(cache=True)
def _events(psi, kind, ax, ix, iy, iz, ang, sig):
    """apply events in the given order. kind 0: link (ax, x) rotation by ang; kind 1: site spin rotation by vector ang3 = sig-weighted."""
    for e in range(kind.shape[0]):
        i, j, k = ix[e], iy[e], iz[e]
        if kind[e] == 0:
            a = ax[e]; i2, j2, k2 = i + (a == 0), j + (a == 1), k + (a == 2); th = ang[e, 0]; c = np.cos(th); s = np.sin(th)
            for m in range(4):
                p1 = psi[i, j, k, m]; p2 = psi[i2, j2, k2, m]; psi[i, j, k, m] = c * p1 + s * p2; psi[i2, j2, k2, m] = -s * p1 + c * p2
        else:
            ph0, ph1, ph2 = ang[e, 0], ang[e, 1], ang[e, 2]; mag = np.sqrt(ph0 * ph0 + ph1 * ph1 + ph2 * ph2)
            if mag == 0:
                continue
            v = psi[i, j, k].copy(); nv = np.zeros(4, dtype=np.complex128)
            for r in range(4):
                acc = 0j
                for q in range(4):
                    acc += (ph0 * sig[0, r, q] + ph1 * sig[1, r, q] + ph2 * sig[2, r, q]) / mag * v[q]
                nv[r] = acc
            for r in range(4):
                psi[i, j, k, r] = np.cos(mag) * v[r] + 1j * np.sin(mag) * nv[r]


def run_async(Ns, Bs, seed):
    rng = np.random.default_rng(seed); psi = probe(); h = DT / 2; C = X.curl_op(Bs); n = NB
    la, li, lj, lk = [], [], [], []
    for a in range(3):
        g = np.indices((n, n, n)).reshape(3, -1).T; g = g[g[:, a] < n - 1]; la.append(np.full(len(g), a)); li.append(g[:, 0]); lj.append(g[:, 1]); lk.append(g[:, 2])
    la = np.concatenate(la); li = np.concatenate(li); lj = np.concatenate(lj); lk = np.concatenate(lk); nl = len(la)
    lang = np.zeros((nl, 3)); lang[:, 0] = Bs[la, li, lj, lk] * h / 2
    s_ = np.indices((n, n, n)).reshape(3, -1).T; sang = (h / 4) * C[:, s_[:, 0], s_[:, 1], s_[:, 2]].T
    kind = np.concatenate([np.zeros(nl, np.int64), np.ones(len(s_), np.int64)]); ax = np.concatenate([la, np.zeros(len(s_), np.int64)]).astype(np.int64)
    ix = np.concatenate([li, s_[:, 0]]).astype(np.int64); iy = np.concatenate([lj, s_[:, 1]]).astype(np.int64); iz = np.concatenate([lk, s_[:, 2]]).astype(np.int64)
    ang = np.concatenate([lang, sang]); sig = np.array(X.SIG, np.complex128); Hm = lambda p: G.Hm(p, Ns, M, 0.0)
    phis, norms = [], []
    def rec(p):
        sx, sy = spin_xy(p); phis.append(np.arctan2(sy, sx)); norms.append(abs(np.linalg.norm(p) - 1))
    rec(psi); t = 0.0; psi = np.ascontiguousarray(psi)
    while t < T - 1e-9:
        o = rng.permutation(len(kind)); _events(psi, kind[o], ax[o], ix[o], iy[o], iz[o], ang[o], sig)
        k1 = -1j * Hm(psi); k2 = -1j * Hm(psi + DT / 2 * k1); k3 = -1j * Hm(psi + DT / 2 * k2); k4 = -1j * Hm(psi + DT * k3)
        psi = np.ascontiguousarray(psi + DT / 6 * (k1 + 2 * k2 + 2 * k3 + k4)); o = o[::-1]; _events(psi, kind[o], ax[o], ix[o], iy[o], iz[o], ang[o], sig); t += DT; rec(psi)
    return dict(phi=phis, norm_dev=float(max(norms)))


def slope(phi):
    t = DT * np.arange(len(phi)); A = np.stack([np.ones_like(t), t], 1); return float(np.linalg.lstsq(A, np.asarray(phi), rcond=None)[0][1])


def pos(n, key):
    c = n // 2; kind, d = key.split("_"); d = int(d); return (c, c, c + d) if kind == "axis" else (c + d, c, c)


def omega_pred(j, kap, n, x0):
    ps = probe(); w = np.sum(np.abs(ps) ** 2, -1); msk = w > 1e-6 * w.max(); idx = np.argwhere(msk)
    pts = idx - NB // 2 + (np.array(x0) - n // 2); om = omega_gr(j, kap, pts.astype(float), n); ww = w[msk] / w[msk].sum()
    return float(np.sum(ww * om)), float(om[np.argmin(np.sum((idx - NB // 2) ** 2, 1))])


# ------------------------------------------------------------------ parts
def part_main(n=97, bc="zero", keys=("axis_12", "axis_16", "plane_16"), full=True):
    t0 = time.time(); src, kap, N, Nd = source_terms(n); res = dict(n=n, bc=bc, kappa=kap, source={k: dict(L=v["L"], Sz=v["Sz"], J_PS=float(T53.angmom(v["PS"])[2])) for k, v in src.items()}, runs={})
    print("59 fields n", n, "kap %.4f" % kap, res["source"], "%.0fs" % (time.time() - t0), flush=True)
    rules = ("PS", "P") if full else ("PS",)
    B = {(r, s): shift(src[s][r], kap, bc) for r in rules for s in ("up", "down")}
    if full:
        B[("PS", "unpol")] = 0.5 * (B[("PS", "up")] + B[("PS", "down")])
    for key in keys:
        x0 = pos(n, key); Ns = cut_site(N, x0); row = {}
        pr, pc = omega_pred(src["up"]["PS"], kap, n, x0); r = np.array(x0, float) - n // 2; d = np.linalg.norm(r); Jz = res["source"]["up"]["J_PS"]; Gk = kap / (4 * np.pi)
        row["Omega_pred"] = pr; row["Omega_pred_centre"] = pc; row["Omega_dipole"] = Gk * Jz * (2 if key.startswith("axis") else -1) / d ** 3
        for s in ("up", "down"):
            row[f"PS_{s}"] = run_ham(Ns, E57.sub(B[("PS", s)], x0), "PS")
        row["Omega_meas"] = slope((np.array(row["PS_up"]["phi"]) - np.array(row["PS_down"]["phi"])) / 2)
        row["ratio"] = row["Omega_meas"] / pr
        if full and key == "axis_16":
            row["zero"] = run_ham(Ns, None, "PS"); row["zero_Ndown"] = run_ham(cut_site(Nd, x0), None, "PS")
            row["unpol"] = run_ham(Ns, E57.sub(B[("PS", "unpol")], x0), "PS")
            for s in ("up", "down"):
                row[f"P_{s}"] = run_ham(Ns, E57.sub(B[("P", s)], x0), "P")
            row["Omega_meas_P"] = slope((np.array(row["P_up"]["phi"]) - np.array(row["P_down"]["phi"])) / 2)
            Bsub = E57.sub(B[("PS", "up")], x0); amp = float(np.abs(Bsub[:, 12:36, 12:36, 12:36]).max()) * 8.0 * np.exp(0.5)
            Bg = E57.grad_field(n, x0, amp); row["grad"] = run_ham(Ns, E57.sub(Bg, x0), "PS"); row["grad_amp"] = amp; del Bg
            # light loop on the same region, calibrated first
            one = np.ones(64); bu = 1e-3 * one; f = (SY.wave(one, bu, +1, -1, width=4.0) - SY.wave(one, bu, -1, -1, width=4.0)) / (2 * bu.sum()); row["light_f"] = f
            st = T53.square_loop(x0, 8); dts = {}
            for s in ("up", "down"):
                Nl, bl = T53.loop_fields(st, N if s == "up" else Nd, B[("PS", s)]); dts[s] = SY.wave(Nl, bl, +1, -1, width=4.0) - SY.wave(Nl, bl, -1, -1, width=4.0)
            dodd = (dts["up"] - dts["down"]) / 2; A = 16.0 ** 2; row["light"] = dict(dt=dts, Omega_light=-(dodd / f) / (4 * A), Omega_pred_loop=-beta_gr_loop(src["up"]["PS"], kap, st, n) / (2 * A))
            row["light"]["ratio"] = row["light"]["Omega_light"] / row["light"]["Omega_pred_loop"]
        res["runs"][key] = row
        print("59", n, bc, key, "pred %.4e (centre %.4e, dipole %.4e) meas %.4e ratio %.5f" % (pr, pc, row["Omega_dipole"], row["Omega_meas"], row["ratio"]),
              ("P %.3e light %.5f f %.4f" % (row["Omega_meas_P"], row["light"]["ratio"], row["light_f"])) if "light" in row else "", "%.0fs" % (time.time() - t0), flush=True)
    return res


def part_async(n=97, key="axis_16"):
    t0 = time.time(); src, kap, N, Nd = source_terms(n); x0 = pos(n, key); Ns = cut_site(N, x0); res = dict(runs={})
    Bs = {s: E57.sub(shift(src[s]["PS"], kap, "zero"), x0) for s in ("up", "down")}
    for seed in (1, 2):
        r = {s: run_async(Ns, Bs[s], 5900 + 10 * seed + i) for i, s in enumerate(("up", "down"))}
        r["Omega_meas"] = slope((np.array(r["up"]["phi"]) - np.array(r["down"]["phi"])) / 2); res["runs"][str(seed)] = r
        print("59 async seed", seed, "Omega %.6e" % r["Omega_meas"], "norm", r["up"]["norm_dev"], "%.0fs" % (time.time() - t0), flush=True)
    return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC59_{name}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    jobs = {"main": lambda: part_main(), "n129": lambda: part_main(129, "zero", ("axis_16",), full=False), "dip": lambda: part_main(97, "dip", ("axis_16",), full=False), "async": part_async}
    save(job, jobs[job](), t0); print("done", job, "%.0fs" % (time.time() - t0), flush=True)
