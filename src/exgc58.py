"""EX-GC58: the link shift as local events. Link event: a real directed transfer (rotation of the pair (psi_x, psi_{x+a}) by
theta = beta dt / 2, spin- and mass-blind) = exp(-i H_beta,link dt). Plaquette event: spinor rotation exp(i (dt/4)(C beta).Sigma) at each
site, C = 57's site curl (mean of the 4 plaquettes touching the site), i.e. H_pl = -1/4 sum s.(C beta), the Belinfante term. Source = the
events' conjugate count, j = P + 1/4 C^T s. Rule labels: P (links only, 55-57) and PS (links + plaquettes).
58A rule tests; 58B drift universality with event stepping; 58C spin drag in 57's ring field; 58D spin sources the shift;
58E energy bookkeeping (55B setup); 58F co-moving pair (55C setup). Plan: spec/EX-GC58_plan.md.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory.
"""
import json, sys, time, platform
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G
import exgc52 as S
import exgc53 as T53
import exgc55 as X55
import exgc56 as X56
import exgc57 as E57

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC58"
SIG = [-1j * G.AL[1] @ G.AL[2], -1j * G.AL[2] @ G.AL[0], -1j * G.AL[0] @ G.AL[1]]
PAIRS = [(1, 2), (2, 0), (0, 1)]                                                 # (a, b) with (a, b, comp) right-handed


# ------------------------------------------------------------------ plaquette operator and its exact transpose
def _A(q, a, b):
    return 0.25 * (q + np.roll(q, 1, a) + np.roll(q, 1, b) + np.roll(np.roll(q, 1, a), 1, b))


def _AT(q, a, b):
    return 0.25 * (q + np.roll(q, -1, a) + np.roll(q, -1, b) + np.roll(np.roll(q, -1, a), -1, b))


def curl_op(B):
    out = np.zeros((3,) + B[0].shape)
    for c, (a, b) in enumerate(PAIRS):
        Q = B[a] + np.roll(B[b], -1, a) - np.roll(B[a], -1, b) - B[b]; out[c] = _A(Q, a, b)
    return out


def curl_T(s):
    out = np.zeros((3,) + s[0].shape)
    for c, (a, b) in enumerate(PAIRS):
        q = _AT(s[c], a, b); out[a] += q - np.roll(q, 1, b); out[b] += np.roll(q, 1, a) - q
    return out


def sdens(psi):
    return np.stack([np.real(np.einsum("...i,ij,...j->...", psi.conj(), Sm, psi)) for Sm in SIG])


def Hpl(psi, C):
    out = np.zeros_like(psi)
    for c in range(3):
        out -= 0.25 * C[c][..., None] * (psi @ SIG[c].T)
    return out


def coupling58(psi, N, m, B, rule):
    if rule == "PS":
        P = T53.link_P(psi); j = P + 0.25 * curl_T(sdens(psi)); return T53.Hb(psi, B) + Hpl(psi, curl_op(B)), j, j
    return _coupling55(psi, N, m, B, rule)


_coupling55 = X55.coupling; X55.coupling = coupling58                            # 55's evolve and Etot now accept rule "PS"


# ------------------------------------------------------------------ events
def link_event(psi, a, x, th):
    y = list(x); y[a] += 1; y = tuple(y); px, py = psi[x].copy(), psi[y].copy(); c, s = np.cos(th), np.sin(th)
    psi[x] = c * px + s * py; psi[y] = -s * px + c * py


def class_sweep(psi, B, dt, reverse=False):
    """all link events, applied in 6 classes of disjoint links (axis, parity of x_a); reverse=True for the mirrored half step."""
    n = psi.shape[0]; cl = [(a, p) for a in range(3) for p in (0, 1)]
    for a, p in (cl[::-1] if reverse else cl):
        idx = np.arange(p, n - 1, 2); sx = [slice(None)] * 3; sy = [slice(None)] * 3; sx[a] = idx; sy[a] = idx + 1; sx, sy = tuple(sx), tuple(sy)
        th = (B[a][sx] * dt / 2)[..., None]; px, py = psi[sx].copy(), psi[sy].copy(); c, s = np.cos(th), np.sin(th)
        psi[sx] = c * px + s * py; psi[sy] = -s * px + c * py
    return psi


def spin_events(psi, C, dt):
    ph = (dt / 4) * C; mag = np.sqrt(np.sum(ph ** 2, 0)); nz = mag > 0; out = psi.copy()
    nS = sum(np.where(nz, ph[c] / np.where(nz, mag, 1), 0)[..., None] * (psi @ SIG[c].T) for c in range(3))
    return np.cos(mag)[..., None] * psi + 1j * np.sin(mag)[..., None] * nS


def expm_apply(Hf, psi, dt, tol=1e-18):
    out = psi.copy(); term = psi.copy(); k = 0
    while True:
        k += 1; term = -1j * dt / k * Hf(term); out = out + term
        if np.linalg.norm(term) < tol or k > 40:
            return out


def part_A(n=13, dt=0.1):
    rng = np.random.default_rng(58); res = dict(rows=[])
    psi0 = rng.standard_normal((n, n, n, 4)) + 1j * rng.standard_normal((n, n, n, 4)); psi0 /= np.linalg.norm(psi0)
    B1 = rng.standard_normal((3, n, n, n)) * 1e-3; links = [(a,) + tuple(int(v) for v in x) for a in range(3) for x in np.ndindex(n, n, n) if x[a] < n - 1]
    orders = [list(range(len(links)))] + [list(rng.permutation(len(links))) for _ in range(4)]
    for s in (1.0, 2.0, 4.0):
        B = s * B1; C = curl_op(B); row = dict(s=s)
        for mode in ("P", "PS"):
            Hf = (lambda p: T53.Hb(p, B)) if mode == "P" else (lambda p: T53.Hb(p, B) + Hpl(p, C))
            ex = expm_apply(Hf, psi0, dt); outs = []; norms = []
            for o in orders:
                p = psi0.copy()
                for i in o:
                    a, x = links[i][0], links[i][1:]; link_event(p, a, x, B[a][x] * dt / 2)
                if mode == "PS":
                    p = spin_events(p, C, dt)
                outs.append(p); norms.append(abs(np.linalg.norm(p) - 1))
            first = float(np.linalg.norm(ex - psi0)); err = float(np.linalg.norm(outs[0] - ex))
            spread = float(max(np.linalg.norm(outs[i] - outs[j]) for i in range(1, 5) for j in range(i + 1, 5)))
            row[mode] = dict(first_order=first, err_sweep=err, spread_orders=spread, norm_dev=float(max(norms)))
        res["rows"].append(row); print("58A", json.dumps(row), flush=True)
    sl = lambda key, mode: float(np.polyfit(np.log([r["s"] for r in res["rows"]]), np.log([r[mode][key] for r in res["rows"]]), 1)[0])
    res["slopes"] = {mode: dict(err=sl("err_sweep", mode), spread=sl("spread_orders", mode)) for mode in ("P", "PS")}
    # C1: beta = 0 events are the identity; C2: H_pl Hermitian and <H_pl> = -sum beta . 1/4 C^T s
    p = psi0.copy(); class_sweep(p, 0 * B1, dt); p = spin_events(p, curl_op(0 * B1), dt); res["C1"] = float(np.linalg.norm(p - psi0))
    ph = rng.standard_normal(psi0.shape) + 1j * rng.standard_normal(psi0.shape); C = curl_op(B1)
    herm = abs(np.vdot(ph, Hpl(psi0, C)) - np.vdot(Hpl(ph, C), psi0)) / abs(np.vdot(ph, Hpl(psi0, C)))
    ex = float(np.real(np.vdot(psi0, Hpl(psi0, C)))); cj = -float(np.sum(B1 * 0.25 * curl_T(sdens(psi0))))
    res["C2"] = dict(herm=float(herm), expect=ex, conj=cj, rel=abs(ex - cj) / abs(cj)); print("58A slopes", res["slopes"], "C1", res["C1"], "C2", res["C2"], flush=True)
    return res


# ------------------------------------------------------------------ 58B
B_CASES = [(m, k0, r, T53.UP) for m in (0.3, 1.2) for r in (0.0, 1.0) for k0 in (0.0, 0.3)] + [(0.6, 0.15, 0.0, T53.UP), (0.6, 0.15, 0.0, T53.DOWN)]


def part_B(i, n=48, T=10.0, dt=0.1, beta=1e-3):
    m, k0, r, chi = B_CASES[i]; c = n // 2; one = np.ones((n, n, n)); psi, K, w = S.packet_k(n, m, r, (c, c, c), S.SIG, (k0, 0, 0), chi=chi)
    B = np.zeros((3, n, n, n)); B[0] = beta; Hm = lambda p: G.Hm(p, one, m, r)
    def rk4(p):
        k1 = -1j * Hm(p); k2 = -1j * Hm(p + dt / 2 * k1); k3 = -1j * Hm(p + dt / 2 * k2); k4 = -1j * Hm(p + dt * k3); return p + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
    p0 = psi.copy(); pe = psi.copy(); t = 0.0
    while t < T - 1e-9:
        p0 = rk4(p0); pe = class_sweep(pe, B, dt / 2); pe = rk4(pe); pe = class_sweep(pe, B, dt / 2, reverse=True); t += dt
    drift_ev = (X56.centroid(pe) - X56.centroid(p0)) / (beta * T); ham = X56.matter_case(m, k0, r, "P", chi=chi, n=n, T=T, dt=dt, beta=beta)
    row = dict(m=m, k0=k0, r=r, chi=list(chi), drift_event=drift_ev, drift_ham=ham["drift"], pred=ham["pred_P"], ratio=drift_ev / ham["pred_P"], rel_ham=drift_ev / ham["drift"] - 1)
    print("58B", json.dumps(row), flush=True); return row


# ------------------------------------------------------------------ 58C
def part_C(n=129, u=0.01):
    c = n // 2; q0, _ = T53.ring_sources(n, 0); kap = 4 * np.pi * 12 * u / float(q0.sum()); res = dict(positions={})
    F = {}
    for ell in (2, -2):
        q, P = T53.ring_sources(n, ell); B = T53.shift_static(P, kap, T53.angmom(P)); F[ell] = (B, curl_op(B))
    spins = {sig: {nm: E57.spin_packet(sig, chi) for nm, chi in (("up", T53.UP), ("down", T53.DOWN))} for sig in (3.0, 5.0)}
    Jz = {}; Sz = {}
    for sig, d in spins.items():
        for nm, p in d.items():
            Sz[(sig, nm)] = float(np.sum(sdens(p)[2])); Jz[(sig, nm)] = float(T53.angmom(T53.link_P(p))[2]) + 0.5 * Sz[(sig, nm)]
    res["spin_J"] = {f"{s:g}_{nm}": dict(Sz=Sz[(s, nm)], Jz=Jz[(s, nm)]) for (s, nm) in Sz}
    probes = {lp: [E57.probe_ring(lp, chi, (0, 0, 1)) for chi in (T53.UP, T53.DOWN)] for lp in (1, -1)}
    Lz = {lp: np.mean([T53.angmom(T53.link_P(p))[2] for p in probes[lp]]) for lp in (1, -1)}
    def E(p, Bs, Cs, rule):
        e = -float(np.sum(Bs * T53.link_P(p)))
        return e - 0.25 * float(np.sum(sdens(p) * Cs)) if rule == "PS" else e
    for d in (24, 32):
        for key, x0 in (("axis", (c, c, c + d)), ("plane", (c + d, c, c))):
            row = {}
            for ell in (2, -2):
                B, C = F[ell]; Bs = E57.sub(B, x0); Cs = E57.sub(C, x0); oc = -0.5 * E57.curl_site(B, x0)[2]; o = dict(Omega_curl=oc)
                for rule in ("P", "PS"):
                    Eo = {lp: np.mean([E(p, Bs, Cs, rule) for p in probes[lp]]) for lp in (1, -1)}; o[f"orb_{rule}"] = (Eo[1] - Eo[-1]) / (Lz[1] - Lz[-1]) / oc
                    if key == "axis":
                        for sig in spins:
                            dE = E(spins[sig]["up"], Bs, Cs, rule) - E(spins[sig]["down"], Bs, Cs, rule)
                            o[f"spin_{sig:g}_{rule}_dEdJ"] = dE / (Jz[(sig, "up")] - Jz[(sig, "down")]) / oc
                            o[f"spin_{sig:g}_{rule}_57def"] = 2 * dE / (Sz[(sig, "up")] - Sz[(sig, "down")]) / oc
                row[str(ell)] = o
            res["positions"][f"{key}_{d}"] = row; print("58C", key, d, json.dumps(row["2"]), flush=True)
    # gradient control: the PS spin energy in beta = grad chi (amplitude as 57 at axis 24)
    x0 = (c, c, c + 24); Bg = E57.grad_field(n, x0, 0.015197109918616345); Cg = curl_op(Bg); Cgs = E57.sub(Cg, x0)
    eg = -0.25 * float(np.sum(sdens(spins[3.0]["up"]) * Cgs) - np.sum(sdens(spins[3.0]["down"]) * Cgs))
    ering = -0.25 * float(np.sum(sdens(spins[3.0]["up"]) * E57.sub(F[2][1], x0)) - np.sum(sdens(spins[3.0]["down"]) * E57.sub(F[2][1], x0)))
    res["grad"] = dict(spin_term=eg, ring_spin_term=ering, rel=abs(eg / ering), max_curl_grad=float(np.abs(Cg).max()), max_grad=float(np.abs(Bg).max()))
    print("58C grad", res["grad"], flush=True); return res


# ------------------------------------------------------------------ 58D
def part_D(n=97, sig=3.0, kap=1.0):
    c = n // 2; psi, _, _ = S.packet_k(n, S.M, 0.0, (c, c, c), sig, (0, 0, 0), chi=T53.UP); P = T53.link_P(psi); s = sdens(psi)
    L = T53.angmom(P); Sv = 0.5 * s.sum(axis=(1, 2, 3)); res = dict(L=L.tolist(), S=Sv.tolist()); Gk = kap / (4 * np.pi); cv = np.array([c, c, c], float)
    for rule, j in (("P", P), ("PS", P + 0.25 * curl_T(s))):
        J = T53.angmom(j); B = T53.shift_static(j, kap, J); row = dict(J=J.tolist(), J_over_L=J[2] / L[2], J_over_LS=J[2] / (L[2] + Sv[2]), loops={})
        for d in (20, 28):
            st = T53.square_loop((c, c, c + d), 12); _, bl = T53.loop_fields(st, np.ones((n, n, n)), B); ro = 2 * float(np.sum(bl)); dp = T53.loop_dipole(st, J, Gk, cv)
            row["loops"][str(d)] = dict(readout=ro, dipole=dp, ratio=ro / dp)
        res[rule] = row; print("58D", rule, json.dumps(row), flush=True)
    return res


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC58_{name}.json").write_text(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    if job == "A":
        r = part_A()
    elif job.startswith("B"):
        r = part_B(int(job[1:]))
    elif job == "C":
        r = part_C()
    elif job == "D":
        r = part_D()
    elif job.startswith("E_"):                      # E_0.1, E_0.05
        r = X55.part_B("PS", dt=float(job[2:]))
    else:                                           # F_PS_0.2, F_P_0.2
        p = job.split("_"); r = X55.run_pair(float(p[2]), p[1])
    save(job, r, t0); print(job, "done %.0fs" % (time.time() - t0), flush=True)
