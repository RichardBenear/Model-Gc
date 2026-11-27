"""EX-GC48 checks by a separate path. C1: L = 8 dense band-basis two-body Hamiltonian from numerically diagonalised
spinors (no analytic projector, no 5x5 trick) vs exgc48.Pair. C2: real-space Wilson-Dirac Hamiltonian at L = 24,
<H_link> of a positive-band packet vs the k-space Q_a formula (W1, W2). Shared dependency: numpy; C1 compares against
exgc48.Pair, C2 against exgc48.Qa1. EXPLORATORY. Our construction; not LaViolette's Model G, not Chain Theory.
"""
import json, time, platform
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import exgc48 as G

sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.array([[1, 0], [0, -1]], complex)
I2 = np.eye(2); O = np.zeros((2, 2))
A = [np.block([[O, s], [s, O]]) for s in (sx, sy, sz)]; B = np.block([[I2, O], [O, -I2]])


def hk(k, m, r, t):
    return sum(t * np.sin(k[j]) * A[j] for j in range(3)) + (m + t * r * sum(1 - np.cos(k[j]) for j in range(3))) * B


def dense(L, m, U, K, r=1.0, t=1.0):
    ks = 2 * np.pi * np.fft.fftfreq(L); pts = [np.array([a, b, c]) for a in ks for b in ks for c in ks]; N = len(pts)
    D = []; C = []
    for k in pts:
        e1, u1 = np.linalg.eigh(hk(k + K, m, r, t)); e2, u2 = np.linalg.eigh(hk(-k, m, r, t))
        for b1 in (2, 3):
            for b2 in (2, 3):
                D.append(e1[b1] + e2[b2]); C.append(np.kron(u1[:, b1], u2[:, b2]))
    C = np.array(C); H = np.diag(D) - U / N * (C.conj() @ C.T)
    return float(np.linalg.eigvalsh(H)[0])


def c1(L=8, m=0.5, U=8.0, h=1e-5):
    K1 = np.array([2 * np.pi / L, 0, 0]); out = {}
    for name, K in (("K0", np.zeros(3)), ("K1", K1)):
        ed = dense(L, m, U, K); ef = G.Pair(L, m, K=K).solve(U); out[name] = dict(dense=ed, fast=ef, diff=ed - ef)
    qd = (dense(L, m, U, np.zeros(3), t=1 + h) - dense(L, m, U, np.zeros(3), t=1 - h)) / (2 * h)
    qf = (G.Pair(L, m, t=1 + h).solve(U) - G.Pair(L, m, t=1 - h).solve(U)) / (2 * h)
    out["Qa"] = dict(dense=qd, fast=qf, diff=qd - qf); out["max_abs_diff_E"] = max(abs(out["K0"]["diff"]), abs(out["K1"]["diff"]))
    return out


def c2(L=24, m=0.3, r=1.0, sig=3.0, n0=2):
    N = L ** 3; idx = np.arange(N).reshape(L, L, L)
    def shift(j):   # (T_j psi)(x) = psi(x + e_j)
        rows = idx.reshape(-1); cols = np.roll(idx, -1, axis=j).reshape(-1); return sp.csr_matrix((np.ones(N), (rows, cols)), shape=(N, N))
    kin = sp.csr_matrix((4 * N, 4 * N), dtype=complex); wil_hop = sp.csr_matrix((4 * N, 4 * N), dtype=complex)
    for j in range(3):
        T = shift(j); kin = kin + sp.kron(T - T.T, A[j] / 2j); wil_hop = wil_hop + sp.kron(T + T.T, -r / 2 * B)
    site_w = sp.kron(sp.identity(N), 3 * r * B); mass = sp.kron(sp.identity(N), m * B)
    HW1_link = kin + wil_hop + site_w; HW2_link = kin + wil_hop
    k = G.grid(L); k0 = np.array([2 * np.pi * n0 / L, 0, 0]); dk = (k - k0 + np.pi) % (2 * np.pi) - np.pi
    g = np.exp(-np.sum(dk ** 2, -1) * sig ** 2 / 2); chi = np.array([1, 0, 0, 0], complex)
    amps = np.zeros((N, 4), complex)
    for i in range(N):
        e, u = np.linalg.eigh(hk(k[i], m, r, 1.0)); P = u[:, 2:] @ u[:, 2:].conj().T; a = P @ chi; amps[i] = g[i] * a / np.linalg.norm(a)
    psi = np.zeros((L, L, L, 4), complex)
    for al in range(4):
        psi[..., al] = np.fft.ifftn(amps[:, al].reshape(L, L, L)) * np.sqrt(N)   # psi(x) = sum_k a(k) e^{ikx}/sqrt(N)
    v = psi.reshape(N, 4).reshape(-1); nrm = np.vdot(v, v).real
    w = g ** 2 / np.sum(g ** 2); out = {}
    for name, Hl in (("W1", HW1_link), ("W2", HW2_link)):
        rs = float(np.vdot(v, Hl @ v).real / nrm); ks = float(np.sum(w * G.Qa1(k, m, r, name))); out[name] = dict(real=rs, kspace=ks, diff=rs - ks)
    Etot = float(np.vdot(v, (HW1_link + mass) @ v).real / nrm); out["E_real"] = Etot; out["E_k"] = float(np.sum(w * G.parts(k, m, r)[2]))
    return out


if __name__ == "__main__":
    t0 = time.time(); r1 = c1(); print("C1", r1, flush=True); r2 = c2(); print("C2", r2, flush=True)
    rec = dict(study="EX-GC48", part="checks", label="EXPLORATORY", platform=platform.platform(), python=platform.python_version(),
               wall_s=time.time() - t0, result=dict(C1=r1, C2=r2))
    (G.OUT / "checks.json").write_text(json.dumps(rec, indent=1, default=float))
