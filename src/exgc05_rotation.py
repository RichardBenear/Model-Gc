"""EX-GC05 (exploratory), Stage A: does a non-round mode of the Gc1 particle lose stability through an oscillation?

Radial particle (EX-GC01 'ref': G feeds X, k = 0.2), species response times tau_s d_t u_s = F_s (steady state
unchanged). For angular index l = 0..3, leading eigenvalues of J_l psi = lambda T psi found by shift-invert at several
shifts along the imaginary axis. Each mode is tagged by its core share (fraction of |psi|^2 in X, Y, W at r < 10; G-diffusion and ether modes
have small shares).
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ex008_longrange as e  # noqa: E402
import exgc02_inertia as g2  # noqa: E402

OUT = os.path.join(e.ROOT, "results", "EX-GC05")
SHIFTS = [0.0, 0.5j, 1j, 2j, 3j, 4j, 5j, 6j, 7j]


def spectrum(M, u, l, tau, k=20, keep=-1.0):
    J = g2.Jl(M, u, l)
    T = g2.tau_vec(M, tau)
    core = np.tile(M.r < 10, 4)
    xyw = core & (np.arange(4 * M.n) < 3 * M.n)     # X, Y, W components inside r < 10
    found = []
    for sg in SHIFTS:
        try:
            vals, vecs = g2.gen_eigs(J, T, sigma=sg + 1e-4, k=k)
        except Exception as ex:          # noqa: BLE001
            print("eigs failed", l, sg, ex, flush=True)
            continue
        for lam, v in zip(vals, vecs.T):
            if lam.real < keep or lam.imag < -1e-9:
                continue
            if any(abs(lam - f["lam"]) < 1e-7 * max(1, abs(lam)) for f in found):
                continue
            w = np.abs(v) ** 2
            found.append(dict(lam=complex(lam), re=float(lam.real), im=float(lam.imag),
                              core=float(w[xyw].sum() / w.sum())))
    found.sort(key=lambda f: -f["re"])
    for f in found:
        f.pop("lam")
    return found


def scan(tauYs, ls=(0, 1, 2, 3), tag="tauY"):
    M, u = g2.load("ref")
    rows = []
    for tY in tauYs:
        for l in ls:
            sp_ = spectrum(M, u, l, {"Y": tY})
            loc = [m for m in sp_ if m["core"] > 0.5]
            rows.append(dict(tauY=tY, l=l, modes=sp_[:6], core_modes=loc[:6]))
            top = [f"{m['re']:+.4f}{m['im']:+.3f}i c{m['core']:.2f}" for m in loc[:4]]
            print("tauY", tY, "l", l, " | ".join(top), flush=True)
        json.dump(rows, open(os.path.join(OUT, f"scan_{tag}.json"), "w"), indent=1)
    return rows


if __name__ == "__main__" and sys.argv[1:2] != ["type"]:
    os.makedirs(OUT, exist_ok=True)
    tY = [float(x) for x in sys.argv[1:]] or [1, 2, 4, 6, 7, 8, 8.5, 9]
    scan(tY, tag="tauY_" + "_".join(f"{x:g}" for x in tY))


def scan_type(cpl_name, k, tauYs, ls=(0, 1, 2, 3)):
    """Same scan for another coupling / k (EX-GC03 particle builder, continued from the reference particle)."""
    import exgc03_coupling as g3
    M, u = g3.particle(g3.COUPLINGS[cpl_name], dict(k=k, c=1.0, dw=0.1))
    assert M is not None, "particle not found"
    tag = f"{cpl_name}_k{k:g}"
    rows = []
    for tY in tauYs:
        for l in ls:
            sp_ = spectrum(M, u, l, {"Y": tY})
            loc = [m for m in sp_ if m["core"] > 0.5]
            rows.append(dict(type=tag, tauY=tY, l=l, modes=sp_[:6], core_modes=loc[:6]))
            top = [f"{m['re']:+.4f}{m['im']:+.3f}i c{m['core']:.2f}" for m in loc[:4]]
            print(tag, "tauY", tY, "l", l, " | ".join(top), flush=True)
        json.dump(rows, open(os.path.join(OUT, f"scan_{tag}.json"), "w"), indent=1)


if __name__ == "__main__" and sys.argv[1:2] == ["type"]:
    os.makedirs(OUT, exist_ok=True)
    scan_type(sys.argv[2], float(sys.argv[3]), [float(x) for x in sys.argv[4:]])
