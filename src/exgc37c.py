"""EX-GC37C: moving sink. Comoving steady carrier field, diffusing vs wave-like, on a periodic Y1 lattice.
Plan: spec/EX-GC37C_plan.md (amendment A1). Our construction, not LaViolette's Model G."""
import json, time
from pathlib import Path
import numpy as np
import scipy.fft as fft

n, D, core, S = 160, 1.0, 1.5, 2.0          # GM = S/(4 pi D) = 0.159
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC37C"


def setup():
    x = np.arange(n) - n // 2
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    s = np.exp(-(X ** 2 + Y ** 2 + Z ** 2) / (2 * core ** 2)); s *= S / s.sum()
    k = 2 * np.pi * np.fft.fftfreq(n)
    kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
    K2 = (2 * np.sin(kx / 2)) ** 2 + (2 * np.sin(ky / 2)) ** 2 + (2 * np.sin(kz / 2)) ** 2
    return s, fft.fftn(fft.ifftshift(s)), kx, K2


def field(sh, kx, K2, v, kind):
    den = (-D * K2 + 1j * v * kx) if kind == "diffusion" else (-D * (K2 - v ** 2 * kx ** 2))
    zero = np.abs(den) == 0
    den[zero] = 1
    c = sh / den
    c[zero] = 0
    return c


def readouts(s, c, kx, v):
    f = fft.fftshift(np.real(fft.ifftn(c)))
    gx = fft.fftshift(np.real(fft.ifftn(1j * kx * c)))
    o = n // 2
    dN = lambda i, j, l: f[o + i, o + j, o + l] - f[o + i + 30 if False else o, o + 60, o + 60] * 0  # deficit relative to box mean (mean is 0)
    ahead, behind, side = f[o + 10, o, o], f[o - 10, o, o], f[o, o + 10, o]
    a40, b40 = f[o + 40, o, o], f[o - 40, o, o]
    self_grad = float(np.sum(s * gx) / np.sum(s))
    g = 1 / np.sqrt(1 - v ** 2)
    return dict(v=v, ahead10=float(ahead), behind10=float(behind), side10=float(side), fore_aft=float(ahead / behind),
                side_axial=float(side / ahead), gamma=float(g), ahead40=float(a40), behind40=float(b40), self_gradient=self_grad)


if __name__ == "__main__":
    t0 = time.process_time()
    s, sh, kx, K2 = setup()
    res = {"diffusion": [], "wave": []}
    for kind in res:
        for v in (0.0, 0.05, 0.1, 0.2, 0.4):
            res[kind].append(readouts(s, field(sh, kx, K2, v, kind), kx, v))
    for kind in res:
        r0 = res[kind][0]
        for r in res[kind]:
            r["ahead40_rel_static"] = r["ahead40"] / r0["ahead40"]
            r["behind40_rel_static"] = r["behind40"] / r0["behind40"]
    # separate check: continuum boosted-Coulomb prediction for the wave case, ahead/side at d = 10 (point-source limit)
    res["continuum_wave_side_axial"] = "gamma (boosted 1/r)"
    res["cpu_seconds"] = time.process_time() - t0
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "results.json").write_text(json.dumps(res, indent=2) + "\n")
    for kind in ("diffusion", "wave"):
        for r in res[kind]:
            print(kind, {k: round(v, 5) for k, v in r.items()})
    print(res["cpu_seconds"])
