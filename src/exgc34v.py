"""EX-GC34V: quadratic source check for one Y1 harmonic interaction."""
import json
import platform
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import exgc34r as R

P = R.P
N, O = P.N, P.O
ETA = N.ETA
IDENTITY = np.eye(4)
L = 8.0


def lower_derivative(k, w):
    return 1j * np.r_[-w, k]


def connection(metric, derivative, inverse):
    # axes are upper index, first lower index, second lower index
    return np.einsum("rs,smn->rmn", inverse,
                     (np.einsum("m,sn->smn", derivative, metric)
                      + np.einsum("n,sm->smn", derivative, metric)
                      - np.einsum("s,mn->smn", derivative, metric)) / 2)


def ricci(gamma, derivative):
    return np.einsum("r,rmn->mn", derivative, gamma) - np.einsum("n,rmr->mn", derivative, gamma)


def einstein(mixed_ricci):
    return mixed_ricci - IDENTITY * np.trace(mixed_ricci) / 2


def stress(u, phi):
    return ETA @ np.outer(u, u) - IDENTITY * (u @ ETA @ u + phi**2) / 2


def source_from_grid(phi, vel, h, k_product):
    sn, sb, sd, so = O.msources2(phi, vel, h)
    coeff = np.array([sn[0, 0, 0], *sb[:, 0, 0, 0], *sd[:, 0, 0, 0], *so[:, 0, 0, 0]])
    coeff *= np.exp(-1j * np.asarray(N.OFFS) @ (k_product * h)) / h**3
    up = np.zeros((4, 4), complex)
    up[0, 0] = -coeff[0]
    up[0, 1:] = up[1:, 0] = coeff[1:4]
    for i, (a, b) in enumerate(N.PAIRS):
        up[a+1, b+1] = up[b+1, a+1] = coeff[4+i] * (2 if a == b else 1)
    return up @ ETA


def run(direction, h=None):
    k = np.asarray(direction) * (2 * np.pi / L)

    def K(j):
        return j*k if h is None else 2*np.sin(j*k*h/2)/h

    w = np.sqrt(1 + K(1) @ K(1))
    d1, d2, d3, d4 = [lower_derivative(K(j), j*w) for j in range(1, 5)]
    if h is None:
        t0 = stress(d1, 1)
    else:
        axis = np.arange(round(L/h)) * h
        phase = np.exp(1j * (k[0]*axis[:, None, None] + k[1]*axis[None, :, None] + k[2]*axis[None, None, :]))
        t0 = source_from_grid(phase, -1j*w*phase, h, 2*k)

    kernel = P.kernel_batch(K(2)[None], 2*w)[0]
    forcing = N.PM.T @ (t0 @ ETA).ravel()
    x = np.zeros(10, complex)
    x[4:] = np.linalg.solve(kernel[4:, 4:], forcing[4:])
    first_resid = np.linalg.norm(kernel @ x - forcing) / np.linalg.norm(forcing)
    metric = (N.PM @ x).reshape(4, 4)
    inv1 = -ETA @ metric @ ETA
    gamma1 = connection(metric, d2, ETA)
    gamma2 = connection(metric, d2, inv1)
    ricci1 = ricci(gamma1, d2)
    ricci2 = (ricci(gamma2, d4) + np.einsum("rrs,smn->mn", gamma1, gamma1)
              - np.einsum("rns,smr->mn", gamma1, gamma1))
    g1 = einstein(ETA @ ricci1)
    g2 = einstein(ETA @ ricci2 + inv1 @ ricci1)

    trace = np.trace(metric) / 2
    if h is None:
        gradient_force = d3[1:] @ ((trace * np.eye(3) - metric[1:, 1:]) @ d1[1:])
    else:
        phases = np.exp(2j * np.asarray(N.OFFS[4:]) @ (k*h))
        metric_grid = (x[4:] * phases)[:, None, None, None] * phase[None]**2
        gradient_force = R.spatial(phase, metric_grid, SimpleNamespace(h=h), True)[0, 0, 0]
    driving = gradient_force + trace * (K(1) @ K(1)) + 2*w*w*trace
    psi = driving / (1 + K(3) @ K(3) - (3*w)**2)

    if h is None:
        grad_psi = d3 * psi
        t_scalar = (ETA @ (np.outer(d1, grad_psi) + np.outer(grad_psi, d1))
                    - IDENTITY * (d1 @ ETA @ grad_psi + psi))
    else:
        perturbation = psi * phase**3
        hi = source_from_grid(phase + perturbation, -1j*w*(phase + 3*perturbation), h, 4*k)
        lo = source_from_grid(phase - perturbation, -1j*w*(phase - 3*perturbation), h, 4*k)
        t_scalar = (hi - lo) / 2
    t_metric = inv1 @ np.outer(d1, d1) - IDENTITY * (d1 @ inv1 @ d1) / 2
    t1 = t_scalar + t_metric
    effective = t1 - g2
    lowered = ETA @ effective
    ward = d4 @ effective
    christoffel = np.einsum("mms,sn->n", gamma1, g1) - np.einsum("smn,ms->n", gamma1, g1)
    scale = np.linalg.norm(t1) + np.linalg.norm(g2)
    div_scale = np.linalg.norm(d4) * scale
    return dict(h=h, direction=list(direction), first_resid=float(first_resid),
                divergence=float(np.linalg.norm(ward) / div_scale),
                symmetry=float(np.linalg.norm(lowered - lowered.T) / scale),
                geometric_bianchi=float(np.linalg.norm(d4 @ g2 + christoffel) / div_scale),
                matter_balance=float(np.linalg.norm(d4 @ t1 + christoffel) / div_scale),
                product_defect=float(np.linalg.norm(K(2) - 2*K(1)) / np.linalg.norm(2*K(1))))


def main():
    start = time.process_time()
    directions = ((1, 0, 0), (1, 1, 1))
    continuum = [run(d) for d in directions]
    lattice = [run(d, h) for d in directions for h in (.5, .25, .125)]
    v0 = all(x["first_resid"] <= 1e-10 for x in continuum + lattice) and all(
        max(x["divergence"], x["symmetry"]) <= 1e-10 for x in continuum)
    v1 = all(max(x["divergence"], x["symmetry"]) <= 1e-10 for x in lattice)
    output = dict(V0=v0, V1=v1, continuum=continuum, lattice=lattice,
                  cpu_seconds=time.process_time()-start, python=platform.python_version(), numpy=np.__version__)
    path = Path(__file__).resolve().parents[1] / "results" / "EX-GC34V"
    path.mkdir(parents=True, exist_ok=True)
    (path / "results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
