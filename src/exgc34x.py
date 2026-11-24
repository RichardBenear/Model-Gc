"""EX-GC34X: solve on-shell gauge conditions for short-range cubic stencils."""
import os

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import json
import platform
import time
from itertools import permutations, product
from pathlib import Path

import numpy as np

import exgc34w as W

AXIS = np.eye(3, dtype=int)
SHIFTS = [tuple(v) for v in np.vstack([np.zeros((1, 3), int), AXIS, -AXIS])]
OFFSETS = np.array([(u, v) for u, v in product(SHIFTS, repeat=2)
                    if (u, v) > (tuple(-i for i in u), tuple(-i for i in v))])
PERMUTATIONS = list(permutations(range(3)))
POLARIZATIONS = list(product(range(2), repeat=2))
X, YZ, DIAGONAL = (1, 0, 0), (0, 1, 1), (1, 1, 1)
BASE_CHANNELS = [((d, d), (p, p)) for d in (X, DIAGONAL) for p in range(2)]
BASE_CHANNELS += [((X, YZ), p) for p in POLARIZATIONS]


def sampled_pairs(rng, count, used):
    pairs = []
    while len(pairs) < count:
        a, b = map(tuple, rng.integers(-3, 4, size=(2, 3)).tolist())
        pair = tuple(sorted((a, b)))
        if any(np.cross(a, b)) and pair not in used:
            used.add(pair)
            pairs.append(pair)
    return pairs


def ordered(legs):
    (dens, inv, _, _), (_, _, low, g1), (_, _, _, g2) = legs
    gamma2 = np.einsum("rs,smn->rmn", inv, low)
    return np.array([W.contract(dens, g1, g2), W.contract(W.ETA, gamma2, g2),
                     W.contract(W.ETA, g2, gamma2)])


def channel(pair, polarizations, h):
    k1, k2 = [np.array(d) * np.pi / 4 for d in pair]
    a, w1, ea = W.wave(k1, h, polarizations[0])
    b, w2, eb = W.wave(k2, h, polarizations[1])
    k3, w3 = -k1-k2, -w1-w2
    q = np.array([k1, k2, k3]) * (h or 0)
    corrections = []
    for i, j, _ in PERMUTATIONS:
        corrections.append(np.cos(OFFSETS[:, 0] @ q[i] + OFFSETS[:, 1] @ q[j]) - 1)
    metric_basis = [basis.reshape(4, 4) for basis in W.N.PM.T]
    baseline, terms, reference = [], [], []
    for m in metric_basis:
        c = W.leg(m, k3, w3, h)
        legs = (a, b, c)
        values = np.array([ordered([legs[i] for i in order]) for order in PERMUTATIONS])
        baseline.append(values.sum())
        reference.append(W.vertex(a, b, c))
        terms.append((values.T @ corrections).ravel())
    baseline, terms, reference = map(np.array, (baseline, terms, reference))
    ward, rows, reference_ward = [], [], []
    p = np.r_[-w3, W.dispersion(k3, h)]
    for i in range(4):
        m = W.gauge_metric(p, i)
        coords = np.r_[-m[0, 0]/2, m[0, 1:], [m[a+1, b+1] for a, b in W.N.PAIRS]]
        scale = np.linalg.norm(baseline) * np.linalg.norm(coords)
        ward.append((baseline @ coords) / scale)
        rows.append((coords @ terms) / scale)
        reference_ward.append((reference @ coords) / scale)
    return dict(pair=pair, polarization=polarizations, h=h, ward=np.array(ward),
                reference=np.array(reference_ward), matrix=np.array(rows), base_source=baseline,
                correction_source=terms, wave_residual=max(ea, eb))


def reported(case, coeffs):
    fitted = case["ward"] + case["matrix"] @ coeffs
    return dict(pair=case["pair"], polarization=case["polarization"], h=case["h"],
                baseline=np.abs(case["ward"]).tolist(), fitted=np.abs(fitted).tolist())


def main():
    start = time.process_time()
    rng = np.random.default_rng(343424)
    used = {tuple(sorted(pair)) for pair, _ in BASE_CHANNELS}
    seen_pairs = sampled_pairs(rng, 40, used)
    fresh_pairs = sampled_pairs(rng, 20, used)
    baseline = [channel(pair, p, .25) for pair, p in BASE_CHANNELS]
    training = baseline + [channel(pair, p, .25) for pair in seen_pairs for p in POLARIZATIONS]
    M = np.vstack([case["matrix"] for case in training])
    y = -np.concatenate([case["ward"] for case in training])
    A, b = np.vstack((M.real, M.imag)), np.r_[y.real, y.imag]
    U, s, Vh = np.linalg.svd(A, full_matrices=False)
    rows = []
    for cutoff in (1e-10, 1e-12, 1e-13):
        rank = int(np.sum(s > cutoff*s[0]))
        retained = U[:, :rank]
        coeffs = Vh[:rank].T @ ((retained.T @ b) / s[:rank])
        residual = A @ coeffs - b
        projection = retained @ (retained.T @ b) - b
        rows.append(dict(cutoff=cutoff, rank=rank, coefficient_norm=float(np.linalg.norm(coeffs)),
                         residual_l2=float(np.linalg.norm(residual)),
                         lower_bound_on_max=float(np.linalg.norm(projection) / np.sqrt(len(y))),
                         left_projection_agreement=float(np.linalg.norm(residual-projection)),
                         max_gauge=float(np.max(abs(M @ coeffs-y)))))
        if cutoff == 1e-12:
            fitted = coeffs

    holdout = []
    for h in (.5, .125):
        for pair in fresh_pairs:
            if np.max(abs((np.array(pair[0])+np.array(pair[1])) * h / 4)) >= 1:
                continue
            holdout.extend(channel(pair, p, h) for p in POLARIZATIONS)
    validation = [reported(case, fitted) for case in holdout]
    convergence = []
    for pair, p in BASE_CHANNELS[4:]:
        reference = channel(pair, p, None)["base_source"]
        for h in (.25, .125, .0625):
            case = channel(pair, p, h)
            source = case["base_source"] + case["correction_source"] @ fitted
            convergence.append(dict(**reported(case, fitted),
                                    source_error=float(np.linalg.norm(source-reference)/np.linalg.norm(reference))))
    reconstruction = max(float(np.max(abs(case["ward"]-case["reference"]))) for case in baseline)
    wave = max(case["wave_residual"] for case in training+holdout)
    x0 = reconstruction <= 1e-11 and wave <= 1e-12
    x1 = rows[1]["max_gauge"] <= 1e-10
    x2 = max(max(case["fitted"]) for case in validation) <= 1e-10
    result = dict(X0=x0, X1=x1, X2=x2, x2_scored=x1, reconstruction=reconstruction, wave=wave,
                  matrix_shape=list(A.shape), singular_values=s.tolist(), ranks=rows,
                  offsets=OFFSETS.tolist(), coefficients=fitted.tolist(),
                  train=[reported(case,fitted) for case in training], holdout=validation,
                  continuum_convergence=convergence, cpu_seconds=time.process_time()-start,
                  python=platform.python_version(), numpy=np.__version__)
    out = Path(__file__).resolve().parents[1] / "results" / "EX-GC34X"
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(result, indent=2)+"\n")
    print(x0,x1,x2, rows[1], "holdout maxima",
          [max(max(case["fitted"]) for case in validation if case["h"] == h) for h in (.5,.125)], flush=True)


if __name__ == "__main__":
    main()
