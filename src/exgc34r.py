"""EX-GC34R: first scalar response on its zeroth-order sourced metric."""
import os

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import json
import platform
import time

import numpy as np

import exgc34p as P
from exgc34q import StandingPacket

P.LBOX, P.SIG, P.NCUT = 16.0, 2.0, 8
EPS = 1e-4


def packet(n, lattice):
    out = StandingPacket(n, "generic", lattice)
    for k in out.KK:
        out.a[np.abs(k * P.LBOX / (2 * np.pi)) > 4.01] = 0
    return out


def forward(x, a):
    return np.roll(x, -1, axis=a)


def backward(x, a):
    return np.roll(x, 1, axis=a)


def average(x, a):
    return (x + forward(x, a)) / 2


def spatial(phi, metric, grid, lattice):
    """Coefficient of the metric in the spatial Euler–Lagrange force."""
    if not lattice:
        gradients = [P.sf.ifftn(1j * k * P.sf.fftn(phi)).real for k in grid.KK]
        currents = [metric[:3].sum(0) / 2 * g for g in gradients]
        for value, (a, b) in zip(metric, P.N.PAIRS):
            currents[a] -= value * gradients[b]
            if a != b:
                currents[b] -= value * gradients[a]
        return sum(P.sf.ifftn(1j * k * P.sf.fftn(j)).real for k, j in zip(grid.KK, currents))
    G = P.N.grads(phi, grid.h)
    trace = sum(metric[:3])
    currents = [0.5 * average(trace - metric[a], a) * G[a]
                - 0.25 * (metric[a] * backward(G[a], a)
                          + forward(metric[a] * G[a], a)) for a in range(3)]
    for value, (a, b) in zip(metric[3:], P.N.PAIRS[3:]):
        flux_a = value * average(G[b], a)
        flux_b = value * average(G[a], b)
        currents[a] -= (flux_a + backward(flux_a, b)) / 2
        currents[b] -= (flux_b + backward(flux_b, a)) / 2
    return sum((j - backward(j, a)) / grid.h for a, j in enumerate(currents))


def lap(phi, grid, lattice):
    if lattice:
        return P.N.lap(phi, grid.h)
    return P.sf.ifftn(-sum(k**2 for k in grid.KK) * P.sf.fftn(phi)).real


def action_check():
    rng = np.random.default_rng(343418)
    shape = (8, 8, 8)
    phi, vel, probe, probe_v = rng.normal(size=(4, *shape))
    hd, ho = rng.normal(size=(2, 3, *shape))
    zero = np.zeros(shape)
    grid = packet(8, True)
    grid.h = 0.5
    def delta(x, y):
        return (P.O.mlagrangian2(x, y, zero, 0 * hd, hd, ho, grid.h)
                - P.O.mlagrangian2(x, y, zero, 0 * hd, 0 * hd, 0 * ho, grid.h))
    trace = hd.sum(0) / 2
    force = spatial(phi, np.concatenate((hd, ho)), grid, True) - P.MASS**2 * trace * phi
    predicted = grid.h**3 * np.sum(force * probe + trace * vel * probe_v)
    step = 1e-4
    measured = (delta(phi + step * probe, vel + step * probe_v)
                - delta(phi - step * probe, vel - step * probe_v)) / (2 * step)
    return float(abs(measured - predicted) / abs(predicted))


def run(n, dt, lattice, finite=False):
    grid = packet(n, lattice)
    modes = P.mask_modes()
    k = modes * (2 * np.pi / P.LBOX)
    K = 2 * np.sin(k * grid.h / 2) / grid.h if lattice else k
    system = P.System(K)
    def force(t):
        if lattice:
            return P.lat_f(grid, t, 1, modes, k)
        return P.cont_f(grid, t, 1, modes)

    ids = tuple((modes % n).T)
    offsets = np.array(P.N.OFFS[4:]) if lattice else np.zeros((6, 3))
    phases = n**3 * np.exp(1j * (k @ (offsets * grid.h).T))
    def field(coeff):
        out = np.zeros((n, n, n, 6), complex)
        out[ids] = coeff * phases
        return P.sf.ifftn(out, axes=(0, 1, 2)).real.transpose(3, 0, 1, 2)

    epsilons = (EPS, -EPS, EPS / 2, -EPS / 2) if finite else ()
    f0 = force(0)
    metric, metric_v, _ = system.init(f0)
    phi0, v0 = grid.field(0)
    tangent = np.zeros_like(phi0)
    states = (metric, metric_v, tangent, tangent.copy(),
              *(a.copy() for _ in epsilons for a in (phi0, v0)))
    maximum_metric = 0.0
    def rhs(t, state, f):
        nonlocal maximum_metric
        metric, metric_v, tangent, tangent_v, *other = state
        h = field(metric)
        if finite:
            maximum_metric = max(maximum_metric, float(np.max(abs(h))))
        trace = sum(h[:3]) / 2
        trace_v = sum(field(metric_v)[:3]) / 2
        phi, phid = grid.field(t)
        base_lap = lap(phi, grid, lattice)
        driving = spatial(phi, h, grid, lattice) - trace * base_lap - trace_v * phid
        result = [metric_v, system.acc(metric, metric_v, f[:, P.IH]), tangent_v,
                  lap(tangent, grid, lattice) - P.MASS**2 * tangent + driving]
        for j, eps in enumerate(epsilons):
            x, v = other[2*j:2*j+2]
            accel = ((lap(x, grid, lattice) + eps * (spatial(x, h, grid, lattice) - trace_v * v))
                     / (1 + eps * trace) - P.MASS**2 * x)
            result.extend((v, accel))
        return result

    def step_add(state, slope, scale):
        return [x + scale * y for x, y in zip(state, slope)]

    rows = []
    max_res = 0.0
    current = f0
    for i in range(round(2 / dt) + 1):
        t = i * dt
        residual, _ = system.cres(states[0], states[1], current)
        max_res = max(max_res, float(np.linalg.norm(residual) / np.linalg.norm(current[:, P.IC])))
        if abs(t - 1) < 1e-10 or abs(t - 2) < 1e-10:
            phi = grid.field(t)[0]
            trace = sum(field(states[0])[:3]) / 2
            volume = grid.h**3
            scalar = 2 * volume * np.sum(phi * states[2])
            measure = volume * np.sum(trace * phi**2)
            row = dict(t=t, scalar=float(scalar), volume=float(measure), response=float(scalar + measure),
                       flat_integral=float(volume * np.sum(phi**2)))
            for j in range(0, len(epsilons), 2):
                ep = epsilons[j]
                positive = volume * np.sum((1 + ep * trace) * states[4 + 2*j]**2)
                negative = volume * np.sum((1 - ep * trace) * states[6 + 2*j]**2)
                row[f"fd_{j//2}"] = float((positive - negative) / (2 * ep))
            rows.append(row)
        if i == round(2 / dt):
            break
        middle, end = force(t + dt/2), force(t + dt)
        a = rhs(t, states, current)
        b = rhs(t + dt/2, step_add(states, a, dt/2), middle)
        c = rhs(t + dt/2, step_add(states, b, dt/2), middle)
        d = rhs(t + dt, step_add(states, c, dt), end)
        states = [x + dt/6 * (y + 2*z + 2*w + u) for x, y, z, w, u in zip(states, a, b, c, d)]
        current = end
    result = dict(n=n, dt=dt, lattice=lattice, max_constraint=max_res, rows=rows,
                  max_metric_times_epsilon=maximum_metric * EPS)
    print(result, flush=True)
    return result


def main():
    start = time.process_time()
    output = P.ROOT / "results" / "EX-GC34R"
    output.mkdir(parents=True, exist_ok=True)
    action = action_check()
    initial = [packet(n, lat).field(0)[0][::n//8, ::n//8, ::n//8]
               for n, lat in ((32, True), (64, True), (40, False))]
    initial_error = max(np.max(abs(a - initial[0])) / np.max(abs(initial[0])) for a in initial)
    print("action check", action, "initial profiles", initial_error, flush=True)
    coarse = run(32, 0.025, True, finite=True)
    fine = run(64, 0.025, True)
    continuum = run(40, 0.025, False)
    half = run(32, 0.0125, True)
    comparisons = []
    for c, f, r, t in zip(coarse["rows"], fine["rows"], continuum["rows"], half["rows"]):
        target = r["response"]
        ce, fe = abs(c["response"] - target), abs(f["response"] - target)
        comparisons.append(dict(t=r["t"], fine_rel=fe / abs(target), ratio=ce / fe,
                                scalar_fine_rel=abs(f["scalar"] / r["scalar"] - 1),
                                temporal_fraction=abs(c["response"] - t["response"]) / ce,
                                finite_rel=abs(c["fd_1"] / c["response"] - 1)))
    gates = dict(R0=bool(action <= 1e-7 and initial_error <= 1e-12
                         and max(r["max_constraint"] for r in (coarse, fine, continuum, half)) <= 1e-6),
                 R1=bool(all(r["fine_rel"] <= .05 and r["ratio"] >= 3 for r in comparisons)
                         and all(r["scalar"] != 0 for r in fine["rows"])),
                 R2=bool(all(r["temporal_fraction"] <= .1 and r["finite_rel"] <= .01 for r in comparisons)))
    results = dict(gates=gates, overall=all(gates.values()), action_check_rel=action,
                   initial_profile_rel=float(initial_error), coarse=coarse, fine=fine,
                   continuum=continuum, half=half, comparisons=comparisons,
                   platform=dict(python=platform.python_version(), numpy=np.__version__, system=platform.system()),
                   cpu_seconds=time.process_time() - start)
    (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    print(gates, comparisons, flush=True)


if __name__ == "__main__":
    main()
