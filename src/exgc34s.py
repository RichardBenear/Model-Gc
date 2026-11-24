"""EX-GC34S: apply the flat scalar source a second time and audit constraints."""
import json
import platform
import time

import numpy as np

import exgc34r as R

P = R.P
STRENGTHS = (1e-4, 5e-5)


def source(state, h):
    return P.O.msources2(*state, h)


def tangent_source(a, b, h):
    """Derivative of the quadratic source at a in direction b."""
    plus = source(tuple(x + y for x, y in zip(a, b)), h)
    minus = source(tuple(x - y for x, y in zip(a, b)), h)
    return tuple((x - y) / 2 for x, y in zip(plus, minus))


def add(a, b, scale=1):
    return tuple(x + scale * y for x, y in zip(a, b))


def norm(x):
    return float(np.linalg.norm(x))


def identities(s, sd, h):
    energy, momentum, *_ = P.N.midentities(s, sd, h)
    return energy, momentum


def audit(phi, vel, acc, tangent, tangent_v, tangent_a, driving, h):
    u, ud = (phi, vel), (vel, acc)
    z, zd = (tangent, tangent_v), (tangent_v, tangent_a)
    base = source(u, h)
    base_dot = tangent_source(u, ud, h)
    s1 = tangent_source(u, z, h)
    sd1 = add(tangent_source(ud, z, h), tangent_source(u, zd, h))
    b = identities(base, base_dot, h)
    r = identities(s1, sd1, h)
    scales = [max(norm(s), np.finfo(float).tiny) for s in base_dot[:2]]
    # Closed-form local force density when the scalar acceleration is off shell.
    gradients = P.N.grads(phi, h)
    product = (h**3 * vel * driving,
               h**3 * np.array([R.average(driving, a) * gradients[a] for a in range(3)]))
    product_error = [norm(x - y) / max(norm(y), np.finfo(float).tiny) for x, y in zip(r, product)]
    finite = []
    for eps in STRENGTHS:
        updated, updated_dot = add(u, z, eps), add(ud, zd, eps)
        imbalance = identities(source(updated, h), tangent_source(updated, updated_dot, h), h)
        measured = [(x - y) / eps for x, y in zip(imbalance, b)]
        finite.append(dict(eps=eps, normalized_defect=[norm(x) * eps / a for x, a in zip(measured, scales)],
                           relative_to_tangent=[norm(x-y) / norm(y) for x, y in zip(measured, r)],
                           signed_alignment=[float(np.vdot(x, y).real / (norm(x) * norm(y)))
                                             for x, y in zip(measured, r)]))
    return dict(flat=[norm(x) / a for x, a in zip(b, scales)],
                tangent=[norm(x) / a for x, a in zip(r, scales)],
                product_error=product_error, finite=finite)


def run(n, dt):
    grid = R.packet(n, True)
    h = grid.h
    modes = P.mask_modes()
    k, K = P.lattice_K(modes, h)
    system = P.System(K)
    ids = tuple((modes % n).T)
    offsets = np.array(P.N.OFFS)
    phases = np.exp(1j * (k @ (offsets * h).T))

    def field(coeff):
        out = np.zeros((n, n, n, 6), complex)
        out[ids] = coeff * phases[:, 4:] * n**3
        return P.sf.ifftn(out, axes=(0, 1, 2)).real.transpose(3, 0, 1, 2)

    def force(s):
        comps = [s[0], *s[1], *s[2], *s[3]]
        transformed = np.column_stack([P.sf.fftn(x)[ids] for x in comps])
        return transformed * phases.conj() * (2 * P.KAP / (n**3 * h**3))

    def flat(t):
        phi, vel = grid.field(t)
        return phi, vel, force(source((phi, vel), h))

    base = flat(0)
    metric, metric_v, _ = system.init(base[2])
    state = (metric, metric_v, np.zeros_like(base[0]), np.zeros_like(base[0]),
             np.zeros_like(metric), np.zeros_like(metric_v))

    def evaluate(state, base):
        m, mv, z, zv, dm, dmv = state
        phi, vel, f0 = base
        g = field(m)
        trace = sum(g[:3]) / 2
        trace_dot = sum(field(mv)[:3]) / 2
        base_lap = R.lap(phi, grid, True)
        drive = R.spatial(phi, g, grid, True) - trace * base_lap - trace_dot * vel
        za = R.lap(z, grid, True) - P.MASS**2 * z + drive
        f1 = force(tangent_source((phi, vel), (z, zv), h))
        slopes = (mv, system.acc(m, mv, f0[:, P.IH]), zv, za,
                  dmv, system.acc(dm, dmv, f1[:, P.IH]))
        return slopes, f1, drive, base_lap

    rows = []
    max_base = max_combined = 0.0
    for i in range(round(2 / dt) + 1):
        t = i * dt
        a, f1, drive, base_lap = evaluate(state, base)
        f0 = base[2]
        r0, _ = system.cres(state[0], state[1], f0)
        r1, _ = system.cres(state[4], state[5], f1)
        scale = norm(f0[:, P.IC])
        max_base = max(max_base, norm(r0) / scale)
        combined = norm(r0 + STRENGTHS[0] * r1) / scale
        max_combined = max(max_combined, combined)
        if abs(t - round(t)) < 1e-10:
            entry = dict(t=t, metric_correction_residual_over_base_source=norm(r1) / scale,
                         combined_constraint_residual=combined)
            if t > 0:
                entry.update(audit(base[0], base[1], base_lap - P.MASS**2 * base[0],
                                   state[2], state[3], a[3], drive, h))
            else:
                s = source(base[:2], h)
                sd = tangent_source(base[:2], (base[1], base_lap - P.MASS**2 * base[0]), h)
                entry["flat"] = [norm(x) / norm(y) for x, y in zip(identities(s, sd, h), sd[:2])]
            rows.append(entry)
        if i == round(2 / dt):
            break
        half, end = flat(t + dt / 2), flat(t + dt)
        b = evaluate(add(state, a, dt / 2), half)[0]
        c = evaluate(add(state, b, dt / 2), half)[0]
        d = evaluate(add(state, c, dt), end)[0]
        state = tuple(x + dt / 6 * (y + 2*z + 2*w + v) for x, y, z, w, v in zip(state, a, b, c, d))
        base = end
    out = dict(n=n, dt=dt, max_base_residual=max_base, max_combined_residual=max_combined, rows=rows)
    print(out, flush=True)
    return out


def main():
    start = time.process_time()
    coarse, fine, half = run(32, .025), run(64, .025), run(32, .0125)
    runs = (coarse, fine, half)
    positive = [row for case in runs for row in case["rows"] if row["t"] > 0]
    s0 = all(case["max_base_residual"] <= 1e-6 and
             all(max(row["flat"]) <= 1e-10 for row in case["rows"]) for case in runs)
    temporal = [abs(a/b - 1) for x, y in zip(coarse["rows"][1:], half["rows"][1:])
                for a, b in zip(x["tangent"], y["tangent"]) if a > 1e-4]
    s1 = ("DIRECT_RESOURCING_RULED_OUT"
          if any(max(row["tangent"]) > 1e-4 for row in fine["rows"][1:]) and max(temporal, default=0) <= .01
          else "COMPATIBLE" if all(max(row["tangent"]) <= 1e-10 for row in positive) else "UNRESOLVED")
    s2 = all(max(row["finite"][1]["relative_to_tangent"]) <= .01 for row in positive)
    result = dict(S0=s0, S1=s1, S2=s2, max_temporal_relative=max(temporal, default=0), runs=runs,
                  cpu_seconds=time.process_time()-start,
                  platform=dict(python=platform.python_version(), numpy=np.__version__, system=platform.system()))
    out = P.ROOT / "results" / "EX-GC34S"
    out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(s0, s1, s2, flush=True)


if __name__ == "__main__":
    main()
