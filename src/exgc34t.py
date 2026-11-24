"""EX-GC34T continuum Bianchi and homogeneous constraint checks."""
import json
import platform
import time
from itertools import product
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


def einstein_coefficients(coords, perturbation):
    """Mixed Einstein coefficients for g = eta + epsilon * perturbation."""
    eta = sp.diag(-1, 1, 1, 1)
    inverse_first = -eta * perturbation * eta
    indices = range(4)
    triples = list(product(indices, repeat=3))

    def connection(inverse):
        return { (r, m, n): sum(
            inverse[r, s] * (sp.diff(perturbation[s, n], coords[m])
                             + sp.diff(perturbation[s, m], coords[n])
                             - sp.diff(perturbation[m, n], coords[s])) / 2
            for s in indices) for r, m, n in triples }

    first, second = connection(eta), connection(inverse_first)

    def ricci(gamma, quadratic):
        def entry(m, n):
            value = sum(sp.diff(gamma[r, m, n], coords[r]) - sp.diff(gamma[r, m, r], coords[n])
                        for r in indices)
            if quadratic:
                value += sum(first[r, r, s] * first[s, m, n] - first[r, n, s] * first[s, m, r]
                             for r, s in product(indices, repeat=2))
            return sp.expand(value)
        return sp.Matrix(4, 4, entry)

    mixed_first = eta * ricci(first, False)
    mixed_second = eta * ricci(second, True) + inverse_first * ricci(first, False)
    return (mixed_first - sp.eye(4) * sp.trace(mixed_first) / 2,
            mixed_second - sp.eye(4) * sp.trace(mixed_second) / 2, first)


def symbolic_check():
    t, x, y, z = coords = sp.symbols("t x y z", real=True)
    h = sp.diag(0, (1 + t) * sp.cos(x), (t + t**2) * sp.sin(x), (t + t**2) * sp.sin(x))
    g1, g2, gamma = einstein_coefficients(coords, h)

    def divergence(tensor):
        return [sp.simplify(sum(sp.diff(tensor[m, n], coords[m]) for m in range(4))) for n in range(4)]

    def connection_correction(tensor):
        return [sp.simplify(sum(gamma[m, m, s] * tensor[s, n] - gamma[s, m, n] * tensor[m, s]
                                for m, s in product(range(4), repeat=2))) for n in range(4)]

    d1 = divergence(g1)
    d2 = divergence(g2)
    correction = connection_correction(g1)
    complete = [sp.simplify(a + b) for a, b in zip(d2, correction)]
    # Use the same G1, subtract a constant positive density from T0 only.
    mismatched = connection_correction(sp.diag(1, 0, 0, 0))
    no_mean = [sp.simplify(v - sp.integrate(v, (x, 0, 2 * sp.pi)) / (2 * sp.pi))
               for v in mismatched]
    return dict(pass_gate=all(v == 0 for v in d1 + complete)
                          and any(v != 0 for v in correction) and any(v != 0 for v in no_mean),
                first_divergence=list(map(str, d1)), second_divergence=list(map(str, d2)),
                connection_correction=list(map(str, correction)), completed=list(map(str, complete)),
                subtraction_control_after_mean_removal=list(map(str, no_mean)))


def homogeneous_check():
    rows = []
    for eps in (1e-4, 2.5e-5):
        rho0 = 1.0
        h0 = np.sqrt(eps * rho0 / 3)

        def rhs(_, state):
            a, H, phi, vel = state
            return (a * H, -eps * vel**2 / 2, vel, -3 * H * vel - phi)

        sol = solve_ivp(rhs, (0, 2), (1, h0, 1, 1), t_eval=np.linspace(0, 2, 101),
                        method="DOP853", rtol=1e-12, atol=1e-14)
        if not sol.success:
            raise RuntimeError(sol.message)
        a, H, phi, vel = sol.y
        rho = (vel**2 + phi**2) / 2
        rel = float(np.max(abs(3 * H**2 / (eps * rho) - 1)))
        rows.append(dict(epsilon=eps, H0=float(h0), max_friedmann_relative=rel,
                         static_initial_relative=abs(0 / (eps * rho0) - 1),
                         final_a=float(a[-1]), final_density=float(rho[-1])))
    ratio = rows[0]["H0"] / rows[1]["H0"]
    return dict(pass_gate=abs(ratio - 2) <= 1e-10 and
                          all(r["max_friedmann_relative"] <= 1e-8 and
                              r["static_initial_relative"] == 1 for r in rows),
                H0_ratio=ratio, rows=rows)


def main():
    start = time.process_time()
    t1, t2 = symbolic_check(), homogeneous_check()
    result = dict(T1=t1, T2=t2, overall=t1["pass_gate"] and t2["pass_gate"],
                  cpu_seconds=time.process_time() - start,
                  platform=dict(python=platform.python_version(), sympy=sp.__version__))
    output = Path(__file__).resolve().parents[1] / "results" / "EX-GC34T"
    output.mkdir(parents=True, exist_ok=True)
    (output / "results.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
