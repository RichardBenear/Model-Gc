"""EX-GC34U: isolated polar-areal Einstein–Klein–Gordon benchmark."""
import json
import platform
import time
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_simpson, cumulative_trapezoid

COUPLINGS = (0.0, 0.01, -0.01, 0.005, -0.005)
TIMES = (1, 2, 4)


class RadialGrid:
    def __init__(self, outer, h):
        self.r = np.linspace(0, outer, round(outer / h) + 1)
        self.h = h

    def derivative(self, x, parity):
        result = np.gradient(x, self.h, edge_order=2)
        result[0] = 0 if parity == "even" else x[1] / self.h
        return result

    def integral(self, x, independent=False):
        fn = cumulative_trapezoid if independent else cumulative_simpson
        return fn(x, dx=self.h, initial=0)

    def over_radius(self, x):
        return np.r_[0.0, x[1:] / self.r[1:]]

    def mass(self, phi, pi, coupling):
        r = self.r
        grad = self.derivative(phi, "even")
        squares = pi**2 + grad**2
        weight = np.exp(self.integral(coupling * r * squares / 2))
        mass = self.integral(weight * r**2 * (squares + phi**2) / 4) / weight
        return mass, grad, weight

    def geometry(self, phi, pi, coupling):
        mass, grad, weight = self.mass(phi, pi, coupling)
        a = (1 - 2 * coupling * self.over_radius(mass))**(-0.5)
        pressure = ((pi**2 + grad**2) / a**2 - phi**2) / 2
        integrand = a**2 * coupling * (self.over_radius(self.over_radius(mass)) + self.r * pressure / 2)
        log_lapse = self.integral(integrand)
        lapse = np.exp(log_lapse - log_lapse[-1]) / a[-1]
        return mass, grad, weight, a, lapse

    def rhs(self, state, geometry):
        phi, pi = state
        _, grad, _, a, lapse = geometry
        # Integrate the radial divergence over shells, including the half-cell at the center.
        transport = lapse / a
        faces = (self.r[:-1] + self.r[1:]) / 2
        flux = faces**2 * (transport[:-1] + transport[1:]) / 2 * np.diff(phi) / self.h
        lower = np.r_[0.0, faces]
        upper = np.r_[faces, self.r[-1] + self.h / 2]
        divergence = (np.r_[flux, 0.0] - np.r_[0.0, flux]) * 3 / (upper**3 - lower**3)
        return np.array([lapse / a * pi, divergence - lapse * a * phi])

    def flux_check(self, state, coupling, geometry):
        phi, pi = state
        mass, grad, weight, a, lapse = geometry
        phid, pid = self.rhs(state, geometry)
        cross = pi * pid + grad * self.derivative(phid, "even")
        source = self.r**2 * (cross + phi * phid) / 2 - coupling * self.r * cross * mass
        differentiated_mass = self.integral(weight * source, independent=True) / weight
        einstein_flux = self.r**2 * lapse * pi * grad / (2 * a**3)
        return differentiated_mass, einstein_flux

    def metrics(self, state, coupling):
        geometry = self.geometry(*state, coupling)
        mass = geometry[0]
        return dict(mass=mass, state=state.copy(),
                    flux=self.flux_check(state, coupling, geometry))

    def run(self, coupling):
        profile = 0.4 * np.exp(-self.r**2 / 8)
        state = np.array([profile, profile])
        dt = self.h / 4
        m0 = self.mass(*state, coupling)[0][-1]
        drift, compactness, edge = 0.0, 0.0, 0.0
        saved = {}
        for i in range(round(4 / dt) + 1):
            t = i * dt
            geometry = self.geometry(*state, coupling)
            mass = geometry[0]
            drift = max(drift, abs(mass[-1] / m0 - 1))
            compactness = max(compactness, float(np.max(2 * coupling * self.over_radius(mass))))
            edge = max(edge, float(np.max(abs(state[:, -1]))))
            if t in TIMES:
                saved[t] = self.metrics(state, coupling)
            if t == 4:
                break
            first = self.rhs(state, geometry)
            mid = state + dt/2 * first
            second = self.rhs(mid, self.geometry(*mid, coupling))
            mid = state + dt/2 * second
            third = self.rhs(mid, self.geometry(*mid, coupling))
            end = state + dt * third
            fourth = self.rhs(end, self.geometry(*end, coupling))
            state += dt/6 * (first + 2*second + 2*third + fourth)
        checks = dict(initial_total_mass=float(coupling*m0), normalized_initial_mass=float(m0),
                      mass_drift=float(drift), compactness=compactness, boundary_amplitude=edge)
        return saved, checks


def relative(observed, reference):
    return float(np.linalg.norm(observed - reference) / np.linalg.norm(reference))


def coupling_derivative(positive, negative, epsilon):
    return (positive - negative) / (2 * epsilon)


def analyze(grid, snapshots):
    rows, signals = [], {}
    inner = grid.r <= 8
    flat = snapshots[0.0]
    for t in TIMES:
        frozen = flat[t]
        low, high = snapshots[-.005][t], snapshots[.005][t]
        measured_flux = coupling_derivative(high["flux"][0], low["flux"][0], .005)
        expected_flux = coupling_derivative(high["flux"][1], low["flux"][1], .005)
        raw = max(relative(*snapshots[c][t]["flux"]) for c in COUPLINGS)
        signal = {}
        for eps in (.01, .005):
            full = coupling_derivative(snapshots[eps][t]["mass"], snapshots[-eps][t]["mass"], eps)
            fixed = coupling_derivative(grid.mass(*frozen["state"], eps)[0],
                                        grid.mass(*frozen["state"], -eps)[0], eps)
            signal[eps] = (full - fixed)[inner]
        ds = coupling_derivative(high["state"], low["state"], .005)
        phi, pi = frozen["state"]
        dphi, dpi = ds
        grad = grid.derivative(phi, "even")
        dgrad = grid.derivative(dphi, "even")
        r = grid.r
        dmass_source = r**2 * (pi*dpi + grad*dgrad + phi*dphi) / 2
        dmass_source -= r * frozen["mass"] * (pi**2 + grad**2) / 2
        predicted = grid.integral(dmass_source, independent=True)
        measured = coupling_derivative(high["mass"], low["mass"], .005)
        rows.append(dict(t=t, raw_flux_relative=raw, coupling_flux_relative=relative(measured_flux, expected_flux),
                         signal_over_mass=float(np.sqrt(grid.h) * np.linalg.norm(signal[.005]) /
                                                (np.sqrt(8) * frozen["mass"][-1])),
                         strength_relative=relative(signal[.01], signal[.005]),
                         r2_feedback=float(signal[.005][round(2/grid.h)]),
                         predicted_mass_derivative_relative=relative(predicted[inner], measured[inner])))
        signals[t] = signal[.005]
    return rows, signals


def main():
    start = time.process_time()
    results, all_snapshots, all_signals = [], [], []
    for h in (.125, .0625, .03125):
        grid = RadialGrid(24, h)
        runs = [grid.run(c) for c in COUPLINGS]
        snapshots = dict(zip(COUPLINGS, [s for s, _ in runs]))
        rows, signals = analyze(grid, snapshots)
        all_snapshots.append(snapshots)
        all_signals.append(signals)
        results.append(dict(h=h, rows=rows, checks=dict(zip(COUPLINGS, [k for _, k in runs]))))
    wide, _ = RadialGrid(48, .0625).run(.01)
    normal = all_snapshots[1][.01][4]
    count = round(8/.0625)+1
    boundary = [relative(wide[4]["mass"][:count], normal["mass"][:count])]
    boundary.extend(relative(a[:count], b[:count]) for a, b in zip(wide[4]["state"], normal["state"]))
    coarse, _, fine = results
    fine_positive = [fine["checks"][c] for c in (.005, .01)]
    u0 = (all(c["initial_total_mass"] > 0 and c["mass_drift"] <= .005 and
              c["boundary_amplitude"] < 1e-8 and c["compactness"] < .05 for c in fine_positive)
          and max(boundary) <= 1e-5)
    u1 = all(b[k] <= .03 and a[k] >= 2*b[k] for a, b in zip(coarse["rows"], fine["rows"])
             for k in ("raw_flux_relative", "coupling_flux_relative"))
    grid_signal = relative(all_signals[1][4], all_signals[2][4][::2])
    end = fine["rows"][-1]
    u2 = end["signal_over_mass"] >= 1e-4 and end["strength_relative"] <= .05 and grid_signal <= .05
    gates = dict(U0=bool(u0), U1=bool(u1), U2=bool(u2))
    output = dict(gates=gates, overall=all(gates.values()), grids=results, boundary_relative=boundary,
                  fine_signal_grid_relative=grid_signal, cpu_seconds=time.process_time()-start,
                  platform=dict(python=platform.python_version(), numpy=np.__version__, system=platform.system()))
    path = Path(__file__).resolve().parents[1] / "results" / "EX-GC34U"
    path.mkdir(parents=True, exist_ok=True)
    (path / "results.json").write_text(json.dumps(output, indent=2) + "\n")
    print(gates, fine["rows"], "boundary", boundary, "grid signal", grid_signal)


if __name__ == "__main__":
    main()
