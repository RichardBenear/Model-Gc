"""Registered EX-GC34Q standing-packet kinetic-energy contrast."""
import os

for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import json
import platform
import time

import numpy as np
from scipy.integrate import solve_ivp

import exgc34p as P


class StandingPacket(P.Packet):
    def field(self, t, s=1):
        wt = self.w * t
        phi = np.cos(wt) + s * np.sin(wt) / self.w
        phid = -self.w * np.sin(wt) + s * np.cos(wt)
        return tuple(P.sf.ifftn(self.a * factor * self.N**3).real for factor in (phi, phid))

    def cgrads(self, t, s=1):
        phi_hat = self.a * (np.cos(self.w * t) + s * np.sin(self.w * t) / self.w) * self.N**3
        return [P.sf.ifftn(1j * k * phi_hat).real for k in self.KK]


def source_components(packet, t, c, lattice):
    if lattice:
        n, b, d, o = P.O.msources2(*packet.field(t, c), packet.h)
        return [n, *b, *d, *o]
    rho, p, stress = P.cont_T(packet, t, c)
    return [-rho, *p, *(0.5 * stress[i][i] for i in range(3)),
            *(stress[i][j] for i, j in P.N.PAIRS[3:])]


def audit(packets):
    rows = []
    cuts = (14, 18, 22, 25)
    for (direction, grid), packet in packets.items():
        lattice = grid != "continuum"
        cell = packet.h**3
        omitted = dict.fromkeys(cuts, 0.0)
        energies = [[], []]
        phi0 = packet.field(0, 0)[0]
        increment = 0.5 * cell * np.sum(phi0**2)
        profile_error = np.max(np.abs(phi0 - packet.field(0, 1)[0]))
        momenta = []
        axis = np.abs(np.fft.fftfreq(packet.N) * packet.N).round()
        mode_radius = np.maximum(np.maximum(axis[:, None, None], axis[None, :, None]), axis[None, None, :])
        for t in range(5):
            for c in (0, 1):
                phi, phid = packet.field(t, c)
                grads = P.N.grads(phi, packet.h) if lattice else packet.cgrads(t, c)
                energies[c].append(0.5 * cell * (np.sum(phid**2 + P.MASS**2 * phi**2)
                                                  + sum(np.sum(g**2) for g in grads)))
                if t == 0:
                    if lattice:
                        momentum = [-cell * np.sum(0.5 * (phid + P.N.fw(phid, a)) * grads[a]) for a in range(3)]
                    else:
                        momentum = [-cell * np.sum(phid * g) for g in grads]
                    momenta.append(float(np.linalg.norm(momentum) / increment))
            low = source_components(packet, t, 0, lattice)
            high = source_components(packet, t, 1, lattice)
            for a, b in zip(low, high):
                spectrum = P.sf.fftn(b - a)
                spectrum[0, 0, 0] = 0
                power = np.abs(spectrum)**2
                total = np.sum(power)
                if total > 0:
                    for cut in cuts:
                        if cut < packet.N // 2:
                            omitted[cut] = max(omitted[cut], float(np.sum(power[mode_radius > cut]) / total))
                        else:
                            omitted[cut] = None
        rows.append(dict(direction=direction, grid=grid, profile_error=profile_error, initial_momentum_rel=max(momenta),
                         initial_energies=[e[0] for e in energies], increment=energies[1][0] - energies[0][0],
                         predicted_increment=increment, drift=max(np.max(np.abs(np.array(e) / e[0] - 1)) for e in energies),
                         omitted=omitted))
    cut = next((c for c in cuts if all(r["omitted"][c] is not None and r["omitted"][c] <= 1e-8 for r in rows)), None)
    for r in rows:
        ref = next(x for x in rows if x["direction"] == r["direction"] and x["grid"] == "continuum")
        r["increment_error_vs_continuum"] = abs(r["increment"] / ref["increment"] - 1)
    passed = cut is not None and all(r["profile_error"] == 0 and r["initial_momentum_rel"] <= 1e-12
                                    and abs(r["increment"] / r["predicted_increment"] - 1) <= 1e-12
                                    and r["increment_error_vs_continuum"] <= 1e-3 and r["drift"] <= 1e-12 for r in rows)
    return dict(pass_gate=passed, cut=cut, rows=rows)


def solve(packet, modes, dt, lattice):
    k = modes * (2 * np.pi / P.LBOX)
    if lattice:
        _, k_eff = P.lattice_K(modes, packet.h)
        def force(t):
            return P.lat_f(packet, t, 1, modes, k) - P.lat_f(packet, t, 0, modes, k)
    else:
        k_eff = k
        def force(t):
            return P.cont_f(packet, t, 1, modes) - P.cont_f(packet, t, 0, modes)
    return P.evolve(force, P.System(k_eff), 4, dt, (0, 2, 4))


def recompute(packet, modes, tidal):
    index = np.argsort(np.sum(P.WT * abs(tidal)**2, axis=1))[-3:][::-1]
    sites = np.arange(packet.N) * packet.h
    offsets = np.array(P.N.OFFS) * packet.h
    rows = []
    for i in index:
        k, k_eff = [x[0] for x in P.lattice_K(modes[i:i + 1], packet.h)]
        phase = np.exp(-1j * (k[0] * sites[:, None, None] + k[1] * sites[None, :, None]
                              + k[2] * sites[None, None, :])) / packet.N**3
        factor = 2 * P.KAP / packet.h**3 * np.exp(-1j * offsets @ k)

        def force(t):
            hi = source_components(packet, t, 1, True)
            lo = source_components(packet, t, 0, True)
            return factor * np.array([np.sum((a - b) * phase) for a, b in zip(hi, lo)])

        def kernel(w):
            vector = np.r_[-w, k_eff]
            cols = [P.N.G_up(col.reshape(4, 4), vector).ravel() for col in P.N.PM.T]
            return P.N.PM.T @ np.column_stack(cols)

        C = kernel(0)
        plus, minus = kernel(1), kernel(-1)
        B = (plus - minus) / 2
        A = (plus + minus) / 2 - C
        ih, ic = P.IH, P.IC
        Ch, Bh, Ah = [m[np.ix_(ih, ih)] for m in (C, B, A)]
        Cc, Bc, Ac = [m[np.ix_(ic, ih)] for m in (C, B, A)]
        ai_c = np.linalg.solve(Ah, Ch)
        ai_b = np.linalg.solve(Ah, 1j * Bh)
        inv = np.linalg.inv(Ah)

        f0 = force(0)
        constraints = np.hstack((Cc - Ac @ ai_c, 1j * Bc - Ac @ ai_b))
        rhs = f0[ic] - Ac @ inv @ f0[ih]
        initial = np.linalg.lstsq(constraints, rhs, rcond=None)[0]

        def ode(t, z):
            return np.r_[z[6:], ai_c @ z[:6] + ai_b @ z[6:] - inv @ force(t)[ih]]

        integrated = solve_ivp(ode, (0, 4), initial, method="DOP853", rtol=1e-10, atol=1e-13)
        if not integrated.success:
            raise RuntimeError(integrated.message)
        independent = -0.5 * ode(4, integrated.y[:, -1])[6:]
        err = P.tidal_err(independent, tidal[i])
        rows.append(dict(mode=modes[i].tolist(), error=err, evaluations=integrated.nfev))
    return dict(pass_gate=all(r["error"] <= 1e-3 for r in rows), rows=rows)


def main():
    start = time.process_time()
    output = P.ROOT / "results" / "EX-GC34Q"
    output.mkdir(parents=True, exist_ok=True)
    directions = ("axis", "generic")
    packets = {(direction, grid): StandingPacket(n, direction, grid != "continuum")
               for direction in directions for grid, n in (("0.5", 64), ("0.25", 128), ("continuum", 96))}
    result = dict(Q0=audit(packets))
    print("Q0", result["Q0"]["pass_gate"], "cut", result["Q0"]["cut"], flush=True)
    if result["Q0"]["pass_gate"]:
        P.NCUT = result["Q0"]["cut"]
        modes = P.mask_modes()
        result["modes"] = len(modes)
        checks, errors = [], []
        for direction in directions:
            continuum, _, _, _ = solve(packets[direction, "continuum"], modes, 0.05, False)
            runs = {}
            for grid in ("0.5", "0.25"):
                tidal, max_res, init_res, rank = solve(packets[direction, grid], modes, 0.05, True)
                checked = dict(direction=direction, grid=grid, max_res=max_res, init_res=init_res, rank=rank)
                checked["pass_gate"] = init_res <= 1e-10 and max_res <= 1e-6
                runs[grid] = tidal
                if direction == "generic" and grid == "0.5":
                    half, half_res, half_init, _ = solve(packets[direction, grid], modes, 0.025, True)
                    time_fraction = P.tidal_err(half[4], tidal[4]) / P.tidal_err(continuum[4], tidal[4])
                    checked.update(half_res=half_res, half_init=half_init, temporal_to_spatial=time_fraction)
                    checked["pass_gate"] &= (half_init <= 1e-10 and half_res <= 1e-6 and time_fraction <= 0.1
                                              and (max_res >= 8 * half_res or max(max_res, half_res) <= 1e-10))
                    result["Q3"] = recompute(packets[direction, grid], modes, tidal[4])
                checks.append(checked)
                print("Q1", direction, grid, checked, flush=True)
            for t in (0, 2, 4):
                coarse, fine = (P.tidal_err(runs[g][t], continuum[t]) for g in ("0.5", "0.25"))
                errors.append(dict(direction=direction, t=t, coarse=coarse, fine=fine, ratio=coarse / fine,
                                   continuum_tidal_norm=float(np.sqrt(np.sum(P.WT * abs(continuum[t])**2)))))
            print("Q2", direction, errors[-3:], flush=True)
        result["Q1"] = dict(pass_gate=all(row["pass_gate"] for row in checks), rows=checks)
        result["Q2"] = dict(pass_gate=all(row["fine"] <= 0.02 and row["ratio"] >= 3 for row in errors), rows=errors)
    result["overall"] = all(result.get(g, {}).get("pass_gate", False) for g in ("Q0", "Q1", "Q2", "Q3"))
    result["cpu_seconds"] = time.process_time() - start
    result["platform"] = dict(python=platform.python_version(), numpy=np.__version__, system=platform.system())
    (output / "results.json").write_text(json.dumps(result, default=P.N.DUMP["default"], indent=2) + "\n")
    print("overall", result["overall"], "CPU seconds", result["cpu_seconds"], flush=True)


if __name__ == "__main__":
    main()
