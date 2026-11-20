"""EX-GC25 (exploratory): short look at a medium with no rest frame, in 1+1 dimensions.

Relativistic phi^4 kinks, V = lam/4 (phi^2 - v^2)^2, stand in for particles (they are not Gc particles).
(1) Fall in a uniform scalar potential Phi = g x, coupled two ways:
    'energy': Nordstrom conformal factor; in 1+1 D only the potential term couples, L = (phi_t^2 - phi_x^2)/2 - e^{2 Phi} V.
              Kink rest energy M(x) = M0 e^{Phi(x)}, so the predicted acceleration is -g for every kink.
    'charge': L = ... - V - Phi alpha (1 - phi^2/v^2), a coupling not tied to energy; charge Q = 2 alpha w,
              alpha set so kink A has Q = M. Predicted acceleration -g Q/M.
(2) Coasting: a boosted kink, Phi = 0, on lattices of spacing h. The continuum has no drag; the lattice may.
Leapfrog, 3-point Laplacian, fixed ends phi = -+v. Our construction.
"""
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC25")
KINKS = dict(A=(1.0, 1.0), B=(4.0, 1.0), C=(1.0, 2.0), D=(0.25, 1.5))


def props(lam, v):
    return dict(M=2 * np.sqrt(2) / 3 * np.sqrt(lam) * v ** 3, w=np.sqrt(2) / (v * np.sqrt(lam)))


def centre(x, phi, v, w, guess):
    sel = np.abs(x - guess) < 6 * w
    d = 1 - (phi[sel] / v) ** 2
    return float((x[sel] * d).sum() / d.sum())


def simulate(lam, v, h, Lbox, T, x0=0.0, u=0.0, g=0.0, mode="energy", alpha=0.0, every=1.0):
    x = np.arange(-Lbox, Lbox + h / 2, h)
    w = props(lam, v)["w"]
    dt = 0.1 * h
    gam = 1 / np.sqrt(1 - u * u)
    prof = lambda t: v * np.tanh(gam * (x - x0 - u * t) / w)
    phi, prev = prof(0.0), prof(-dt)
    Phi = g * x
    fac = np.exp(2 * Phi) if mode == "energy" else np.ones_like(x)
    extra = 2 * Phi * alpha / v ** 2 if mode == "charge" else 0.0
    n = int(round(T / dt))
    ke = int(round(every / dt))
    ts, xs = [], []
    guess = x0
    for i in range(n + 1):
        if i % ke == 0:
            guess = centre(x, phi, v, w, guess)
            ts.append(i * dt)
            xs.append(guess)
        if i == n:
            break
        acc = np.zeros_like(phi)
        acc[1:-1] = (phi[2:] - 2 * phi[1:-1] + phi[:-2]) / h ** 2
        acc -= fac * lam * phi * (phi * phi - v * v)
        acc += extra * phi
        acc[0] = acc[-1] = 0.0
        new = 2 * phi - prev + dt * dt * acc
        prev, phi = phi, new
    return np.array(ts), np.array(xs)


def fall():
    g, T, h = 2e-4, 150.0, 0.05
    A = props(*KINKS["A"])
    alpha = A["M"] / (2 * A["w"])
    rows = []
    for mode in ("energy", "charge"):
        for name, (lam, v) in KINKS.items():
            p = props(lam, v)
            t, xc = simulate(lam, v, h, 200.0, T, g=g, mode=mode, alpha=alpha)
            sel = t >= 5
            a = 2 * np.polyfit(t[sel], xc[sel], 2)[0]
            pred = -g if mode == "energy" else -g * (2 * alpha * p["w"]) / p["M"]
            rows.append(dict(mode=mode, kink=name, lam=lam, v=v, M=p["M"], w=p["w"], a=float(a), a_pred=float(pred),
                             ratio=float(a / pred), a_over_g=float(-a / g)))
            print(f"fall {mode:6s} {name} M={p['M']:.3f} w={p['w']:.3f}  a/g={-a / g:+.4f}  pred {-pred / g:+.4f}"
                  f"  meas/pred {a / pred:.4f}", flush=True)
    return dict(g=g, T=T, h=h, alpha=alpha, rows=rows)


def coast():
    rows = []
    for name in ("A", "B"):
        lam, v = KINKS[name]
        p = props(lam, v)
        for h in (0.4, 0.2, 0.1):
            t, xc = simulate(lam, v, h, 150.0, 200.0, x0=-60.0, u=0.5)
            e = t <= 40
            l = t >= 160
            sel_e = e & (t >= 10)
            u0 = np.polyfit(t[sel_e], xc[sel_e], 1)[0]
            u1 = np.polyfit(t[l], xc[l], 1)[0]
            rows.append(dict(kink=name, h=h, w_over_h=p["w"] / h, u_early=float(u0), u_late=float(u1),
                             loss=float(1 - u1 / u0)))
            print(f"coast {name} h={h:<4g} w/h={p['w'] / h:5.2f}  u early {u0:.6f} late {u1:.6f}"
                  f"  loss {1 - u1 / u0:+.2e}", flush=True)
    return dict(u=0.5, T=200.0, rows=rows)


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    t0 = time.process_time()
    rec = dict(fall=fall(), coast=coast())
    rec.update(cpu_seconds=time.process_time() - t0, platform=sys.version.split()[0] + " numpy " + np.__version__)
    json.dump(rec, open(os.path.join(OUT, "EX-GC25_runs.json"), "w"), indent=1)
    print("cpu s", round(rec["cpu_seconds"], 1))
