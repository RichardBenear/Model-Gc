"""EX-GC37D: the 37B wave rule with stochastically realised coincidences.
Each face transfer at each site event is accepted with probability p = p0 * w_f (w_f = N_xN_y handshake, or N_f/N_x lapse-only),
and an accepted transfer carries weight 1/p0, so the mean equals the deterministic 37B rule.
Plan: spec/EX-GC37D_plan.md. Our construction, not LaViolette's Model G."""
import json, sys, time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.ndimage import map_coordinates
import scipy.fft as fft
import exgc35k as K


@njit
def events(u, old, period, lapse, nx, ny, step, h, ticks, rule, p0, seed):
    np.random.seed(seed)
    n = nx * ny
    history = np.zeros((n, 5), np.complex128)
    history[:, 0] = u
    history[:, 4] = old
    size = np.max(period) + 1
    head = np.full(size, -1, np.int64)
    links = np.empty(n, np.int64)
    done = np.zeros(n, np.int64)
    for i in range(n):
        at = period[i] % size
        links[i] = head[at]
        head[at] = i
    accepted = 0
    tried = 0
    for when in range(1, ticks + 1):
        at = when % size
        i = head[at]
        head[at] = -1
        while i >= 0:
            next_i = links[i]
            pi = period[i]
            earlier = when - pi
            v = done[i]
            a, b = i // ny, i % ny
            z = history[i, v % 5]
            new = 0j
            if 0 < a < nx - 1 and 0 < b < ny - 1:
                flux = 0j
                for j in (i - ny, i + ny, i - 1, i + 1):
                    nb = K.read(history, period, done, j, earlier, when, 2)
                    if rule == 1:
                        w = lapse[i] * lapse[j]
                    else:
                        w = 2 * lapse[i] * lapse[j] / (lapse[i] + lapse[j]) / lapse[i]
                    tried += 1
                    if np.random.random() < p0 * w:
                        flux += (nb - z) / p0
                        accepted += 1
                new = 2 * z - history[i, (v - 1) % 5] + step[i] * step[i] * flux / (h * h)
            history[i, (v + 1) % 5] = new
            done[i] += 1
            if when + pi <= ticks:
                nxt = (when + pi) % size
                links[i] = head[nxt]
                head[nxt] = i
            i = next_i
    output = np.empty(n, np.complex128)
    for i in range(n):
        output[i] = K.read(history, period, done, i, ticks, ticks + 1, 2)
    return output, accepted, tried


def run(pot, hg, mode, p0, seed, h=.16, step=.04):
    grid, kx, ky, u, v, a = K.phase_data(h)
    nx, ny = u.shape
    units = step / K.CLOCK_TICKS
    ticks = round(K.STOP / units)
    p = map_coordinates(pot, [z / hg + pot.shape[0] // 2 for z in grid], order=3, mode="wrap")
    target = 1 + p if mode != "flat" else np.ones_like(p)
    period = np.rint(K.CLOCK_TICKS / target).astype(np.int64)
    lapse = K.CLOCK_TICKS / period
    dt = units * period
    speed = lapse ** 2 if mode == "handshake" else lapse
    old = u - dt * (v * speed) + .5 * dt ** 2 * (a * speed ** 2)
    rule = 0 if mode == "lapse" else 1
    out, acc, tried = events(u.ravel(), old.ravel(), period.ravel(), lapse.ravel(), nx, ny, (dt * lapse).ravel(), h, ticks, rule, p0, seed)
    out = out.reshape(u.shape)
    phase, group = K.directions(h, kx, ky, out)
    P = abs(fft.fft2(out)) ** 2
    coherent = float(P[np.hypot(kx - K.CARRIER, ky) < 1].sum() / P.sum())
    P0 = abs(fft.fft2(u)) ** 2
    coherent0 = float(P0[np.hypot(kx - K.CARRIER, ky) < 1].sum() / P0.sum())
    return dict(mode=mode, p0=p0, seed=seed, group_angle=group, phase_angle=phase,
                norm_ratio=float(np.linalg.norm(out) / np.linalg.norm(u)), coherent_fraction=coherent / coherent0,
                acceptance=acc / tried, finite=bool(np.isfinite(out).all()), out_hash=hash(out.tobytes()))


if __name__ == "__main__":
    t0 = time.process_time()
    p0 = float(sys.argv[1]); mode = sys.argv[2]; seeds = [int(x) for x in sys.argv[3].split(",")]
    pot, hg, _ = K.metric()
    rows = [run(pot, hg, mode, p0, s) for s in seeds]
    out = Path(__file__).resolve().parents[1] / "results" / "EX-GC37D"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{mode}_p{p0}.json").write_text(json.dumps(dict(rows=rows, cpu_seconds=time.process_time() - t0), indent=2) + "\n")
    for r in rows: print(r, flush=True)
