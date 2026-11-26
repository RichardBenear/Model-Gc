"""EX-GC37B B1b: wave version of the event-counting rule, on the EX-GC35K packet, metric source and readouts.

Each site fires at its clock rate N (proper step fixed). Per event, the face weight is
  handshake: N_x N_y                (coordinate speed^2 = N^4, speed ~ N^2; no separate spatial metric)
  lapse:     harmonic N_f / N_x     (35I/35K lapse-only control, speed ~ N)
Plan: spec/EX-GC37B_plan.md (amendment A1). Our construction, not LaViolette's Model G.
"""
import argparse, json, time
from pathlib import Path
import numpy as np
from numba import njit
from scipy.ndimage import map_coordinates
import exgc35k as K


@njit
def events(u, old, period, lapse, nx, ny, step, h, ticks, rule):
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
    count = 0
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
                        face = lapse[i] * lapse[j]
                    else:
                        face = 2 * lapse[i] * lapse[j] / (lapse[i] + lapse[j]) / lapse[i]
                    flux += face * (nb - z)
                new = 2 * z - history[i, (v - 1) % 5] + step[i] * step[i] * flux / (h * h)
            history[i, (v + 1) % 5] = new
            done[i] += 1
            count += 1
            if when + pi <= ticks:
                nxt = (when + pi) % size
                links[i] = head[nxt]
                head[nxt] = i
            i = next_i
    output = np.empty(n, np.complex128)
    for i in range(n):
        output[i] = K.read(history, period, done, i, ticks, ticks + 1, 2)
    return output, count


def run(pot, hg, h, step, mode):
    grid, kx, ky, u, v, a = K.phase_data(h)
    nx, ny = u.shape
    units = step / K.CLOCK_TICKS
    ticks = round(K.STOP / units)
    p = map_coordinates(pot, [z / hg + pot.shape[0] // 2 for z in grid], order=3, mode="wrap")
    target = 1 + p if mode != "flat" else np.ones_like(p)
    period = np.rint(K.CLOCK_TICKS / target).astype(np.int64)
    lapse = K.CLOCK_TICKS / period
    dt = units * period
    speed = lapse ** 2 if mode == "handshake" else lapse            # coordinate speed (dt below is a coordinate step)
    old = u - dt * (v * speed) + .5 * dt ** 2 * (a * speed ** 2)
    rule = 1 if mode in ("handshake", "flat") else 0
    out, count = events(u.ravel(), old.ravel(), period.ravel(), lapse.ravel(), nx, ny, (dt * lapse).ravel(), h, ticks, rule)
    out = out.reshape(u.shape)
    phase, group = K.directions(h, kx, ky, out)
    start, end = K.center(grid, u), K.center(grid, out)
    cspeed = target ** 2 if mode == "handshake" else target
    pred = K.ray_prediction(grid, cspeed, np.ones_like(p), [start[0], end[0]])
    return dict(mode=mode, h=h, proper_step=step, events=count, group_angle=group, phase_angle=phase,
                eikonal_angle=pred, initial_centroid=start, final_centroid=end,
                norm_ratio=float(np.linalg.norm(out) / np.linalg.norm(u)), finite=bool(np.isfinite(out).all()))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("tag"); ap.add_argument("--h", type=float, default=.16); ap.add_argument("--step", type=float, default=.04)
    ap.add_argument("--mode", choices=("flat", "lapse", "handshake"), default="handshake")
    a = ap.parse_args()
    t0 = time.process_time()
    pot, hg, audit = K.metric()
    row = run(pot, hg, a.h, a.step, a.mode)
    row.update(GM=K.GM, impact=K.IMPACT, ideal_full_deg=float(np.degrees(4 * K.GM / K.IMPACT)), cpu_seconds=time.process_time() - t0)
    out = Path(__file__).resolve().parents[1] / "results" / "EX-GC37B"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"wave_{a.tag}.json").write_text(json.dumps(row, indent=2) + "\n")
    print(row, flush=True)
