"""EX-GC37B B1a + B2: token-level event counting and a carrier-sourced firing rate.

Sites fire as independent Poisson processes at rate N(x).
Two-site rule: a token on x bound for y hops when one of x, y fires and the other fired within the last TAU.
One-site control: a token hops when x fires.
Plan: spec/EX-GC37B_plan.md (amendment A1). Our construction, not LaViolette's Model G.
"""
import json
import sys
import time
from pathlib import Path
import numpy as np
from numba import njit

TAU = 0.01
TOKENS = 2000
OUT = Path(__file__).resolve().parents[1] / "results" / "EX-GC37B"


@njit(cache=True)
def walk(N, two, tau, seed):
    """Time for one token to traverse sites 0..len(N)-1 (len(N)-1 hops)."""
    np.random.seed(seed)
    t = 0.0
    nx = N[0]
    tx = -np.random.exponential(1.0 / nx)          # stationary backward recurrence
    for k in range(len(N) - 1):
        if not two:
            t += np.random.exponential(1.0 / N[k])
            continue
        ny = N[k + 1]
        ty = t - np.random.exponential(1.0 / ny)
        while True:
            t += np.random.exponential(1.0 / (nx + ny))
            if np.random.random() < nx / (nx + ny):
                if t - ty <= tau:
                    break
                tx = t
            else:
                if t - tx <= tau:
                    ty = t
                    break
                ty = t
        nx, tx = ny, ty                              # the token now sits on y; y's own history carries over
    return t


def traverse(N, two, base):
    return np.array([walk(N, two, TAU, base + s) for s in range(TOKENS)])


def rate(a, b, two):
    """Stationary coincidence rate (two-site) or firing rate (one-site) for one link."""
    if not two:
        return a
    return a * (1 - np.exp(-b * TAU)) + b * (1 - np.exp(-a * TAU))


def predicted_time(N, two):
    return float(np.sum(1.0 / rate(N[:-1], N[1:], two)))


def b1a():
    rows = {}
    L = 400
    flat = np.ones(L + 1)
    norm = {}
    for two in (True, False):
        T = traverse(flat, two, 1000 * two)
        norm[two] = T.mean() / L                    # flat time per hop: sets the unit of speed
        rows[f"flat_{'two' if two else 'one'}"] = dict(mean_time=float(T.mean()), sem=float(T.std(ddof=1) / np.sqrt(TOKENS)),
                                                         predicted=predicted_time(flat, two))
    # speed exponent at uniform N
    levels = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    for two in (True, False):
        v, pv = [], []
        for i, n in enumerate(levels):
            T = traverse(np.full(L + 1, n), two, 5000 + 100 * i + 10000 * two)
            v.append(L * norm[two] / T.mean())
            pv.append(predicted_time(flat, two) / predicted_time(np.full(L + 1, n), two))
        p = np.polyfit(np.log(levels), np.log(v), 1)[0]
        pp = np.polyfit(np.log(levels), np.log(pv), 1)[0]
        rows[f"speed_{'two' if two else 'one'}"] = dict(levels=levels, speed=[float(x) for x in v],
                                                        predicted_speed=[float(x) for x in pv], exponent=float(p),
                                                        predicted_exponent=float(pp))
    # Shapiro delay past a point mass, N = 1 - GM/r on the line y = b
    GM, Ls = 2.0, 400
    x = np.arange(-Ls, Ls + 1, dtype=float)
    for b in (8, 16, 32):
        N = 1 - GM / np.sqrt(x ** 2 + b ** 2)
        res = {}
        for two in (True, False):
            base = 20000 + b * 7 + 50000 * two
            T = traverse(N, two, base)
            dT = T - traverse(np.ones_like(N), two, base)   # paired (common random numbers) flat path
            d = dT.mean()
            P1 = predicted_time(np.ones_like(N), two)
            pred = (predicted_time(N, two) - P1) / (P1 / (2 * Ls))   # in flat hop times
            res["two" if two else "one"] = dict(delay=float(d / norm[two]), sem=float(dT.std(ddof=1) / np.sqrt(TOKENS) / norm[two]),
                                                 predicted=float(pred))
        lead = 2 * GM * np.log(2 * Ls / b)            # leading one-site (lapse-only) delay: GM*ln(4L^2/b^2)
        res["ratio"] = res["two"]["delay"] / res["one"]["delay"]
        res["ratio_err"] = res["ratio"] * np.hypot(res["two"]["sem"] / res["two"]["delay"], res["one"]["sem"] / res["one"]["delay"])
        res["predicted_ratio"] = res["two"]["predicted"] / res["one"]["predicted"]
        res["leading_one_site"] = float(lead)
        rows[f"shapiro_b{b}"] = res
    return rows


def b2():
    """Steady carrier with a prescribed sink: D lap c = s(x); N = c / c_inf (c_inf = 1).
    Periodic lattice FFT solve; the periodic box adds a uniform compensating term, so fit N = A - GM/r + B r^2."""
    n, h, S, D = 128, 1.0, 2.0, 1.0
    k = 2 * np.sin(np.pi * np.fft.fftfreq(n)) / h
    K2 = k[:, None, None] ** 2 + k[None, :, None] ** 2 + k[None, None, :] ** 2
    zero = K2 == 0
    K2[zero] = 1
    x = (np.arange(n) - n // 2) * h
    r = np.sqrt(x[:, None, None] ** 2 + x[None, :, None] ** 2 + x[None, None, :] ** 2)
    out = {}
    for core in (1.5, 3.0):
        s = np.exp(-r ** 2 / (2 * core ** 2))
        s *= S / s.sum()
        f = np.fft.fftn(np.fft.ifftshift(s))
        f[zero] = 0
        c = np.fft.fftshift(np.real(np.fft.ifftn(-f / (D * K2))))   # lattice Laplacian symbol is -K2
        rr = x[n // 2:]
        line = c[n // 2:, n // 2, n // 2]
        sel = (rr >= 3 * core) & (rr <= 24)
        A = np.column_stack([np.ones(sel.sum()), -1 / rr[sel], rr[sel] ** 2])
        coef, *_ = np.linalg.lstsq(A, line[sel], rcond=None)
        resid = line[sel] - A @ coef
        out[f"core{core}"] = dict(GM_fit=float(coef[1]), GM_continuum=S / (4 * np.pi * D),
                                  box_term_B=float(coef[2]), box_term_predicted=S / (6 * D * (n * h) ** 3),
                                  max_rel_fit_residual=float(np.max(abs(resid)) / np.max(abs(coef[1] / rr[sel]))))
    a, b = out["core1.5"]["GM_fit"], out["core3.0"]["GM_fit"]
    out["universality_rel_diff"] = abs(a - b) / abs(a)
    return out


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.process_time()
    part = sys.argv[1] if len(sys.argv) > 1 else "all"
    res = {}
    if part in ("all", "b1a"):
        res["b1a"] = b1a()
    if part in ("all", "b2"):
        res["b2"] = b2()
    res["tau"] = TAU
    res["tokens"] = TOKENS
    res["cpu_seconds"] = time.process_time() - t0
    (OUT / f"tokens_{part}.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=1))
