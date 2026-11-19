"""EX-GC13 lattice part: can the lattice knots perform the test loops through allowed moves?

Each continuous loop of EX-GC13 (src/exgc13_loopclass.py) is sampled at K times, every configuration is snapped to
the nearest of the 48 lifts on an n^3 lattice (as EX-GC12), and each transition between consecutive configurations
is replayed as single-node moves under EX-GC08's three-node lift rule (greedy order, repeated passes). A transition is
'clean' if every changed node can be moved by an allowed step. Clean loops are homotopic to the continuous loops
(straight-line interpolation of nearby lifts never passes through zero), so they carry the same class.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc08_belt as b  # noqa: E402
import exgc12_soliton as s12  # noqa: E402
import exgc13_loopclass as lc  # noqa: E402

OUT = lc.OUT
STEP = set(int(i) for i in np.nonzero([b.cls(i) in ("90", "120") for i in range(b.N)])[0])


def snap(U, t, n, scale):
    L = (n - 1) // 2
    ax = (np.arange(n) - L) * scale
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    xs = np.stack([X.ravel(), Y.ravel(), Z.ravel()], -1)
    q = U(t, xs)
    return np.argmax(q @ s12.Q.T, axis=-1).reshape(n, n, n)


def replay(old, new, A):
    """Apply the changed nodes one by one under the lift rule; return the number that could not be applied."""
    cur = old.copy()
    todo = [tuple(i) for i in np.argwhere(old != new)]
    bad_step = 0
    for idx in todo:                                           # frame step must be one allowed step (or 0)
        g, h = old[idx] // 2, new[idx] // 2
        if g != h and b.MUL[b.INV[g], h] not in STEP:
            bad_step += 1
    n = old.shape[0]
    progress = True
    while todo and progress:
        progress = False
        rest = []
        for idx in todo:
            g, h = cur[idx] // 2, new[idx] // 2
            ok = True
            if g != h:
                if b.MUL[b.INV[g], h] not in STEP:
                    ok = False
                else:
                    for d in s12.NB:
                        j = tuple(np.add(idx, d))
                        if min(j) < 0 or max(j) >= n:
                            continue
                        if not A[cur[j] // 2, g, h]:
                            ok = False
                            break
            if ok and (s12.Q[cur[idx]] @ s12.Q[new[idx]]) > 0:   # the lift must follow a short step
                cur[idx] = new[idx]
                progress = True
            else:
                rest.append(idx)
        todo = rest
    return len(todo), bad_step


def run(kind, K, n=29, scale=0.25, d=3.0, R=2.0):
    _, A = b.move_table(("90", "120"), "lift")
    U = lc.make_loop(kind, d=d, R=R)
    rng = np.random.default_rng(7)
    pts = [p / np.linalg.norm(p) for p in rng.normal(size=(2, 4))]
    t0 = time.time()
    prev = snap(U, 0.0, n, scale)
    first = prev.copy()
    stuck_total, bad_total, changed_total, worst = 0, 0, 0, 0
    degs = [[s12.degree(prev, p) for p in pts]]
    mind = s12.min_dot(prev)
    for k in range(1, K + 1):
        cur = snap(U, k / K, n, scale)
        changed_total += int((cur != prev).sum())
        st, bad = replay(prev, cur, A)
        stuck_total += st
        bad_total += bad
        worst = max(worst, st)
        mind = min(mind, s12.min_dot(cur))
        if k % max(1, K // 8) == 0:
            degs.append([s12.degree(cur, p) for p in pts])
        prev = cur
    st, bad = replay(prev, first, A)                           # closing transition (tie-breaking differences)
    stuck_total += st
    bad_total += bad
    rec = dict(kind=kind, K=K, n=n, scale=scale, d=d, R=R, changed=changed_total, stuck=stuck_total,
               bad_step=bad_total, worst_transition=worst, closes=int((prev != first).sum()), degrees=degs,
               min_neighbour_dot=mind, wall=time.time() - t0)
    json.dump(rec, open(os.path.join(OUT, f"lattice_{kind}_K{K}_h{scale:g}.json"), "w"), indent=1)
    print("RESULT", kind, "K", K, "changed", changed_total, "stuck", stuck_total, "bad_step", bad_total,
          "closes", rec["closes"], "degrees", degs[0], "...", degs[-1], "min_dot", mind, "wall", round(rec["wall"]),
          flush=True)
    return rec


if __name__ == "__main__":
    run(sys.argv[1], int(sys.argv[2]), n=int(sys.argv[3]) if len(sys.argv) > 3 else 29,
        scale=float(sys.argv[4]) if len(sys.argv) > 4 else 0.25)
