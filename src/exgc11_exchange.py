"""EX-GC11 (exploratory): exchange of two tethered cores in the frame-field tether of EX-GC10 (frames only).

Two rigid cores (3x3x3 blocks, frames fixed at the identity) sit in the 25^3 lattice of EX-GC10 with the outer face
fixed. They are exchanged by pure translation (orientation kept) along a half circle, one lattice step at a time. A
core may step into a node only if setting that node to the identity is an allowed move under the three-node lift
rule (no cutting); otherwise the tether relaxes first. Between steps the tether relaxes by Metropolis at T_f.
Diagnostics: mismatched edges and the Z2 index of each core (product of short lifts along +z to the boundary).
Our construction.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exgc08_belt as b  # noqa: E402
import exgc10_tether as t10  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC11")
NB = t10.NB


def lift_z(f, c):
    """Product of short lifts from node c along +z to the outer face (sign = Z2 index)."""
    x, y, z = c
    p = np.array([1.0, 0, 0, 0])
    for k in range(z, f.shape[2] - 1):
        r = b.MUL[b.INV[f[x, y, k]], f[x, y, k + 1]]
        q = b.quat_short(r) if r != b.E_ID else np.array([1.0, 0, 0, 0])
        if q is None:
            return float("nan")
        p = b.qmul(p, q)
    return float(np.round(p[0], 6))


class Pair:
    def __init__(self, n=25, sep=4, Tf=0.6, seed=0, rule="lift"):
        self.fr = t10.Frames(n, 0, rule=rule, core=-1, Tf=Tf, seed=seed)
        self.n = n
        L = (n - 1) // 2
        self.L = L
        self.cen = [np.array([L - sep, L, L]), np.array([L + sep, L, L])]
        ax = np.arange(n)
        X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
        self.XYZ = (X, Y, Z)
        self.outer = (X == 0) | (X == n - 1) | (Y == 0) | (Y == n - 1) | (Z == 0) | (Z == n - 1)
        self.refree()

    def block(self, c):
        X, Y, Z = self.XYZ
        return (abs(X - c[0]) <= 1) & (abs(Y - c[1]) <= 1) & (abs(Z - c[2]) <= 1)

    def refree(self):
        core = self.block(self.cen[0]) | self.block(self.cen[1])
        self.fr.free = ~core & ~self.outer
        self.fr.f[core] = b.E_ID

    def twist(self, which, turns, R=6):
        """Spread `turns` full turns about z radially between core `which` (|d|inf <= 1) and |d|inf = R."""
        X, Y, Z = self.XYZ
        c = self.cen[which]
        lay = np.maximum(np.maximum(abs(X - c[0]), abs(Y - c[1])), abs(Z - c[2]))
        for k in range(2, R):
            q = int(round(4 * turns * (R - k) / (R - 2)))
            sel = lay == k
            self.fr.f[sel] = b.rot_z(q)

    def indices(self):
        return [lift_z(self.fr.f, c) for c in self.cen]

    def try_step(self, which, d):
        """Move core `which` by unit vector d if every entering node may be set to the identity."""
        new = self.cen[which] + np.array(d)
        enter = self.block(new) & ~self.block(self.cen[which])
        f = self.fr.f
        trial = f.copy()
        for x, y, z in np.argwhere(enter):
            g = trial[x, y, z]
            if g == b.E_ID:
                continue
            for a, bb, cc in NB:
                fn = trial[x + a, y + bb, z + cc]
                if not self.fr.A[fn, g, b.E_ID]:
                    return False
            trial[x, y, z] = b.E_ID
        self.fr.f = trial
        self.cen[which] = new
        self.refree()
        return True

    def relax(self, sweeps):
        one = np.ones((self.n,) * 3)
        for _ in range(sweeps):
            self.fr.update(one)


def exchange_path(sep, nphi=24, sense=1):
    """Unit-step sequences for both cores along a half circle of radius sep (A starts at -x, B at +x)."""
    pts = [np.array([int(round(-sep * np.cos(p))), int(round(-sense * sep * np.sin(p)))])
           for p in np.linspace(0, np.pi, nphi + 1)]
    steps = []
    for a0, a1 in zip(pts[:-1], pts[1:]):
        d = a1 - a0
        for ax_ in (0, 1):
            for _ in range(abs(d[ax_])):
                u = [0, 0, 0]
                u[ax_] = int(np.sign(d[ax_]))
                steps.append(u)
    return steps


def run(case, sep=4, Tf=0.6, seed=0, relax_step=30, relax_end=4000, rule="lift"):
    os.makedirs(OUT, exist_ok=True)
    p = Pair(sep=sep, Tf=Tf, seed=seed, rule=rule)
    if "twistA" in case:
        p.twist(0, 1)
    if "twistB" in case:
        p.twist(1, 1)
    p.relax(200)
    hist = [dict(stage="start", E=p.fr.energy(), idx=p.indices(), cen=[c.tolist() for c in p.cen])]
    blocked = 0
    t0 = time.time()
    if "exch" in case:
        nex = 2 if "exch2" in case else 1
        for rep in range(nex):
            steps = exchange_path(sep)
            for u in steps:
                for which, sgn in ((0, 1), (1, -1)):
                    d = [sgn * x for x in u]
                    tries = 0
                    while not p.try_step(which, d):
                        p.relax(10)
                        blocked += 1
                        tries += 1
                        if tries > 500:
                            raise RuntimeError("core stuck")
                    p.relax(relax_step)
            hist.append(dict(stage=f"after exchange {rep + 1}", E=p.fr.energy(), idx=p.indices(),
                             cen=[c.tolist() for c in p.cen], blocked=blocked))
    for k in range(relax_end // 500):
        p.relax(500)
        hist.append(dict(stage=f"relax {(k + 1) * 500}", E=p.fr.energy(), idx=p.indices()))
    rec = dict(case=case, sep=sep, Tf=Tf, seed=seed, rule=rule, relax_step=relax_step, relax_end=relax_end,
               hist=hist, wall=time.time() - t0)
    json.dump(rec, open(os.path.join(OUT, f"{case}_{rule}_s{seed}.json"), "w"), indent=1,
              default=lambda o: o.item() if hasattr(o, "item") else o)
    print("RESULT", case, rule, seed, "start", hist[0]["E"], hist[0]["idx"], "| end", hist[-1]["E"], hist[-1]["idx"],
          "blocked", blocked, "wall", round(rec["wall"]), flush=True)
    return rec


if __name__ == "__main__":
    run(sys.argv[1], seed=int(sys.argv[2]) if len(sys.argv) > 2 else 0,
        rule=sys.argv[3] if len(sys.argv) > 3 else "lift")
