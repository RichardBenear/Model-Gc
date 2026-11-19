"""EX-GC08 (exploratory): route 3, stage 1 - a discrete belt trick on a cubic ribbon graph.

On a cubic lattice a node's rotation system (which incident edge is its +x, +y, +z) is one of the 24 rotations of
the cube group O. "No-cut" rule: neighbouring frames may differ only by a step in S (a union of conjugacy classes).
Part A: topology of the move space = flag complex of the Cayley graph Cay(O, S); is the 360 degree loop nontrivial
(H1 over GF(2)), and is the 720 degree loop a boundary?  Part B: Metropolis annealing of a twisted shell.
"""
import itertools
import json
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC08")

# ---------------------------------------------------------------- the cube group O
G = []
for p in itertools.permutations(range(3)):
    for s in itertools.product((1, -1), repeat=3):
        m = np.zeros((3, 3), int)
        for i in range(3):
            m[i, p[i]] = s[i]
        if round(np.linalg.det(m)) == 1:
            G.append(m)
N = len(G)                                   # 24
KEY = {m.tobytes(): i for i, m in enumerate(G)}
MUL = np.array([[KEY[(G[a] @ G[b]).tobytes()] for b in range(N)] for a in range(N)])
INV = np.array([KEY[G[a].T.copy().tobytes()] for a in range(N)])
E_ID = KEY[np.eye(3, dtype=int).tobytes()]


def cls(i):
    m = G[i]
    tr = int(np.trace(m))
    if tr == 3:
        return "id"
    if tr == 1:
        return "90"
    if tr == 0:
        return "120"
    return "180f" if np.count_nonzero(m - np.diag(np.diag(m))) == 0 else "180e"


CLASSES = ["90", "180f", "120", "180e"]


def rot_z(q):
    c = {0: (1, 0), 1: (0, 1), 2: (-1, 0), 3: (0, -1)}[q % 4]
    return KEY[np.array([[c[0], -c[1], 0], [c[1], c[0], 0], [0, 0, 1]]).tobytes()]


# ---------------------------------------------------------------- Part A: flag complex homology
def gf2_rank(M):
    M = M.copy() % 2
    r = 0
    rows, cols = M.shape
    for c in range(cols):
        piv = np.nonzero(M[r:, c])[0]
        if len(piv) == 0:
            continue
        p = r + piv[0]
        M[[r, p]] = M[[p, r]]
        others = np.nonzero(M[:, c])[0]
        others = others[others != r]
        M[others] ^= M[r]
        r += 1
        if r == rows:
            break
    return r


def complex_for(S):
    Sset = set(S)
    edges = sorted({tuple(sorted((g, MUL[g, s]))) for g in range(N) for s in S})
    eidx = {e: i for i, e in enumerate(edges)}
    adj = {g: {MUL[g, s] for s in S} for g in range(N)}
    tris = sorted({tuple(sorted((a, b, c))) for a in range(N) for b in adj[a] for c in adj[a] & adj[b]})
    d1 = np.zeros((N, len(edges)), np.uint8)
    for i, (a, b) in enumerate(edges):
        d1[a, i] = d1[b, i] = 1
    d2 = np.zeros((len(edges), len(tris)), np.uint8)
    for j, (a, b, c) in enumerate(tris):
        for e in ((a, b), (a, c), (b, c)):
            d2[eidx[e], j] = 1
    return edges, eidx, tris, d1, d2, Sset


def loop_vec(path, eidx, nE):
    v = np.zeros(nE, np.uint8)
    for a, b in zip(path[:-1], path[1:]):
        v[eidx[tuple(sorted((a, b)))]] ^= 1
    return v


def turn_path(S, turns):
    """A closed path of 'turns' full rotations about one axis using steps in S (None if S cannot do it)."""
    cs = {cls(s) for s in S}
    if "90" in cs:
        step = [rot_z(1)]
    elif "120" in cs:
        step = [next(s for s in S if cls(s) == "120")]
    elif "180f" in cs:
        step = [rot_z(2)]
    else:
        return None
    s = step[0]
    order = next(k for k in range(1, 7) if np.array_equal(np.linalg.matrix_power(G[s], k), np.eye(3, dtype=int)))
    path, g = [E_ID], E_ID
    for _ in range(order * turns):
        g = MUL[g, s]
        path.append(g)
    return path


def part_a():
    rows = []
    for k in range(1, 5):
        for combo in itertools.combinations(CLASSES, k):
            S = [i for i in range(N) if cls(i) in combo]
            edges, eidx, tris, d1, d2, _ = complex_for(S)
            r1, r2 = gf2_rank(d1), gf2_rank(d2)
            h1_gf2 = len(edges) - r1 - r2
            r1q = np.linalg.matrix_rank(d1.astype(float) * np.where(np.arange(N)[:, None] == np.array([a for a, b in edges])[None, :], -1, 1))
            # rational H1 with oriented boundary maps
            D1 = np.zeros((N, len(edges)))
            for i, (a, b) in enumerate(edges):
                D1[a, i], D1[b, i] = -1, 1
            D2 = np.zeros((len(edges), len(tris)))
            for j, (a, b, c) in enumerate(tris):
                D2[eidx[(a, b)], j] += 1
                D2[eidx[(b, c)], j] += 1
                D2[eidx[(a, c)], j] -= 1
            h1_q = len(edges) - np.linalg.matrix_rank(D1) - (np.linalg.matrix_rank(D2) if len(tris) else 0)
            res = dict(S="+".join(combo), n_steps=len(S), edges=len(edges), triangles=len(tris), H1_GF2=int(h1_gf2),
                       H1_Q=int(h1_q))
            for turns in (1, 2):
                p = turn_path(S, turns)
                if p is None:
                    res[f"loop{turns}"] = None
                    continue
                c = loop_vec(p, eidx, len(edges))
                if not c.any():
                    res[f"loop{turns}"] = "trivial (backtrack)"
                    continue
                bd = gf2_rank(np.column_stack([d2, c])) == r2 if len(tris) else not c.any()
                res[f"loop{turns}"] = "boundary" if bd else "NOT boundary"
            rows.append(res)
            print(res, flush=True)
    return rows


# ---------------------------------------------------------------- Part B: annealing a twisted shell
def part_b(Sclasses, turns, L=10, sweeps=1500, T0=2.0, seed=0):
    S = [i for i in range(N) if cls(i) in Sclasses]
    ok = np.zeros((N, N), bool)                   # ok[g, h]: relative rotation g^-1 h allowed
    for g in range(N):
        for h in range(N):
            r = MUL[INV[g], h]
            ok[g, h] = r == E_ID or r in S
    mism = np.array([[MUL[INV[g], h] != E_ID for h in range(N)] for g in range(N)], int)
    n = 2 * L + 1
    ax = np.arange(n) - L
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    layer = np.maximum(np.maximum(abs(X), abs(Y)), abs(Z))
    # initial twist: layer k in 1..L, angle steps of 90 degrees spread over layers 1..L
    steps = 4 * turns
    if steps > L - 1:
        raise ValueError("not enough layers")
    q = np.zeros_like(layer)
    for k in range(1, L + 1):
        q[layer == k] = int(round(steps * (L - k) / (L - 1)))     # core side rotated by 'turns', boundary 0
    frame = np.vectorize(rot_z)(q)
    frame[layer <= 1] = rot_z(0)                   # core held at identity after the twist (full turns)
    frame[layer == L] = E_ID
    free = np.argwhere((layer > 1) & (layer < L))
    nb = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

    def energy():
        e = 0
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            a = frame[:n - d[0], :n - d[1], :n - d[2]]
            b = frame[d[0]:, d[1]:, d[2]:]
            e += mism[a, b].sum()
        return int(e)

    def valid_all():
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            a = frame[:n - d[0], :n - d[1], :n - d[2]]
            b = frame[d[0]:, d[1]:, d[2]:]
            if not ok[a, b].all():
                return False
        return True

    assert valid_all(), "initial twist violates the step rule"
    rng = np.random.default_rng(seed)
    E = energy()
    hist = [dict(sweep=0, T=T0, E=E)]
    t0 = time.time()
    for sw in range(1, sweeps + 1):
        T = T0 * (1 - sw / sweeps) + 1e-3
        order = rng.permutation(len(free))
        props = rng.integers(0, N, len(free))
        us = rng.random(len(free))
        for j, idx in enumerate(order):
            x, y, z = free[idx]
            g = frame[x, y, z]
            h = props[j]
            if h == g:
                continue
            nbr = [frame[x + a, y + b, z + c] for a, b, c in nb]
            if not all(ok[h, f] for f in nbr):
                continue
            dE = sum(mism[h, f] - mism[g, f] for f in nbr)
            if dE <= 0 or us[j] < np.exp(-dE / T):
                frame[x, y, z] = h
                E += dE
        if sw % 50 == 0:
            hist.append(dict(sweep=sw, T=T, E=E))
            print(Sclasses, turns, "sweep", sw, "T %.3f E %d" % (T, E), "wall %.0f" % (time.time() - t0), flush=True)
    assert valid_all() and E == energy()
    return dict(S="+".join(Sclasses), turns=turns, L=L, sweeps=sweeps, E_init=hist[0]["E"], E_final=E, hist=hist,
                frame=frame)


def z2_index(frame, L, S):
    """Product of short SU(2) lifts along the +x radial path from the core (layer 1) to the boundary."""
    def quat(i):
        m = G[i].astype(float)
        w = np.sqrt(max(0, 1 + np.trace(m))) / 2
        if w > 1e-9:
            v = np.array([m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1]]) / (4 * w)
            return np.array([w, *v])
        ax_ = np.sqrt(np.maximum(0, (np.diag(m) + 1) / 2))
        for i2 in range(3):                        # fix relative signs from off-diagonals
            for j2 in range(i2 + 1, 3):
                if m[i2, j2] < 0 and ax_[i2] > 0 and ax_[j2] > 0:
                    ax_[j2] = -ax_[j2]
        return np.array([0.0, *ax_])

    def qmul(a, b):
        w1, x1, y1, z1 = a
        w2, x2, y2, z2 = b
        return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                         w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])
    c = L
    p = np.array([1.0, 0, 0, 0])
    for k in range(1, L):
        g, h = frame[c + k, c, c], frame[c + k + 1, c, c]
        p = qmul(p, quat(MUL[INV[g], h]))
    return float(p[0])


if __name__ == "__main__" and sys.argv[1] != "B2":
    os.makedirs(OUT, exist_ok=True)
    if sys.argv[1] == "A":
        json.dump(part_a(), open(os.path.join(OUT, "partA.json"), "w"), indent=1)
    else:
        Scl = sys.argv[2].split("+")
        turns = int(sys.argv[3])
        sweeps = int(sys.argv[4]) if len(sys.argv) > 4 else 1500
        seed = int(sys.argv[5]) if len(sys.argv) > 5 else 0
        r = part_b(Scl, turns, sweeps=sweeps, seed=seed)
        S = [i for i in range(N) if cls(i) in Scl]
        r["z2_final"] = z2_index(r["frame"], r["L"], S)
        fr = r.pop("frame")
        tag = f"B_{sys.argv[2]}_t{turns}_s{seed}"
        json.dump(r, open(os.path.join(OUT, f"{tag}.json"), "w"), indent=1)
        print("RESULT", tag, "E", r["E_init"], "->", r["E_final"], "z2", r["z2_final"], flush=True)


# ---------------------------------------------------------------- Part A2: lift-consistent ("no-cut") 2-cells
def quat_short(i):
    m = G[i].astype(float)
    w = np.sqrt(max(0.0, 1 + np.trace(m))) / 2
    if w < 1e-9:
        return None                                # 180 degree: no unique short lift
    v = np.array([m[2, 1] - m[1, 2], m[0, 2] - m[2, 0], m[1, 0] - m[0, 1]]) / (4 * w)
    return np.array([w, *v])


def qmul(a, b):
    w1, x1, y1, z1 = a
    w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def consistent(f, g, h):
    """Triangle f-g-h is 'small': short lifts compose, q(f^-1 g) q(g^-1 h) = +q(f^-1 h)."""
    a, b, c = quat_short(MUL[INV[f], g]), quat_short(MUL[INV[g], h]), quat_short(MUL[INV[f], h])
    if a is None or b is None or c is None:
        return False
    return np.allclose(qmul(a, b), c)


def part_a2(Scl=("90", "120")):
    S = [i for i in range(N) if cls(i) in Scl]
    edges, eidx, tris_all, d1, d2_all, _ = complex_for(S)
    tris = [t for t in tris_all if consistent(*t)]
    D1 = np.zeros((N, len(edges)))
    for i, (a, b) in enumerate(edges):
        D1[a, i], D1[b, i] = -1, 1
    D2 = np.zeros((len(edges), len(tris)))
    d2 = np.zeros((len(edges), len(tris)), np.uint8)
    for j, (a, b, c) in enumerate(tris):
        D2[eidx[(a, b)], j] += 1
        D2[eidx[(b, c)], j] += 1
        D2[eidx[(a, c)], j] -= 1
        for e_ in ((a, b), (a, c), (b, c)):
            d2[eidx[e_], j] = 1
    r2 = gf2_rank(d2)
    h1g = len(edges) - gf2_rank(d1) - r2
    h1q = len(edges) - np.linalg.matrix_rank(D1) - np.linalg.matrix_rank(D2)
    p1 = turn_path(S, 1)
    c = loop_vec(p1, eidx, len(edges))
    res = dict(S="+".join(Scl), edges=len(edges), triangles_all=len(tris_all), triangles_consistent=len(tris),
               H1_GF2=int(h1g), H1_Q=int(h1q), loop1="boundary" if gf2_rank(np.column_stack([d2, c])) == r2 else
               "NOT boundary")
    # product of lifts around the 1-turn and 2-turn loops (the Z2 index)
    for turns in (1, 2):
        p = turn_path(S, turns)
        q = np.array([1.0, 0, 0, 0])
        for a, b in zip(p[:-1], p[1:]):
            q = qmul(q, quat_short(MUL[INV[a], b]))
        res[f"lift_product_{turns}turn"] = float(q[0])
    print(res, flush=True)
    return res


def move_table(Scl, rule):
    """allowed[f, g, h]: node at g may move to h next to a neighbour at f.
    rule 'flag': only f^-1 h in S+e (pairwise no-cut).  rule 'lift': also the triangle (f, g, h) is lift-consistent."""
    S = [i for i in range(N) if cls(i) in Scl]
    Se = set(S) | {E_ID}
    A = np.zeros((N, N, N), bool)
    for f in range(N):
        for g in range(N):
            for h in range(N):
                if MUL[INV[f], h] not in Se or MUL[INV[g], h] not in Se:
                    continue
                if rule == "lift" and len({f, g, h}) == 3 and not consistent(f, g, h):
                    continue
                if rule == "lift" and len({f, g, h}) == 2 and f != g and g != h and f == h:
                    pass                      # moving back onto the neighbour's frame: always small
                A[f, g, h] = True
    return S, A


def part_b2(Scl, rule, turns, L=10, sweeps=1500, T0=1.0, seed=0):
    S, A = move_table(Scl, rule)
    steps_out = [MUL[g, s] for g in range(N) for s in S]
    nbrs = np.array([[MUL[g, s] for s in S] for g in range(N)])
    mism = np.array([[MUL[INV[g], h] != E_ID for h in range(N)] for g in range(N)], int)
    n = 2 * L + 1
    ax_ = np.arange(n) - L
    X, Y, Z = np.meshgrid(ax_, ax_, ax_, indexing="ij")
    layer = np.maximum(np.maximum(abs(X), abs(Y)), abs(Z))
    q = np.zeros_like(layer)
    for k in range(1, L + 1):
        q[layer == k] = int(round(4 * turns * (L - k) / (L - 1)))
    frame = np.vectorize(rot_z)(q)
    frame[layer <= 1] = E_ID
    frame[layer == L] = E_ID
    free = np.argwhere((layer > 1) & (layer < L))
    nb = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]

    def energy():
        e_ = 0
        for d in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
            e_ += mism[frame[:n - d[0], :n - d[1], :n - d[2]], frame[d[0]:, d[1]:, d[2]:]].sum()
        return int(e_)
    rng = np.random.default_rng(seed)
    E = energy()
    z0 = z2_lift(frame, L)
    hist = [dict(sweep=0, T=T0, E=E, z2=z0)]
    t0 = time.time()
    for sw in range(1, sweeps + 1):
        T = T0 * (1 - sw / sweeps) + 1e-3
        order = rng.permutation(len(free))
        pick = rng.integers(0, len(S), len(free))
        us = rng.random(len(free))
        for j, idx in enumerate(order):
            x, y, z = free[idx]
            g = frame[x, y, z]
            h = nbrs[g, pick[j]]
            nbr = [frame[x + a, y + b, z + c] for a, b, c in nb]
            if not all(A[f, g, h] for f in nbr):
                continue
            dE = sum(mism[h, f] - mism[g, f] for f in nbr)
            if dE <= 0 or us[j] < np.exp(-dE / T):
                frame[x, y, z] = h
                E += dE
        if sw % 100 == 0:
            z = z2_lift(frame, L)
            hist.append(dict(sweep=sw, T=T, E=E, z2=z))
            print(Scl, rule, turns, seed, "sweep", sw, "T %.3f E %d z2 %+.0f" % (T, E, z), "wall %.0f" % (
                time.time() - t0), flush=True)
    assert E == energy()
    return dict(S="+".join(Scl), rule=rule, turns=turns, L=L, sweeps=sweeps, seed=seed, E_init=hist[0]["E"],
                E_final=E, z2_init=z0, z2_final=z2_lift(frame, L), hist=hist, wall=time.time() - t0)


def z2_lift(frame, L):
    """Product of short lifts along the +x radial path from the core to the boundary (sign = Z2 index)."""
    c = L
    p = np.array([1.0, 0, 0, 0])
    for k in range(0, L):
        r = MUL[INV[frame[c + k, c, c]], frame[c + k + 1, c, c]]
        qq = quat_short(r) if r != E_ID else np.array([1.0, 0, 0, 0])
        if qq is None:
            return float("nan")
        p = qmul(p, qq)
    return float(np.round(p[0], 6))


if __name__ == "__main__" and sys.argv[1] == "B2":
    os.makedirs(OUT, exist_ok=True)
    Scl, rule, turns = sys.argv[2].split("+"), sys.argv[3], int(sys.argv[4])
    sweeps = int(sys.argv[5]) if len(sys.argv) > 5 else 1500
    seed = int(sys.argv[6]) if len(sys.argv) > 6 else 0
    r = part_b2(Scl, rule, turns, sweeps=sweeps, seed=seed)
    tag = f"B2_{sys.argv[2]}_{rule}_t{turns}_s{seed}"
    json.dump(r, open(os.path.join(OUT, f"{tag}.json"), "w"), indent=1, default=lambda o: o.item())
    print("RESULT", tag, "E", r["E_init"], "->", r["E_final"], "z2", r["z2_init"], "->", r["z2_final"], flush=True)
