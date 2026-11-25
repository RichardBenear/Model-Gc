"""EX-GC36A: causal-invariance audit of CRR Rule 55744 (Claude, 2026-11-25). EXPLORATORY.

Generic multiway engine with relation instances and provenance, a bounded equal-length join test under two conventions
(U unlabelled, A anchored), Stage 1 calibration rules, Stage 2 vacuum forks and Stage 3 particle fork.
Run: python src/exgc36a.py   (writes results/EX-GC36A/results.json)
"""
import json, sys, time, platform, pathlib, collections, itertools
import networkx as nx
from networkx.algorithms.isomorphism import categorical_node_match, categorical_multiedge_match

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "external" / "crr"))
import x6_generic as X6

K = 3          # extension depth after the fork
CAP = 20000    # histories per side per level


# ---------------------------------------------------------------- rules
class Rule:
    """lhs/rhs: tuples of relations over variable names. arity 2 = ordered pair; arity 3 = T(a;{b,c}).
    read: indices of lhs relations that are read, not consumed. injective: distinct variables -> distinct vertices.
    setsem: an output already present is not added again. order: list of (var1,var2) requiring var1 < var2 (by id)."""
    def __init__(self, name, lhs, rhs, arity, read=(), injective=False, setsem=False, order=()):
        self.name, self.lhs, self.rhs, self.arity = name, lhs, rhs, arity
        self.read, self.injective, self.setsem, self.order = set(read), injective, setsem, list(order)
        lv = [v for r in lhs for v in r]
        self.vars = list(dict.fromkeys(lv))
        self.fresh = [v for v in dict.fromkeys(x for r in rhs for x in r) if v not in self.vars]

    def norm(self, r):
        if self.arity == 3 and r[1] > r[2]:
            return (r[0], r[2], r[1])
        return tuple(r)

    def orientations(self, r):
        if self.arity == 3:
            return [r, (r[0], r[2], r[1])] if r[1] != r[2] else [r]
        return [r]


def R55744(version):
    # match: T(x;{y,z}), T(y;{z,x}), T(x;{u,v}) consumed; T(v;{y,z}) read (gate). Roles 1..5 = y,z,u,v,w.
    tpl = [[1, 4, 5], [2, 1, 3], [2, 4, 5], [3, 1, 4], [3, 2, 5], [3, 4, 5], [4, 2, 5], [5, 2, 3]]
    role = {1: "y", 2: "z", 3: "u", 4: "v", 5: "w"}
    rhs = tuple(tuple(role[i] for i in t) for t in tpl)
    lhs = (("x", "y", "z"), ("y", "z", "x"), ("x", "u", "v"), ("v", "y", "z"))
    return Rule("55744" + version, lhs, rhs, 3, read=(3,), injective=True, setsem=True,
                order=[("u", "v")] if version == "L" else [])


# ---------------------------------------------------------------- states
# A history: (rels, events, nextid). rels: tuple of (relation, producer) with producer = event index or -(i+1) for
# initial relation i. events: tuple of (consumed producers..., read producers...) as two tuples.

def make(seed):
    return (tuple((tuple(r), -(i + 1)) for i, r in enumerate(seed)), (), 1 + max(v for r in seed for v in r))


def matches(rule, hist):
    rels, events, nid = hist
    found = {}
    n = len(rule.lhs)
    by_src = collections.defaultdict(list)
    for k, (r, _) in enumerate(rels):
        by_src[r[0]].append(k)

    def rec(i, bind, used):
        if i == n:
            for a, b in rule.order:
                if not bind[a] < bind[b]:
                    return
            if rule.injective and len(set(bind.values())) < len(bind):
                return
            cons = tuple(sorted(used[j] for j in range(n) if j not in rule.read))
            read = tuple(sorted(used[j] for j in range(n) if j in rule.read))
            b2 = dict(bind)
            for j, f in enumerate(rule.fresh):
                b2[f] = nid + j
            out = tuple(sorted(rule.norm(tuple(b2[v] for v in r)) for r in rule.rhs))
            found[(cons, read, out)] = True
            return
        pat = rule.lhs[i]
        cand = by_src.get(bind[pat[0]], ()) if pat[0] in bind else range(len(rels))
        for k in cand:
            r = rels[k][0]
            if k in used:
                continue
            for o in rule.orientations(r):
                b = dict(bind); ok = True
                for var, val in zip(pat, o):
                    if b.setdefault(var, val) != val:
                        ok = False; break
                if ok:
                    rec(i + 1, b, used + [k])

    rec(0, {}, [])
    return list(found)


def apply(rule, hist, m):
    rels, events, nid = hist
    cons, read, out = m
    e = len(events)
    keep = [x for k, x in enumerate(rels) if k not in cons]
    present = {r for r, _ in keep}
    for r in out:
        if rule.setsem and r in present:
            continue
        keep.append((r, e)); present.add(r)
    ev = (tuple(rels[k][1] for k in cons), tuple(rels[k][1] for k in read))
    return (tuple(keep), events + (ev,), nid + len(rule.fresh))


def graph(rule, hist, seedverts, conv, base_events=0):
    """Combined graph of events after the fork root and live relations. conv 'U' or 'A'."""
    rels, events, _ = hist
    g = nx.MultiDiGraph()
    def pnode(p):
        if p >= base_events:
            return ("e", p)
        return None   # produced before (or at) the root: treated like an initial relation
    for i in range(base_events, len(events)):
        g.add_node(("e", i), lab="ev")
        cons, read = events[i]
        for kind, ps in (("cons", cons), ("read", read)):
            for p in ps:
                if p >= base_events:
                    g.add_edge(("e", p), ("e", i), kind=kind)
                elif conv == "A":
                    root = ("i", p); g.add_node(root, lab=("i", p)); g.add_edge(root, ("e", i), kind=kind)
    for j, (r, p) in enumerate(rels):
        rn = ("r", j); g.add_node(rn, lab="rel")
        names = ("src", "tgt", "tgt") if rule.arity == 3 else ("tail", "head")
        for nm, v in zip(names, r):
            vn = ("v", v)
            g.add_node(vn, lab=("v", v) if (conv == "A" and v in seedverts) else "v")
            g.add_edge(rn, vn, kind=nm)
        if p >= base_events:
            g.add_edge(("e", p), rn, kind="prod")
        elif conv == "A":
            root = ("i", p); g.add_node(root, lab=("i", p)); g.add_edge(root, rn, kind="prod")
    return g


def iso(a, b):
    return nx.is_isomorphic(a, b, node_match=categorical_node_match("lab", None),
                            edge_match=categorical_multiedge_match("kind", None))


def whash(g, rounds=3):
    """Colour-refinement invariant for a labelled multigraph (a prefilter only; equality is decided by iso)."""
    col = {n: str(d["lab"][0] if isinstance(d["lab"], tuple) else d["lab"]) if False else repr(d["lab"]) for n, d in g.nodes(data=True)}
    for _ in range(rounds):
        new = {}
        for n in g.nodes:
            outs = sorted((d["kind"], col[b]) for _, b, d in g.out_edges(n, data=True))
            ins = sorted((d["kind"], col[a]) for a, _, d in g.in_edges(n, data=True))
            new[n] = repr((col[n], outs, ins))
        col = {n: str(hash(c)) for n, c in new.items()}
    return hash(tuple(sorted(col.values())))


def relabel_root(hist):
    """Make the fork root's live relations count as 'initial' (negative producers) so graphs start at the root."""
    rels, events, nid = hist
    return (tuple((r, -(k + 1)) for k, (r, _) in enumerate(rels)), (), nid)


def join(rule, root, m1, m2, conv, k=K):
    root = relabel_root(root)
    seedverts = {v for r, _ in root[0] for v in r}
    L = [apply(rule, root, m1)]; R = [apply(rule, root, m2)]
    for level in range(0, k + 1):
        buckets = collections.defaultdict(list)
        for h in L:
            g = graph(rule, h, seedverts, conv); buckets[whash(g)].append(g)
        for h in R:
            g = graph(rule, h, seedverts, conv); key = whash(g)
            if key in buckets and any(iso(g, x) for x in buckets[key]):
                return ("join", level)
        if level == k:
            break
        L2 = [apply(rule, h, m) for h in L for m in matches(rule, h)]
        R2 = [apply(rule, h, m) for h in R for m in matches(rule, h)]
        if not L2 and not R2:
            return ("terminal-distinct", level)
        if not L2 or not R2:
            return ("terminal-distinct", level)
        if max(len(L2), len(R2)) > CAP:
            return ("cap", level)
        L, R = L2, R2
    return ("unresolved", k)


def forks(rule, hist):
    ms = matches(rule, hist)
    return ms, [(a, b) for a, b in itertools.combinations(ms, 2)]


def shared(a, b):
    return bool(set(a[0]) & (set(b[0]) | set(b[1]))) or bool(set(b[0]) & set(a[1]))


# ---------------------------------------------------------------- stage 1
def stage1():
    out = []
    pos1 = Rule("xy->xy,yz", (("x", "y"),), (("x", "y"), ("y", "z")), 2)
    pos2 = Rule("xy,zy->xw,yw,zw", (("x", "y"), ("z", "y")), (("x", "w"), ("y", "w"), ("z", "w")), 2)
    neg = Rule("xy,xz->xy", (("x", "y"), ("x", "z")), (("x", "y"),), 2)
    cases = [(pos1, [(0, 0), (0, 0)], "two loops"),
             (pos2, [(0, 3), (1, 3), (2, 3), (4, 0)], "3-in star + leaf edge"),
             (pos2, [(0, 4), (1, 4), (2, 4), (3, 4)], "4-in star"),
             (neg, [(0, 1), (0, 2), (2, 3)], "2-out star + leaf edge")]
    for rule, seed, label in cases:
        h = make(seed); ms, fs = forks(rule, h)
        rows = []
        for a, b in fs:
            rows.append(dict(shared=shared(a, b), U=join(rule, h, a, b, "U"), A=join(rule, h, a, b, "A")))
        out.append(dict(rule=rule.name, seed=seed, label=label, matches=len(ms), forks=rows))
    return out


# ---------------------------------------------------------------- stage 2
START = [(0, 1, 2), (1, 2, 0), (0, 3, 4), (4, 1, 2)]


def stage2a(version, depth=4):
    rule = R55744(version)
    level = [make([rule.norm(r) for r in START])]
    rows = []; seen = set()
    for d in range(depth + 1):
        nxt = []
        for h in level:
            ms, fs = forks(rule, h)
            key = (tuple(sorted(r for r, _ in h[0])), h[2])
            if key not in seen:
                seen.add(key)
                for a, b in fs:
                    rows.append(dict(depth=d, shared=shared(a, b), U=join(rule, h, a, b, "U"),
                                     A=join(rule, h, a, b, "A")))
            nxt += [apply(rule, h, m) for m in ms]
        # merge exact duplicates (same relations and counter) to keep the tree finite
        uniq = {}
        for h in nxt:
            uniq.setdefault((tuple(sorted(r for r, _ in h[0])), h[2]), h)
        level = list(uniq.values())
    return rows


def stage2b(version, events=2000, npairs=50):
    T = [[1, 4, 5], [2, 1, 3], [2, 4, 5], [3, 1, 4], [3, 2, 5], [3, 4, 5], [4, 2, 5], [5, 2, 3]]
    sim = X6.Sim(X6.RuleSpec(55744, T, []), 0); sim.run(events)
    gates = sorted(sim.gates_set)
    cons = {g: {X6.rel(g[0], g[1], g[2]), X6.rel(g[1], g[2], g[0]), X6.rel(g[0], g[3], g[4])} for g in gates}
    gate = {g: X6.rel(g[4], g[1], g[2]) for g in gates}
    sh, rw = [], []
    for a, b in itertools.combinations(gates, 2):
        if cons[a] & cons[b]:
            if len(sh) < npairs: sh.append((a, b))
        elif gate[a] in cons[b] or gate[b] in cons[a]:
            if len(rw) < npairs: rw.append((a, b))
        if len(sh) >= npairs and len(rw) >= npairs:
            break
    A2 = sim.twosection()
    rule = R55744(version)
    rows = []
    for kind, pairs in (("shared", sh), ("readwrite", rw)):
        for a, b in pairs:
            core = set(a) | set(b)
            ball = set(core)
            for v in core:
                ball |= A2.get(v, set())
            ext = sorted(r for r in sim.rels if set(r) <= ball)
            if len(ext) > 400:
                rows.append(dict(kind=kind, pair=[a, b], extract=len(ext), note="skipped: extract > 400 (D2)")); continue
            h = make(ext)
            h = (h[0], h[1], sim.n)     # fresh ids continue the simulator counter
            ms = matches(rule, h)
            def find(g):
                x, y, z, u, v = g
                for m in ms:
                    c = {h[0][k][0] for k in m[0]}
                    if c == {X6.rel(x, y, z), X6.rel(y, z, x), X6.rel(x, u, v)} and h[0][m[1][0]][0] == X6.rel(v, y, z):
                        return m
                return None
            ma, mb = find(a), find(b)
            if ma is None or mb is None:
                rows.append(dict(kind=kind, pair=[a, b], note="not a match in this version")); continue
            rows.append(dict(kind=kind, pair=[a, b], extract=len(ext), matches=len(ms),
                             U=join(rule, h, ma, mb, "U"), A=join(rule, h, ma, mb, "A")))
    return rows


# ---------------------------------------------------------------- stage 3
def stage3(version):
    spec = X6.load_specs(str(ROOT / "external" / "crr" / "W21_CLOSED_PARTICLE_HITS.json"))[55744]
    B = sorted(spec.phase_rels[1])
    rule = R55744(version)
    h = make(B)
    ms, fs = forks(rule, h)
    closure = []
    for m in ms:
        h2 = apply(rule, h, m)
        cons, read, out = m
        w = h[2]
        rels = {r for r, _ in h2[0]}
        phases = []
        for sub in itertools.combinations(sorted({v for r in rels for v in r}), 5):
            rr = X6.normalize_subset(rels, sub)
            for q, p in enumerate(spec.phase_rels):
                if rr == p and w in sub:
                    phases.append((list(sub), q))
        closure.append(dict(match=[[list(h[0][k][0]) for k in cons], list(h[0][read[0]][0])], successor_phases=phases))
    rows = [dict(shared=shared(a, b), U=join(rule, h, a, b, "U"), A=join(rule, h, a, b, "A")) for a, b in fs]
    return dict(version=version, matches=len(ms), closure=closure, forks=rows)


def main():
    t0 = time.process_time()
    res = dict(platform=dict(python=platform.python_version(), networkx=nx.__version__, machine=platform.machine()),
               K=K, CAP=CAP)
    res["stage1"] = stage1(); print("stage1 done", time.process_time() - t0, flush=True)
    res["stage3"] = [stage3(v) for v in "LC"]; print("stage3 done", time.process_time() - t0, flush=True)
    res["stage2a"] = {v: stage2a(v) for v in "LC"}; print("stage2a done", time.process_time() - t0, flush=True)
    res["stage2b"] = {v: stage2b(v) for v in "LC"}; print("stage2b done", time.process_time() - t0, flush=True)
    res["cpu_seconds"] = time.process_time() - t0
    out = ROOT / "results" / "EX-GC36A"; out.mkdir(parents=True, exist_ok=True)
    (out / "results.json").write_text(json.dumps(res, indent=1, default=list) + "\n")
    print("cpu", res["cpu_seconds"])


if __name__ == "__main__":
    main()
