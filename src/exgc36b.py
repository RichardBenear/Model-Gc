"""EX-GC36B (Claude, 2026-11-25, EXPLORATORY). See spec/EX-GC36B_plan.md.
Step 1: order-free closure screen of the 80 CRR closed-particle rules. Step 2: the four odd 36A vacuum forks.
Step 3: state-frontier delayed-join comparator F, calibrated on typed-chain witnesses, then applied.
Run: python src/exgc36b.py  (writes results/EX-GC36B/results.json)"""
import json, sys, time, platform, pathlib, itertools, collections
import networkx as nx
from networkx.algorithms.isomorphism import categorical_node_match, categorical_multiedge_match

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src")); sys.path.insert(0, str(ROOT / "external" / "crr"))
import exgc36a as A

D = 6
CAP = 20000
PERMS = list(itertools.permutations(range(5)))


def rel(a, b, c):
    return (a, b, c) if b < c else (a, c, b)


def decode(m):
    U = [(a, b, c) for a in range(5) for b, c in itertools.combinations([i for i in range(5) if i != a], 2)]
    return {U[i] for i in range(30) if (int(m) >> i) & 1}


def canon5(R, S):
    """Canonical form of the induced relations of R on the 5-vertex support S (brute force over 120 bijections)."""
    S = list(S); idx = {v: i for i, v in enumerate(S)}
    sub = [(idx[a], idx[b], idx[c]) for a, b, c in R if a in idx and b in idx and c in idx]
    best = None
    for p in PERMS:
        f = tuple(sorted(rel(p[a], p[b], p[c]) for a, b, c in sub))
        if best is None or f < best:
            best = f
    return best


def moves(R, ordered):
    """All x6-pattern matches in relation set R. ordered=True: u<v by vertex number (version L); else both (C)."""
    by = collections.defaultdict(set)
    for r in R:
        by[r[0]].add(r)
    out = set()
    for x, xr in by.items():
        for r1 in xr:
            _, p, q = r1
            for y, z in ((p, q), (q, p)):
                if len({x, y, z}) < 3 or rel(y, z, x) not in R:
                    continue
                for r3 in xr:
                    if r3 == r1:
                        continue
                    _, a, b = r3
                    for u, v in ((a, b),) if ordered else ((a, b), (b, a)):
                        if len({x, y, z, u, v}) < 5:
                            continue
                        if rel(v, y, z) in R:
                            out.add((x, y, z, u, v))
    return sorted(out)


def fire(R, g, tpl, w):
    x, y, z, u, v = g
    R2 = set(R) - {rel(x, y, z), rel(y, z, x), rel(x, u, v)}
    roles = {1: y, 2: z, 3: u, 4: v, 5: w}
    for a, b, c in tpl:
        R2.add(rel(roles[a], roles[b], roles[c]))
    return R2, (y, z, u, v, w)


def step1():
    data = json.load(open(ROOT / "external" / "crr" / "W21_CLOSED_PARTICLE_HITS.json"))
    rows = []
    for e in data:
        tpl = e["template"]
        phases = [decode(m) for m in e["state_masks"]]
        classes = {canon5(P, range(5)) for P in phases}
        prow = []
        for i, P in enumerate(phases):
            # (a) L order dependence
            succ_sets = collections.Counter()
            for p in PERMS:
                Q = {rel(p[a], p[b], p[c]) for a, b, c in P}
                s = frozenset(canon5(fire(Q, g, tpl, 5)[0], fire(Q, g, tpl, 5)[1]) for g in moves(Q, True))
                succ_sets[s] += 1
            # (b) C closure
            cm = moves(P, False)
            land = [canon5(*fire(P, g, tpl, 5)) in classes for g in cm]
            lm = moves(P, True)
            land_L = [canon5(*fire(P, g, tpl, 5)) in classes for g in lm]
            prow.append(dict(phase=i, L_moves=len(lm), L_closed=all(land_L) and len(lm) > 0,
                             L_distinct_successor_sets=len(succ_sets),
                             L_successor_set_sizes=sorted(collections.Counter(len(s) for s in succ_sets.elements()).items()),
                             C_moves=len(cm), C_in_class=sum(land), C_closed=len(cm) > 0 and all(land)))
        rows.append(dict(rule=e["template_index"], period=e["period"], phases=prow,
                         passes=all(r["C_closed"] for r in prow),
                         order_free_L=all(r["L_distinct_successor_sets"] == 1 for r in prow),
                         L_closed_published=all(r["L_closed"] for r in prow)))
    return rows


# ---------------------------------------------------------------- step 2
def step2():
    rule = A.R55744("L")
    level = [A.make([rule.norm(r) for r in A.START])]
    seen = set(); out = []
    for d in range(5):
        nxt = []
        for h in level:
            ms = A.matches(rule, h)
            key = (tuple(sorted(r for r, _ in h[0])), h[2])
            if key not in seen:
                seen.add(key)
                for a, b in itertools.combinations(ms, 2):
                    if A.shared(a, b):
                        continue
                    if A.join(rule, h, a, b, "U")[0] == "join":
                        continue
                    out.append(diagnose(rule, h, a, b, d))
            nxt += [A.apply(rule, h, m) for m in ms]
        uniq = {}
        for h in nxt:
            uniq.setdefault((tuple(sorted(r for r, _ in h[0])), h[2]), h)
        level = list(uniq.values())
    return out


def diagnose(rule, h, a, b, depth):
    rels = [r for r, _ in h[0]]
    def desc(m):
        return dict(consumed=[rels[k] for k in m[0]], read=[rels[k] for k in m[1]], outputs=list(m[2]))
    res = dict(depth=depth, e1=desc(a), e2=desc(b))
    ha, hb = A.apply(rule, h, a), A.apply(rule, h, b)
    def find(hh, m):
        want_c = {rels[k] for k in m[0]}; want_r = {rels[k] for k in m[1]}
        for mm in A.matches(rule, hh):
            if {hh[0][k][0] for k in mm[0]} == want_c and {hh[0][k][0] for k in mm[1]} == want_r:
                return mm
        return None
    b_after_a, a_after_b = find(ha, b), find(hb, a)
    res["second_legal"] = dict(b_after_a=b_after_a is not None, a_after_b=a_after_b is not None)
    dup = {}
    for name, m, base in (("e1", a, h), ("e2", b, h)):
        present = {r for r, _ in base[0]} - {rels[k] for k in m[0]}
        dup[name] = [r for r in m[2] if r in present]
    res["outputs_already_present"] = dup
    if b_after_a is not None and a_after_b is not None:
        x, y = A.apply(rule, ha, b_after_a), A.apply(rule, hb, a_after_b)
        res["same_live_set_after_both_orders"] = sorted(r for r, _ in x[0]) == sorted(r for r, _ in y[0])
        res["live_counts"] = [len(x[0]), len(y[0])]
    cause = "iv"
    if b_after_a is None or a_after_b is None:
        # order sensitivity through numbering? test eligibility of the same pattern with renumbered fresh vertex
        cause = "i"
    elif dup["e1"] or dup["e2"]:
        cause = "ii"
    res["class"] = cause
    return res


# ---------------------------------------------------------------- step 3: comparator F
def hyper_graph(rels, typed=False):
    g = nx.MultiDiGraph()
    for j, r in enumerate(sorted(rels)):
        rn = ("r", j)
        if typed:
            g.add_node(rn, lab=r[0]); g.add_node(("v", r[1]), lab="v"); g.add_edge(rn, ("v", r[1]), kind="at")
            continue
        g.add_node(rn, lab="r")
        names = ("src", "tgt", "tgt") if len(r) == 3 else ("tail", "head")
        for nm, v in zip(names, r):
            g.add_node(("v", v), lab="v"); g.add_edge(rn, ("v", v), kind=nm)
    return g


def inv(g):
    return (g.number_of_nodes(), g.number_of_edges(), A.whash(g))


class Frontier:
    def __init__(self):
        self.b = collections.defaultdict(list)
    def add(self, g):
        k = inv(g)
        for x in self.b[k]:
            if A.iso(x, g):
                return False
        self.b[k].append(g); return True
    def has(self, g):
        return any(A.iso(x, g) for x in self.b.get(inv(g), ()))


def frontier_join(step_fn, graph_fn, s1, s2, depth=D, cap=CAP):
    """BFS each side with isomorphism dedupe; join if any class at any depth pair is shared. Returns
    ('join', d1, d2) with smallest d1+d2, ('terminal-distinct'/'unresolved'/'cap', ...)."""
    sides = []
    for s in (s1, s2):
        levels = []; cur = [s]; seen = Frontier()
        for d in range(depth + 1):
            keep = []
            for st in cur:
                g = graph_fn(st)
                if seen.add(g):
                    keep.append((st, g))
            levels.append(keep)
            if d == depth or not keep:
                break
            cur = [n for st, _ in keep for n in step_fn(st)]
            if len(cur) > cap:
                levels.append("cap"); break
        sides.append(levels)
    best = None
    L, R = sides
    for i, li in enumerate(L):
        if li == "cap": continue
        for j, rj in enumerate(R):
            if rj == "cap": continue
            fr = collections.defaultdict(list)
            for _, g in li:
                fr[inv(g)].append(g)
            if any(any(A.iso(g, x) for x in fr.get(inv(g), ())) for _, g in rj):
                if best is None or i + j < best[1] + best[2]:
                    best = ("join", i, j)
    if best:
        return best
    if any(l == "cap" for l in L + R):
        return ("cap", len(L) - 1, len(R) - 1)
    termL = len(L) < depth + 1 and not L[-1]; termR = len(R) < depth + 1 and not R[-1]
    return ("terminal-distinct" if (termL and termR) else "unresolved", len(L) - 1, len(R) - 1)


# typed witness engine: states are tuples of (type, vertex); rules map type->type
def witness_rules(m, n, neg=False):
    rules = [("R", "P1"), ("R", "Q1")]
    for k in range(1, m): rules.append((f"P{k}", f"P{k+1}"))
    for k in range(1, n): rules.append((f"Q{k}", f"Q{k+1}"))
    rules.append((f"P{m}", "Z")); rules.append((f"Q{n}", "Zp" if neg else "Z"))
    return rules


def w_step(rules):
    def f(st):
        out = []
        rels, ev = st
        for i, (t, v) in enumerate(rels):
            for a, b in rules:
                if a == t:
                    out.append((rels[:i] + ((b, v),) + rels[i + 1:], ev + 1))
        return out
    return f


def w_graph(st):
    return hyper_graph(st[0], typed=True)


def w_history_join(rules, m_level):
    """36A-style equal-length history comparison for the witness: events form chains; compare chain length + live types."""
    s = ((("R", 0),), 0)
    step = w_step(rules)
    first = step(s)
    L, R = [first[0]], [first[1]]
    for level in range(0, D + 1):
        if any(w_graph(a) is not None and A.iso(w_graph(a), w_graph(b)) and a[1] == b[1] for a in L for b in R):
            return ("join", level)
        L = [n for x in L for n in step(x)]; R = [n for x in R for n in step(x)]
        if not L or not R:
            return ("terminal-distinct", level)
    return ("unresolved", D)


def step3():
    out = dict(witnesses=[], A1=[], A2=None, A3=[])
    cases = [(m, m, False) for m in range(1, 7)] + [(2, 4, False), (3, 5, False), (3, 3, True)]
    for m, n, neg in cases:
        rules = witness_rules(m, n, neg)
        s = ((("R", 0),), 0)
        b1, b2 = w_step(rules)(s)
        F = frontier_join(w_step(rules), w_graph, b1, b2, depth=D)
        H = w_history_join(rules, m)
        out["witnesses"].append(dict(m=m, n=n, negative=neg, F=F, history_equal_length=H))
    # A1
    pos2 = A.Rule("xy,zy->xw,yw,zw", (("x", "y"), ("z", "y")), (("x", "w"), ("y", "w"), ("z", "w")), 2)
    h = A.make([(0, 3), (1, 3), (2, 3), (4, 0)])
    def st_fn(rule):
        return lambda hh: [A.apply(rule, hh, m) for m in A.matches(rule, hh)]
    def g_fn(hh):
        return hyper_graph([r for r, _ in hh[0]])
    for a, b in itertools.combinations(A.matches(pos2, h), 2):
        out["A1"].append(frontier_join(st_fn(pos2), g_fn, A.apply(pos2, h, a), A.apply(pos2, h, b)))
    # A2
    rule = A.R55744("L")
    spec_B = sorted(decode(75924531))   # phase 1 (B) mask of 55744
    hB = A.make(spec_B)
    ms = A.matches(rule, hB)
    out["A2"] = dict(matches=len(ms), F=frontier_join(st_fn(rule), g_fn, A.apply(rule, hB, ms[0]), A.apply(rule, hB, ms[1])))
    # A3
    level = [A.make([rule.norm(r) for r in A.START])]; seen = set()
    for d in range(5):
        nxt = []
        for hh in level:
            ms = A.matches(rule, hh)
            key = (tuple(sorted(r for r, _ in hh[0])), hh[2])
            if key not in seen:
                seen.add(key)
                for a, b in itertools.combinations(ms, 2):
                    if A.shared(a, b):
                        t0 = time.process_time()
                        res = frontier_join(st_fn(rule), g_fn, A.apply(rule, hh, a), A.apply(rule, hh, b))
                        out["A3"].append(dict(depth=d, F=res, cpu=round(time.process_time() - t0, 1)))
                        print("A3", d, res, flush=True)
            nxt += [A.apply(rule, hh, m) for m in ms]
        uniq = {}
        for hh in nxt:
            uniq.setdefault((tuple(sorted(r for r, _ in hh[0])), hh[2]), hh)
        level = list(uniq.values())
    return out


# ---------------------------------------------------------------- step 4 (conditional)
def rule_from_template(tpl, version, idx):
    role = {1: "y", 2: "z", 3: "u", 4: "v", 5: "w"}
    rhs = tuple(tuple(role[i] for i in t) for t in tpl)
    lhs = (("x", "y", "z"), ("y", "z", "x"), ("x", "u", "v"), ("v", "y", "z"))
    return A.Rule(f"{idx}{version}", lhs, rhs, 3, read=(3,), injective=True, setsem=True,
                  order=[("u", "v")] if version == "L" else [])


def step4(passing=(66427,)):
    data = {e["template_index"]: e for e in json.load(open(ROOT / "external" / "crr" / "W21_CLOSED_PARTICLE_HITS.json"))}
    out = []
    for idx in passing:
        e = data[idx]; rule = rule_from_template(e["template"], "C", idx)
        st_fn = lambda hh: [A.apply(rule, hh, m) for m in A.matches(rule, hh)]
        g_fn = lambda hh: hyper_graph([r for r, _ in hh[0]])
        for i, m in enumerate(e["state_masks"]):
            h = A.make(sorted(decode(m)))
            ms = A.matches(rule, h)
            for a, b in itertools.combinations(ms, 2):
                t0 = time.process_time()
                res = frontier_join(st_fn, g_fn, A.apply(rule, h, a), A.apply(rule, h, b))
                out.append(dict(rule=idx, phase=i, moves=len(ms), shared=A.shared(a, b), F=res,
                                cpu=round(time.process_time() - t0, 1)))
                print("step4", idx, i, res, flush=True)
    return out


def main(which="123"):
    out = ROOT / "results" / "EX-GC36B"; out.mkdir(parents=True, exist_ok=True)
    path = out / "results.json"
    res = json.loads(path.read_text()) if path.exists() else {}
    res["platform"] = dict(python=platform.python_version(), networkx=nx.__version__, machine=platform.machine(),
                           where="cloud container")
    for s, fn in (("1", step1), ("2", step2), ("3", step3), ("4", step4)):
        if s in which:
            t0 = time.process_time()
            res["step" + s] = fn(); res["cpu_step" + s] = time.process_time() - t0
            path.write_text(json.dumps(res, indent=1, default=list) + "\n")
            print("step", s, "cpu", res["cpu_step" + s], flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "123")
