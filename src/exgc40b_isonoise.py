"""EX-GC40B descriptive add-on (after P6 was scored): is the deterministic sector spread a fixed pattern or noise?
Repeats part_isotropy on a second, disjoint set of 20 graph seeds and correlates the two sector patterns."""
import json, numpy as np
import exgc40b as B
from scipy.sparse.csgraph import dijkstra
SKIP = []
def run(seeds):
    V = []
    for s in seeds:
        X, pairs = B.graph(s); n = np.ones(len(X))
        A = B.build(X, pairs, n, True, False, 0)
        src = B.nearest(X, (0, 0))[0]; T = dijkstra(A, indices=src)
        d = np.hypot(*(X - X[src]).T); m = (d > 95) & (d < 105) & np.isfinite(T)
        ang = np.arctan2(*(X[m] - X[src])[:, ::-1].T); k = ((ang + np.pi) / (2 * np.pi) * 16).astype(int) % 16
        if m.sum() < 1000:          # source site not in the giant component: skip this graph (recorded)
            SKIP.append(int(s)); continue
        v = d[m] / T[m]; V.append([v[k == q].mean() for q in range(16)])
    return np.array(V)
A, Bv = run(range(20)), run(range(20, 40))
a, b = A.mean(0), Bv.mean(0)
se = A.std(0, ddof=1).mean() / np.sqrt(len(A)) / a.mean()
out = dict(spread_A=float((a.max() - a.min()) / a.mean()), spread_B=float((b.max() - b.min()) / b.mean()),
           pattern_correlation=float(np.corrcoef(a, b)[0, 1]), sector_se_rel=float(se),
           expected_spread_from_noise_approx=float(3.5 * se), skipped_seeds=SKIP, n_graphs=[len(A), len(Bv)])
json.dump(out, open('../results/EX-GC40B/isotropy_noise_check.json', 'w'), indent=1); print(out)
