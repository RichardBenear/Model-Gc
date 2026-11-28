"""EX-GC63: 54's directed signal-count rule counted per local tick (D2-E1): source 4 kappa * (1/2){H1, v1} with H1, v1 = 54's
Hop/vop at N = 1 (flat weights). 63A co-moving pair (53E/54B), 63B resting polarised source magnetic moment vs the PS rule, C1-C3.
Plan: spec/EX-GC63_plan.md. EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys, time, platform
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc54 as E54, exgc53 as T53, exgc59 as E59, exgc58 as X

ROOT = Path(__file__).resolve().parents[1]; OUT = ROOT / "results" / "EX-GC63"
_orig_source = E54.source


def d2e1(psi, m):
    return E54.rule_current(psi, np.ones(psi.shape[:3]), m, "D2E")


def source(ps, N, ms, rule):
    if rule == "D2E1":
        return sum(d2e1(p, mm) for p, mm in zip(ps, ms))
    return _orig_source(ps, N, ms, rule)


E54.source = source          # evolve() and run_pair() look the name up in exgc54 at call time


def moment(j, n):
    LP = T53.link_pos(n); m = np.zeros(3)
    for a in range(3):
        e = np.eye(3)[a]
        for i in range(3):
            m[i] += 0.5 * np.sum(np.cross(np.stack([LP[a][k] for k in range(3)], -1), e)[..., i] * j[a])
    return m


def part_B(n=97):
    res = {}; one = np.ones((n, n, n))
    for nm, chi in (("up", E59.UP), ("down", E59.DOWN)):
        p = E59.source(n, chi); P = T53.link_P(p); PS = P + 0.25 * X.curl_T(X.sdens(p)); s = d2e1(p, E59.M)
        res[nm] = dict(m_D2E1=moment(s, n).tolist(), m_PS=moment(PS, n).tolist(), m_P=moment(P, n).tolist(),
                       ratio=float(np.linalg.norm(moment(s, n)) / np.linalg.norm(moment(PS, n))),
                       tot_D2E1=s.sum(axis=(1, 2, 3)).tolist(), maxabs_D2E1_over_PS=float(np.abs(s).max() / np.abs(PS).max()))
    return res


def part_checks(n=65):
    m = E54.S.M; p1, p2, kap = E54.packets(0.2, n); q1, _, _ = E54.packets(-0.2, n); N1 = np.ones((n,) * 3)
    s1 = d2e1(p1, m); sm = d2e1(q1, m); P1 = T53.link_P(p1)
    c1 = float(np.abs(d2e1(np.zeros_like(p1), m)).max())
    c2 = float(np.abs(s1.sum(axis=(1, 2, 3)) + sm.sum(axis=(1, 2, 3))).max() / np.abs(s1.sum(axis=(1, 2, 3))).max())
    c3 = float(np.abs(s1 - E54.rule_current(p1, N1, m, "D2E")).max() / np.abs(s1).max())
    f, g, info = E54.S.initial_field([p1, p2], [m, m], 0.0, kap, n, n // 2, 0.2); N = E54.Nf(f)
    sp = source([p1, p2], N, [m, m], "D2E1"); Pp = sum(T53.link_P(p) for p in (p1, p2)); sE = source([p1, p2], N, [m, m], "D2E")
    return dict(C1_zero=c1, C2_flip=c2, C3_flatN=c3, P1_strength_pair=float(4 * sp[0].sum() / Pp[0].sum()),
                D2E_inplace_strength_pair=float(4 * sE[0].sum() / Pp[0].sum()), single_strength_N1=float(4 * s1[0].sum() / P1[0].sum()))


def save(name, res, t0):
    res["_meta"] = dict(wall_s=time.time() - t0, platform=platform.platform(), python=platform.python_version(), numpy=np.__version__)
    (OUT / f"EX-GC63_{name}.json").write_text(json.dumps(res, indent=1))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True); job = sys.argv[1]; t0 = time.time()
    if job == "B":
        save("B", part_B(), t0)
    elif job == "checks":
        save("checks", part_checks(), t0)
    else:                        # A_0.2_pair, A_0.0_one, ...
        _, V, lab = job.split("_"); save(job, E54.run_pair(float(V), "D2E1", only=(lab,)), t0)
