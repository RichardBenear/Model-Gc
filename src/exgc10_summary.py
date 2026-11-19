"""EX-GC10 summary: twist energy, Z2 index and G consumption of the tethered particle (reads records only)."""
import glob
import json
import os

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "EX-GC10")
FIG = os.path.join(ROOT, "paper", "figures", "EX-GC10_tether.png")


def load(pat):
    return {os.path.basename(f)[:-5]: json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, pat)))}


def late(rows, key, t0):
    v = [r[key] for r in rows if r["t"] >= t0]
    return float(np.mean(v)), float(np.std(v))


def main():
    summ = dict(frames={}, det={}, tok={})
    for k, d in load("frames_*.json").items():
        h = d["hist"]
        summ["frames"][k] = dict(E0=h[0][1], E_final=h[-1][1], E_min=min(x[1] for x in h), z2_0=h[0][2],
                                 z2_final=h[-1][2])
    for kind in ("det", "tok"):
        for k, d in load(f"{kind}_*.json").items():
            r = d["rows"]
            tl = 0.7 * r[-1]["t"]
            Mm, Ms = late(r, "M", tl)
            Em, _ = late(r, "E", tl)
            rec = dict(turns=d["turns"], eps=d["eps"], rule=d["rule"], seed=d["seed"], T=r[-1]["t"], E0=r[0]["E"],
                       E_final=r[-1]["E"], E_late=Em, z2_0=r[0]["z2"], z2_final=r[-1]["z2"],
                       z2_constant=all(x["z2"] == r[0]["z2"] for x in r), M0=r[0]["M"], M_late=Mm, M_late_sd=Ms,
                       E_by_layer=d.get("E_by_layer"))
            if "Wcen" in r[-1]:
                rec["Wcen_max_dist"] = float(max(np.linalg.norm(x["Wcen"]) for x in r))
            summ[kind][k] = rec
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    for kind in ("det", "tok"):
        print(kind)
        for k, v in summ[kind].items():
            print(f"  {k:40s} E {v['E0']:6d} -> {v['E_final']:5d}  z2 {v['z2_0']:+.0f}->{v['z2_final']:+.0f} "
                  f"const {v['z2_constant']}  M_late {v['M_late']:.3f} ({v['M_late_sd']:.3f})"
                  + (f"  drift {v['Wcen_max_dist']:.2f}" if "Wcen_max_dist" in v else ""))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    det = load("det_*.json")
    sty = {0: ("k", "0 turns"), 1: ("C3", "1 turn (360°)"), 2: ("C0", "2 turns (720°)")}
    for k, d in det.items():
        r = d["rows"]
        t = [x["t"] for x in r]
        if d["rule"] == "flag":
            col, lab, ls = "C2", "1 turn, pairwise rule (control)", ":"
        elif d["eps"] == 0:
            col, lab, ls = "C1", "1 turn, ε = 0 (control)", "--"
        else:
            col, lab = sty[d["turns"]]
            ls = "-"
        lab = lab if d["seed"] == 0 else None
        ax[0].plot(t, [max(x["E"], 1) for x in r], ls, color=col, label=lab, lw=1.2)
        ax[1].plot(t, [x["M"] for x in r], ls, color=col, label=lab, lw=1.2)
    ax[0].set_yscale("log")
    ax[0].set_xlabel("time")
    ax[0].set_ylabel("mismatched edges (twist energy)")
    ax[0].set_title("Tether twist, deterministic chemistry")
    ax[0].legend(fontsize=7)
    ax[1].set_xlabel("time")
    ax[1].set_ylabel("G consumption M")
    ax[1].set_title("Particle consumption (ε = 0.5: twisted edges halve G passage)")
    for k, d in load("tok_*.json").items():
        r = d["rows"]
        col, lab = sty[d["turns"]]
        ax[2].plot([x["t"] for x in r], [x["M"] for x in r], color=col, alpha=1.0 if d["seed"] == 0 else 0.5,
                   label=f"{lab}, seed {d['seed']}", lw=1)
    ax[2].set_xlabel("time")
    ax[2].set_ylabel("G consumption M")
    ax[2].set_title("Token chemistry, Ω = 3000")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG, dpi=130)


if __name__ == "__main__":
    main()
