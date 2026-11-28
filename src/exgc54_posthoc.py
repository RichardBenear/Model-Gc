"""EX-GC54 post-hoc diagnostic (not scored): is the smooth-gauge loop delay (54D) set by the shift at the clock site?
Same pattern with phase chosen so that beta at the launch/arrival site is maximal (sin) or near zero (cos)."""
import json
from pathlib import Path
import numpy as np
import exgc53 as T
L = 256; one = np.ones(L); s = np.arange(L); out = {}
for Lam in (32, 64, 128):
    for ph, lab in ((0.0, "sin"), (np.pi / 2, "cos")):
        chi = 1e-3 * np.sin(2 * np.pi * s / Lam + ph); bl = np.roll(chi, -1) - chi; d = T.loop_dt(one, bl)
        out[f"{Lam}_{lab}"] = dict(beta_at_clock=float(bl[0]), dt=float(d)); print(Lam, lab, out[f"{Lam}_{lab}"], flush=True)
(Path(__file__).resolve().parents[1] / "results/EX-GC54/posthoc_gauge_clock.json").write_text(json.dumps(out, indent=1))
