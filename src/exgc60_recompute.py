"""EX-GC60 recomputation -> results/EX-GC60/recompute.json. sigma = 3, xy: rotation and shear rates rerun at dt = 0.05 with the spin read
from X.sdens sums and slopes by np.polyfit (separate readout and fit); the lattice curl of the shear field (zero in the interior, the
periodic-roll wrap at the box faces is where PS and P can differ). Shared: probe builder, field builder, Hamiltonian.
EXPLORATORY. Our construction; not LaViolette's Model G and not Chain Theory."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import exgc51 as G, exgc53 as T53, exgc58 as X, exgc60 as P


def run(Bs, dt=0.05):
    one = np.ones((P.NB,) * 3); psi = P.probe(3.0, "x"); C = X.curl_op(Bs) if Bs is not None else None; ph = []; t = 0.0
    H = lambda p: G.Hm(p, one, P.M, 0.0) + ((T53.Hb(p, Bs) + X.Hpl(p, C)) if Bs is not None else 0)
    def rec(p):
        s = X.sdens(p); ph.append(np.arctan2(s[1].sum(), s[0].sum()))
    rec(psi)
    while t < P.T - 1e-9:
        k1 = -1j * H(psi); k2 = -1j * H(psi + dt / 2 * k1); k3 = -1j * H(psi + dt / 2 * k2); k4 = -1j * H(psi + dt * k3); psi = psi + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4); t += dt; rec(psi)
    return np.array(ph), dt


if __name__ == "__main__":
    z, dt = run(None); res = {}
    for kind in ("rot", "shear"):
        p, _ = run(P.field(kind, "xy")); res[kind] = float(np.polyfit(dt * np.arange(len(p)), p - z, 1)[0]) / P.S0
    A3 = json.loads((P.OUT / "EX-GC60_A3.json").read_text())["rows"]; res["rel_rot"] = res["rot"] / A3["xy_rot_PS"]["S"] - 1; res["rel_shear"] = res["shear"] / A3["xy_shear_PS"]["S"] - 1
    Cs = X.curl_op(P.field("shear", "xy")); inner = (slice(None),) + (slice(2, -2),) * 3
    res["curl_shear_interior_max"] = float(np.abs(Cs[inner]).max()); res["curl_shear_face_max"] = float(np.abs(Cs).max())
    print(res); (P.OUT / "recompute.json").write_text(json.dumps(res, indent=1, default=float))
