"""Regression checks for the Model Gc solvers against values recorded in ModelG (EX008, EX009, MG195).

python tests/test_regression.py          radial checks (a few seconds)
python tests/test_regression.py --pair   also the axisymmetric pair regression (EX009 k = 0.5, d = 20; about 1 minute)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
import ex008_longrange as e  # noqa: E402


def check(name, got, want, rtol):
    ok = abs(got - want) <= rtol * abs(want)
    print(f"{'PASS' if ok else 'FAIL'}  {name}: got {got:.7g}, recorded {want:.7g} (rtol {rtol:g})")
    return ok


def main():
    ok = True
    # MG195 P1 radial point (G feeds X, k = 0.2, R = 120, dr = 0.05): continue from the recorded profile.
    prof = np.load(os.path.join(e.ROOT, "results", "MG195_inverse_square", "radial_P1_R120.npz"))
    M, u, res, conv = e.particle(120.0, 0.05, 0.2, 1, 0, prev=(prof["r"], prof["u"]))
    ok &= conv and res < 1e-8
    ok &= check("MG195 P1 consumption M", M.integral(0.2 * u[3] * u[2]), 5.317804244620182, 1e-6)
    ok &= check("MG195 P1 mobility mu", M.drift(u)["mu"], -0.02553789034253153, 1e-6)
    # Far field: G = G_inf - M (1/r - 1/R) / (4 pi) at r = 40.
    Mc = M.integral(0.2 * u[3] * u[2])
    i = np.argmin(abs(M.r - 40.0))
    ok &= check("1/r far field at r = 40", e.GINF - u[3][i], Mc / (4 * np.pi) * (1 / M.r[i] - 1 / 120.0), 1e-3)
    if "--pair" in sys.argv:
        import ex009_pair2d as p2
        p2.K, p2.CX, p2.CY = 0.5, 1, 0
        pr = np.load(os.path.join(e.ROOT, "results", "EX008_longrange", "particle_GX_k0.5_R120.npz"))
        Mrad = e.Model(120.0, 0.05, 0.5, 1, 0).integral(0.5 * pr["u"][3] * pr["u"][2])
        a = p2.Axi(20.0)
        _, v, res, conv, _ = a.solve(a.guess(pr["r"], pr["u"], Mrad))
        ok &= conv and check("EX009 pair speed k = 0.5, d = 20", v, -9.287891e-05, 1e-5)
    print("ALL PASS" if ok else "SOME FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
