"""EX-GC34W: cubic Einstein–Hilbert vertex-average action on Y1."""
import json
import platform
import time
from itertools import permutations
from pathlib import Path

import numpy as np

import exgc34n as N

ETA = N.ETA
OFFSETS = np.zeros((4, 4, 3), int)
for a in range(4):
    for b in range(4):
        for index in (a, b):
            if index:
                OFFSETS[a, b, index-1] ^= 1


def average_factor(parity, q):
    return np.prod(np.where(parity, np.cos(q/2), 1), axis=-1)


def leg(metric, k, omega, h):
    q = np.zeros(3) if h is None else h*k
    K = k if h is None else 2*np.sin(q/2)/h
    p = np.r_[-omega, K]
    derivative = np.empty((4, 4, 4), complex)
    for i in range(4):
        parity = OFFSETS.copy()
        if i:
            parity[:, :, i-1] ^= 1
        derivative[i] = 1j*p[i] * average_factor(parity, q) * metric
    h_center = average_factor(OFFSETS, q)*metric
    delta_inverse = -ETA @ h_center @ ETA
    density_inverse = delta_inverse + 0.5*np.trace(ETA @ h_center)*ETA
    lower = (derivative.transpose(1, 0, 2) + derivative.transpose(1, 2, 0) - derivative) / 2
    gamma = np.einsum("rs,smn->rmn", ETA, lower)
    return density_inverse, delta_inverse, lower, gamma


def contract(density, first, second):
    first_product = np.einsum("mn,amb,bna->", density, first, second)
    second_product = np.einsum("mn,amn,bab->", density, first, second)
    return first_product - second_product


def vertex(a, b, c):
    total = 0j
    for (dens, inv, _, _), (_, _, low, g1), (_, _, _, g2) in permutations((a, b, c)):
        gamma2 = np.einsum("rs,smn->rmn", inv, low)
        total += contract(dens, g1, g2) + contract(ETA, gamma2, g2) + contract(ETA, g2, gamma2)
    return total


def direct_gauge_leg(a, b, c):
    # Independent coefficient extraction: evaluate the cubic polynomial of sums.
    def action(terms):
        dens, inv, low, g = [sum(term[j] for term in terms) for j in range(4)]
        g2 = np.einsum("rs,smn->rmn", inv, low)
        return contract(dens, g, g) + contract(ETA, g2, g) + contract(ETA, g, g2)

    return action((a,b,c)) - action((a,b)) - action((a,c)) - action((b,c)) + action((a,)) + action((b,)) + action((c,))


def dispersion(k, h):
    return k if h is None else 2*np.sin(h*k/2)/h


def wave(k, h, pol):
    K = dispersion(k, h)
    omega = np.linalg.norm(K)
    axis = K / omega
    reference = np.eye(3)[np.argmin(abs(axis))]
    u = np.cross(axis, reference)
    u /= np.linalg.norm(u)
    v = np.cross(axis, u)
    metric = np.zeros((4,4))
    metric[1:,1:] = ((np.outer(u,u)-np.outer(v,v)) if pol == 0
                     else (np.outer(u,v)+np.outer(v,u))) / np.sqrt(2)
    G = N.G_up(metric, np.r_[-omega, K])
    return leg(metric, k, omega, h), omega, float(np.linalg.norm(G)/omega**2)


def gauge_metric(p, index):
    xi = np.eye(4)[index]
    return 1j * (np.outer(p, xi) + np.outer(xi, p))


def run(pair, pol, h):
    k1, k2 = (np.array(d)*2*np.pi/8 for d in pair)
    a, w1, ea = wave(k1,h,pol[0])
    b, w2, eb = wave(k2,h,pol[1])
    k, w = -k1-k2, -w1-w2
    basis = [v.reshape(4,4) for v in N.PM.T]
    source = np.array([vertex(a,b,leg(m,k,w,h)) for m in basis])
    contractions, separate = [], []
    p = np.r_[-w, dispersion(k,h)]
    for gauge in range(4):
        m = gauge_metric(p, gauge)
        coords = np.r_[-m[0,0]/2,m[0,1:],[m[i+1,j+1] for i,j in N.PAIRS]]
        scale = np.linalg.norm(source)*np.linalg.norm(coords)
        value = source @ coords
        direct = direct_gauge_leg(a,b,leg(m,k,w,h))
        contractions.append(float(abs(value)/scale))
        separate.append(float(abs(value-direct)/scale))
    return dict(h=h, pair=pair, polarizations=pol, input_residual=max(ea,eb),
                source_norm=float(np.linalg.norm(source)), gauge=contractions,
                separate_agreement=separate)


def main():
    start = time.process_time()
    x, diagonal, yz = (1,0,0),(1,1,1),(0,1,1)
    channels = [((d,d),(p,p)) for d in (x,diagonal) for p in range(2)]
    channels += [((x,yz),(p,q)) for p in range(2) for q in range(2)]
    continuum = [run(pair, p, None) for pair,p in channels]
    lattice = [run(pair, p, h) for h in (.5,.25,.125) for pair,p in channels]
    w0 = (max(x["input_residual"] for x in continuum+lattice) <= 1e-12
          and max(max(x["gauge"]) for x in continuum) <= 1e-10)
    w1 = max(max(x["gauge"]) for x in lattice) <= 1e-10
    result = dict(W0=w0,W1=w1, continuum=continuum,lattice=lattice,
                  cpu_seconds=time.process_time()-start, python=platform.python_version(), numpy=np.__version__)
    path = Path(__file__).resolve().parents[1] / "results" / "EX-GC34W"
    path.mkdir(parents=True,exist_ok=True)
    (path / "results.json").write_text(json.dumps(result,indent=2)+"\n")
    print("W0",w0,"W1",w1)
    for group in (continuum,*[lattice[i:i+8] for i in (0,8,16)]):
        print(group[0]["h"], max(max(x["gauge"]) for x in group),
              max(max(x["separate_agreement"]) for x in group))


if __name__ == "__main__":
    main()
