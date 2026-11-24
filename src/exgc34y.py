"""EX-GC34Y: majority-parity contractions and the full cubic Ward test."""
import json
import platform
import time
from itertools import permutations, product
from pathlib import Path

import numpy as np

import exgc34w as W

ETA = W.ETA
OFF = W.OFFSETS
TRIPLE = OFF[:, :, None] ^ OFF[0, None, None, :]
PERMS = tuple(permutations(range(3)))
IDX = range(4)


def monomials():
    groups = [[] for _ in range(3)]
    for m, n, a, b in product(IDX, repeat=4):
        for sign, first, second in (
            (1, (a,m,b), (b,n,a)),
            (-1, (a,m,n), (b,a,b)),
        ):
            groups[0].append((sign, (m,n), first, second))
            if m != n:
                continue
            for s in IDX:
                groups[1].append((sign*ETA[m,n], (first[0],s), (s,*first[1:]), second))
                groups[2].append((sign*ETA[m,n], (second[0],s), (s,*second[1:]), first))
    return groups


def assemble_group(values):
    signs, first, second, third = zip(*values)
    indices = tuple(tuple(np.array(x).T) for x in (first, second, third))
    parities = [lookup[idx] for lookup, idx in zip((OFF, TRIPLE, TRIPLE), indices)]
    majority = np.sum(parities, axis=0) >= 2
    return np.array(signs), indices, parities, majority


GROUPS = [assemble_group(v) for v in monomials()]


def derivative(k, frequency, h):
    return 1j*np.r_[-frequency, W.dispersion(k,h)]


def leg(metric, k, frequency, h):
    q = np.zeros(3) if h is None else k*h
    d = derivative(k,frequency,h)
    raw = d[:,None,None]*metric
    lower = (raw.transpose(1,0,2)+raw.transpose(1,2,0)-raw)/2
    inverse = -ETA @ metric @ ETA
    density = inverse + .5*np.trace(ETA @ metric)*ETA
    gamma = np.einsum("rs,smn->rmn",ETA,lower)
    return (density,inverse,lower,gamma), q


def coefficient(legs, majority):
    value = 0j
    for a,b,c in PERMS:
        (arrays_a, qa), (arrays_b, qb), (arrays_c, qc) = [legs[i] for i in (a,b,c)]
        for i, (signs, indices, parities, common) in enumerate(GROUPS):
            anchor = common if majority else False
            factors = (arrays_a[0 if i == 0 else 1], arrays_b[3 if i == 0 else 2], arrays_c[3])
            items = [f[idx] * W.average_factor(parity ^ anchor, q)
                     for f, idx, parity, q in zip(factors, indices, parities, (qa,qb,qc))]
            value += np.sum(signs * items[0]*items[1]*items[2])
    return value


def polarization(k, h, pol):
    k_eff = W.dispersion(k,h)
    omega = np.linalg.norm(k_eff)
    axis = k_eff/omega
    u = np.cross(axis,np.eye(3)[np.argmin(abs(axis))])
    u /= np.linalg.norm(u)
    v = np.cross(axis,u)
    out = np.zeros((4,4))
    out[1:,1:] = ((np.outer(u,u)-np.outer(v,v)) if pol == 0
                  else (np.outer(u,v)+np.outer(v,u)))/np.sqrt(2)
    return out,omega


def gauge_coords(m):
    return np.r_[-m[0,0]/2,m[0,1:],[m[a+1,b+1] for a,b in W.N.PAIRS]]


def on_shell(pair, pol, h, majority):
    k1,k2 = [np.array(d)*np.pi/4 for d in pair]
    m1,w1 = polarization(k1,h,pol[0])
    m2,w2 = polarization(k2,h,pol[1])
    a,b = leg(m1,k1,w1,h),leg(m2,k2,w2,h)
    k,w = -k1-k2,-w1-w2
    source = np.array([coefficient((a,b,leg(v.reshape(4,4),k,w,h)),majority) for v in W.N.PM.T])
    p = np.r_[-w,W.dispersion(k,h)]
    norms = []
    for i in IDX:
        z = gauge_coords(W.gauge_metric(p,i))
        norms.append(float(abs(source @ z)/(np.linalg.norm(source)*np.linalg.norm(z))))
    return norms


def nonlinear(xi, kg, wg, metric, k, w, h, majority):
    d = derivative(k,w,h)
    dg = derivative(kg,wg,h)
    q = np.zeros(3) if h is None else k*h
    qg = np.zeros(3) if h is None else kg*h
    xi_up = ETA @ xi
    result = np.zeros((4,4),complex)
    for m,n in product(IDX,repeat=2):
        target = OFF[m,n] if majority else False
        for c in IDX:
            s = OFF[0,c]
            result[m,n] += (xi_up[c]*d[c]*metric[m,n] *
                            W.average_factor(s ^ target,qg)*W.average_factor((s ^ OFF[m,n]) ^ target,q))
            for i,j in ((m,n),(n,m)):
                result[m,n] += (dg[i]*xi_up[c]*metric[c,j] *
                                W.average_factor(OFF[i,c] ^ target,qg)*W.average_factor(OFF[c,j] ^ target,q))
    return result


def full_terms(pair, frequencies, metrics, gauge, h, majority):
    k1,k2 = pair
    w1,w2 = frequencies
    a,b = metrics
    kg,wg = -k1-k2,-w1-w2
    xi = np.eye(4)[gauge]
    gm = W.gauge_metric(np.r_[-wg,W.dispersion(kg,h)],gauge)
    v = coefficient((leg(a,k1,w1,h),leg(b,k2,w2,h),leg(gm,kg,wg,h)),majority)
    ga = W.N.G_up(a,np.r_[-w1,W.dispersion(k1,h)])
    gb = W.N.G_up(b,np.r_[-w2,W.dispersion(k2,h)])
    return np.array([v, -np.sum(gb*nonlinear(xi,kg,wg,a,k1,w1,h,majority)),
                     -np.sum(ga*nonlinear(xi,kg,wg,b,k2,w2,h,majority))])


def main():
    start = time.process_time()
    x,d,yz = (1,0,0),(1,1,1),(0,1,1)
    channels = [((k,k),(p,p)) for k in (x,d) for p in range(2)]
    channels += [((x,yz),(p,q)) for p,q in product(range(2),repeat=2)]
    continuum_on = [on_shell(k,p,None,True) for k,p in channels]
    control = max(np.max(abs(np.asarray(on_shell(k,p,.25,False))-W.run(k,p,.25)["gauge"]))
                  for k,p in channels)
    rng = np.random.default_rng(343425)
    triads = []
    for _ in range(16):
        k = rng.integers(-2,3,(2,3))*np.pi/4
        w = rng.normal(size=2)
        m = rng.normal(size=(2,4,4))
        triads.append((k,w,(m+m.transpose(0,2,1))/2))
    reference = [[full_terms(*case,i,None,True) for i in IDX] for case in triads]
    continuum_off = [[float(abs(np.sum(v))/np.sum(abs(v))) for v in row] for row in reference]
    rows = []
    for h in (.5,.25,.125):
        old_on = [on_shell(k,p,h,False) for k,p in channels]
        new_on = [on_shell(k,p,h,True) for k,p in channels]
        old_off,new_off = [],[]
        for case, norms in zip(triads,reference):
            scales = [np.sum(abs(terms)) for terms in norms]
            old_off.append([float(abs(np.sum(full_terms(*case,i,h,False)))/scales[i]) for i in IDX])
            new_off.append([float(abs(np.sum(full_terms(*case,i,h,True)))/scales[i]) for i in IDX])
        rows.append(dict(h=h, vertex_on=old_on, majority_on=new_on,
                         vertex_off=old_off, majority_off=new_off))
    maximum = lambda group: float(np.max(group))
    gates = dict(Y0=bool(control <= 1e-10 and max(maximum(continuum_on),maximum(continuum_off)) <= 1e-10),
                 Y1=all(maximum(row["majority_on"]) <= 1e-10 for row in rows),
                 Y2=all(maximum(row["majority_off"]) <= 1e-10 for row in rows))
    result = dict(gates=gates, vertex_reconstruction=float(control), continuum_on=continuum_on,
                  continuum_off=continuum_off, rows=rows, cpu_seconds=time.process_time()-start,
                  python=platform.python_version(), numpy=np.__version__)
    path = Path(__file__).resolve().parents[1]/"results"/"EX-GC34Y"
    path.mkdir(parents=True,exist_ok=True)
    (path/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(gates,control,maximum(continuum_on),maximum(continuum_off),flush=True)
    for row in rows:
        print(row["h"], *(maximum(row[k]) for k in ("vertex_on","majority_on","vertex_off","majority_off")))


if __name__ == "__main__":
    main()
