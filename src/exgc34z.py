"""EX-GC34Z: joint short-range cubic action and gauge identity fit."""
import os

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[key] = "1"

import json
import platform
import time
from itertools import product
from pathlib import Path

import numpy as np

import exgc34x as X
import exgc34y as Y

W = Y.W
OFFSETS = X.OFFSETS
ETA = W.ETA


def translations(q1, q2):
    return np.cos(OFFSETS[:, 0] @ q1 + OFFSETS[:, 1] @ q2) - 1


def action(legs):
    base = 0j
    features = np.zeros((3, len(OFFSETS)), complex)
    for a, b, c in Y.PERMS:
        (arr_a, qa), (arr_b, qb), (arr_c, qc) = (legs[i] for i in (a, b, c))
        shift = translations(qa, qb)
        for i, (signs, indices, parities, anchor) in enumerate(Y.GROUPS):
            vals = (arr_a[0 if i == 0 else 1], arr_b[3 if i == 0 else 2], arr_c[3])
            factors = [val[idx] * W.average_factor(parity ^ anchor, q)
                       for val, idx, parity, q in zip(vals, indices, parities, (qa, qb, qc))]
            coefficient = np.sum(signs * factors[0] * factors[1] * factors[2])
            base += coefficient
            features[i] += coefficient * shift
    return base, features.ravel()


def gauge_parts(xi, kg, wg, metric, k, w, h):
    d, dg = Y.derivative(k,w,h), Y.derivative(kg,wg,h)
    q, qg = k*(h or 0), kg*(h or 0)
    xi_up = ETA @ xi
    transport = np.zeros((4,4), complex)
    stretch = transport.copy()
    for m,n in product(Y.IDX, repeat=2):
        target = Y.OFF[m,n]
        for c in Y.IDX:
            s = Y.OFF[0,c]
            transport[m,n] += (xi_up[c]*d[c]*metric[m,n] *
                               W.average_factor(s ^ target,qg)*W.average_factor(s,q))
            for i,j in ((m,n),(n,m)):
                stretch[m,n] += (dg[i]*xi_up[c]*metric[c,j] *
                                 W.average_factor(Y.OFF[i,c] ^ target,qg) *
                                 W.average_factor(Y.OFF[c,j] ^ target,q))
    return transport, stretch


def off_shell(case, gauge, h):
    (k1,k2), (w1,w2), (a,b) = case
    kg, wg = -k1-k2, -w1-w2
    gm = W.gauge_metric(np.r_[-wg,W.dispersion(kg,h)],gauge)
    base, action_features = action((Y.leg(a,k1,w1,h), Y.leg(b,k2,w2,h), Y.leg(gm,kg,wg,h)))
    reference = Y.full_terms(*case,gauge,h,True)
    gauge_features = np.zeros((2,len(OFFSETS)),complex)
    transforms = []
    for k,w,metric,other_k,other_w,other in ((k1,w1,a,k2,w2,b),(k2,w2,b,k1,w1,a)):
        G = W.N.G_up(other,np.r_[-other_w,W.dispersion(other_k,h)])
        terms = np.array([-np.sum(G*part) for part in gauge_parts(np.eye(4)[gauge],kg,wg,metric,k,w,h)])
        transforms.append(terms.sum())
        gauge_features += terms[:,None]*translations(kg*(h or 0),k*(h or 0))
    scale = np.sum(abs(Y.full_terms(*case,gauge,None,True)))
    discrepancy = np.max(abs(np.r_[base,transforms]-reference))/scale
    return ((base+sum(transforms))/scale,
            np.r_[action_features,gauge_features.ravel()]/scale,float(discrepancy))


def on_shell(pair, polarizations, h):
    k1,k2 = [np.asarray(d)*np.pi/4 for d in pair]
    a,w1 = Y.polarization(k1,h,polarizations[0])
    b,w2 = Y.polarization(k2,h,polarizations[1])
    k,w = -k1-k2,-w1-w2
    legs = Y.leg(a,k1,w1,h), Y.leg(b,k2,w2,h)
    varied = [action((*legs,Y.leg(m.reshape(4,4),k,w,h))) for m in W.N.PM.T]
    base, features = (np.array(z) for z in zip(*varied))
    baseline, rows = [], []
    for i in Y.IDX:
        z = Y.gauge_coords(W.gauge_metric(np.r_[-w,W.dispersion(k,h)],i))
        scale = np.linalg.norm(base)*np.linalg.norm(z)
        baseline.append(base @ z / scale)
        rows.append(np.r_[z @ features, np.zeros(48)]/scale)
    reference = Y.on_shell(pair,polarizations,h,True)
    return np.array(baseline),np.array(rows),float(np.max(abs(abs(np.array(baseline))-reference)))


def triads():
    rng = np.random.default_rng(343426)
    result = []
    for _ in range(36):
        k = rng.integers(-2,3,(2,3))*np.pi/4
        w = rng.normal(size=2)
        m = rng.normal(size=(2,4,4))
        result.append((k,w,(m+m.transpose(0,2,1))/2))
    return result[:24],result[24:]


def flatten(cases,h):
    evaluated = [off_shell(case,i,h) for case in cases for i in Y.IDX]
    b,rows,checks = (np.array(x) for x in zip(*evaluated))
    return b,rows,float(np.max(checks))


def fit(M,b):
    A, target = np.vstack((M.real,M.imag)), -np.r_[b.real,b.imag]
    U,s,V = np.linalg.svd(A,full_matrices=False)
    rows = []
    for cutoff in (1e-10,1e-12,1e-13):
        rank = int(np.sum(s>cutoff*s[0]))
        uk = U[:,:rank]
        x = V[:rank].T @ ((uk.T @ target)/s[:rank])
        projection = uk @ (uk.T @ target)-target
        rows.append(dict(cutoff=cutoff,rank=rank,coeff_norm=float(np.linalg.norm(x)),
                         max_residual=float(np.max(abs(M@x+b))),
                         max_lower_bound=float(np.linalg.norm(projection)/np.sqrt(len(b))),
                         projection_agreement=float(np.linalg.norm(A@x-target-projection))))
        if cutoff == 1e-12:
            chosen = x
    return rows,chosen,s.tolist()


def main():
    start = time.process_time()
    training,unseen = triads()
    evaluated = [on_shell(pair,p,.25) for pair,p in X.BASE_CHANNELS]
    ob = np.concatenate([r[0] for r in evaluated])
    oM = np.vstack([r[1] for r in evaluated])
    b,M,check = flatten(training,.25)
    controls = max(check,max(r[2] for r in evaluated))
    cb,_,_ = flatten(training+unseen,None)
    continuum = max(float(max(abs(cb))),max(float(np.max(abs(on_shell(pair,p,None)[0]))) for pair,p in X.BASE_CHANNELS))
    joint,coeff,spectrum = fit(np.vstack((oM,M)),np.r_[ob,b])
    necessary,on_coeff,on_singular = fit(oM[:,:72],ob)
    holdout = []
    for h in (.5,.25,.125,.0625):
        hb,hM,_ = flatten(unseen,h)
        holdout.append(dict(h=h, baseline_max=float(max(abs(hb))),fitted_max=float(max(abs(hb+hM@coeff))),
                            baseline=abs(hb).tolist(), fitted=abs(hb+hM@coeff).tolist()))
    verdict = joint[1]["max_residual"] <= 1e-10
    result = dict(Z0=controls <= 1e-10 and continuum <= 1e-10,Z1=verdict,
                  Z2_scored=verdict,Z2=max(row["fitted_max"] for row in holdout if row["h"] in (.5,.125)) <= 1e-10,
                  baseline_reconstruction=controls,continuum=continuum,
                  joint_baseline=float(max(abs(np.r_[ob,b]))), joint_fit=joint,on_shell_fit=necessary,
                  joint_coefficients=coeff.tolist(),on_shell_coefficients=on_coeff.tolist(),
                  joint_singular_values=spectrum,on_shell_singular_values=on_singular,
                  fitted_on_shell=abs(ob+oM@coeff).tolist(),fitted_off_shell=abs(b+M@coeff).tolist(),
                  holdout=holdout,cpu_seconds=time.process_time()-start,
                  python=platform.python_version(),numpy=np.__version__)
    out = Path(__file__).resolve().parents[1]/"results"/"EX-GC34Z"
    out.mkdir(parents=True,exist_ok=True)
    (out/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result["Z0"],result["Z1"],controls,continuum,joint,necessary,
          [(r["h"],r["fitted_max"]) for r in holdout],flush=True)



def posthoc():
    """Add previously held-out equations without rescoring the registered fit."""
    training, unseen = triads()
    waves = [on_shell(pair,p,.25) for pair,p in X.BASE_CHANNELS]
    b0 = np.concatenate([r[0] for r in waves])
    m0 = np.vstack([r[1] for r in waves])
    checks = {}
    for name, spacings in (("all_at_training_spacing",(.25,)),("all_at_three_spacings",(.5,.25,.125))):
        blocks = [flatten(training+unseen,h) for h in spacings]
        b = np.r_[b0,*[r[0] for r in blocks]]
        M = np.vstack([m0,*[r[1] for r in blocks]])
        rows,_,spectrum = fit(M,b)
        checks[name] = dict(number_complex_equations=len(b),ranks=rows,singular_values=spectrum)
    out = Path(__file__).resolve().parents[1]/"results"/"EX-GC34Z"/"posthoc.json"
    out.write_text(json.dumps(checks,indent=2)+"\n")
    print({name: v["ranks"][1] for name,v in checks.items()},flush=True)


if __name__ == "__main__":
    import sys

    (posthoc if sys.argv[1:] == ["posthoc"] else main)()
