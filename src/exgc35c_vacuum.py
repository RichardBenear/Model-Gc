"""Restricted and unrestricted synchronous Y1 vacuum spectrum."""
import json
from pathlib import Path
import exgc35c as C
import numpy as np


def spectrum(n):
    exp=C.Experiment(n,True,True)
    s=exp.system
    rows=[]
    poly_max=0.
    for k,m,M,RC,RB in zip(exp.K,exp.modes,s.AiC,s.RC,s.RB):
        k2=k@k
        Q=M/k2
        eye=np.eye(6)
        zero=np.zeros((6,6))
        G=np.block([[zero,eye],[Q,zero]])
        # Scale velocity by |K| so the operator is dimensionless.
        constraints=np.concatenate([RC,RB*np.sqrt(k2)],axis=-1)
        _,sing,vh=np.linalg.svd(constraints,full_matrices=True)
        rank=int(np.sum(sing>sing[0]*1e-12))
        V=vh[rank:].conj().T
        W=V.conj().T@G@V
        vals=np.linalg.eigvals(W)
        leakage=np.linalg.norm(G@V-V@W)
        defect=np.linalg.norm(Q@Q@(Q+eye))
        rank_q=np.linalg.matrix_rank(Q,tol=1e-11)
        rank_q2=np.linalg.matrix_rank(Q@Q,tol=1e-11)
        rows.append((defect,leakage,max(abs(vals.real)),rank,rank_q,rank_q2))
        poly_max=max(poly_max,defect)
    return dict(n=n,checked_modes=len(rows),max_polynomial_defect=poly_max,
                max_constraint_subspace_leakage=max(r[1] for r in rows),
                max_restricted_real_dimensionless=max(r[2] for r in rows),
                unique_ranks=sorted({r[3:] for r in rows}))


if __name__=="__main__":
    result=[spectrum(n) for n in (24,48)]
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35C"
    path.mkdir(parents=True,exist_ok=True)
    (path/"vacuum.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
