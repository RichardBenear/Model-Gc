"""GC32B, registered blind response tomography and nonlinear controls.
Run: OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python exgc32b_response.py
Outputs JSON to stdout. No optimization or per-species metric fit.
"""
import json, time
import numpy as np
from scipy.linalg import eigh
from scipy.sparse.csgraph import shortest_path


def rel(a, b):
    return float(np.linalg.norm(a-b)/max(np.linalg.norm(b), 1e-15))


def lap(adj):
    return np.diag(adj.sum(axis=1))-adj


def cycle(n):
    a=np.zeros((n,n))
    for i in range(n):
        a[i,(i+1)%n]=a[(i+1)%n,i]=1
    return a


def torus(n):
    return np.kron(cycle(n),np.eye(n))+np.kron(np.eye(n),cycle(n))


def make_graph(B,n):
    blocks=B.reshape(2,n,2,n).transpose(1,3,0,2)
    norm=np.sqrt(np.sum(blocks**2,axis=(2,3)))
    np.fill_diagonal(norm,0)
    return (norm>1e-4*norm.max()).astype(float)


def stable_spectral(A,B):
    a,U=eigh(A)
    sqrtA=(U*np.sqrt(a))@U.T
    w2,V=eigh(sqrtA@B@sqrtA)
    assert min(w2)>0
    w=np.sqrt(w2)
    F=sqrtA@V
    return w,F


def kernel(w,F,t):
    return (F*(np.sin(w*t)/w))@F.T


def estimate_AB(samples, h1=.01, h2=.005):
    # Input is only impulse samples. Odd extension implied by oscillator response.
    def deriv(h):
        q,r=samples[h]
        return (8*q-r)/(6*h),(r-2*q)/h**3
    a1,d1=deriv(h1); a2,d2=deriv(h2)
    A=(4*a2-a1)/3
    D3=(4*d2-d1)/3
    IA=np.linalg.inv(A)
    return A,-IA@D3@IA


def harmonic_case(tag,adj,b,seed):
    n=len(adj)
    p=np.random.default_rng(seed).permutation(n)
    perm=np.r_[p,p+n]; inv=np.argsort(p)
    L=lap(adj)+.01*np.eye(n)
    M=np.array([[1,b],[b,1]])
    A=np.kron(M,np.eye(n))[np.ix_(perm,perm)]
    B=np.kron(M,L)[np.ix_(perm,perm)]
    w,F=stable_spectral(A,B)
    samples={h:(kernel(w,F,h),kernel(w,F,2*h)) for h in [.02,.01,.005]}
    ah,bh=estimate_AB(samples)
    ac,bc=estimate_AB(samples,.02,.01)
    g=make_graph(bh,n)
    canonical=g[np.ix_(inv,inv)]
    # Inference frozen before the dispersion analysis.
    lg=lap(g); ev=eigh(lg,eigvals_only=True)
    distances=shortest_path(g,directed=False,unweighted=True)
    truth_dist=shortest_path(adj,directed=False,unweighted=True)
    qcov=(F*(.5/w))@F.T
    C1=qcov[:n,:n]
    inverse=np.linalg.inv(2*C1)
    staticL=inverse@inverse
    staticG=make_graph(np.kron(np.eye(2),staticL),n)
    # Retarded-derived generator, with predefined exchange parity channels.
    D=ah@bh
    I=np.eye(n)
    slopes=[]
    for sign in (1,-1):
        P=np.vstack((I,sign*I))/np.sqrt(2)
        W=P.T@D@P
        z=eigh((W+W.T)/2,eigvals_only=True)
        # Remove gap and match ordered spectra, including degeneracies.
        slope=np.sqrt(np.dot(ev,z-z[0])/np.dot(ev,ev))
        slopes.append(float(slope))
    ratio=slopes[0]/slopes[1]
    # Off-diagonal covariance is not a sparse adjacency.
    raw=C1.copy(); np.fill_diagonal(raw,0)
    lvals,lvecs=eigh(L)
    exact_cov=(lvecs*(.5/np.sqrt(lvals)))@lvecs.T
    exact_cov=exact_cov[np.ix_(p,p)]
    return dict(graph=tag,n=n,b=b,seed=seed,Aerr=rel(ah,A),Berr=rel(bh,B),
                A_coarse=rel(ac,A), B_coarse=rel(bc,B),
                edge_ok=bool(np.array_equal(canonical,adj)),
                distance_ok=bool(np.array_equal(distances[np.ix_(inv,inv)],truth_dist)),
                static_edge_ok=bool(np.array_equal(staticG,g)),
                raw_cov_density=float(np.sum(abs(raw)>1e-4*np.max(abs(raw)))/(n*(n-1))),
                cov_error=rel(C1,exact_cov),
                qcov_fingerprint=[float(np.trace(C1)),float(np.linalg.norm(C1))],
                E0=float(w.sum()/2),speed_plus=slopes[0],speed_minus=slopes[1],ratio=ratio,
                ratio_error=abs(ratio/((1+b)/(1-b))-1),common=bool(abs(ratio-1)<.01),
                gap=float(ev[1]),volume1=float(np.mean(np.sum(distances<=1,axis=1))),
                volume2=float(np.mean(np.sum(distances<=2,axis=1))))


def coherent_control():
    n=12; L=lap(cycle(n))+.01*np.eye(n)
    M=np.array([[1,.25],[.25,1]])
    A=np.kron(M,np.eye(n)); B=np.kron(M,L)
    w,F=stable_spectral(A,B)
    Fi=np.linalg.inv(F)
    kick=np.zeros(2*n); kick[1]=1
    base=kernel(w,F,1.3)@kick
    out=[]
    for amp in [0,1,2]:
        q=np.zeros(2*n); q[0]=amp
        unforced=(F*np.cos(w*1.3))@Fi@q
        forced=unforced+kernel(w,F,1.3)@kick
        out.append(dict(amplitude=amp,energy=float(q@B@q/2),
                        differential_error=rel(forced-unforced,base)))
    return out


def source_states(n):
    z=np.zeros(n,complex); z[0]=1
    u=np.ones(n,complex)/np.sqrt(n)
    k=np.exp(2j*np.pi*np.arange(n)/n)/np.sqrt(n)
    standing=(k+k.conj())/np.sqrt(2)
    return [
        ("vacuum",0,0,None),("local1",1,0,z),("local2",2,0,z),
        ("uniform1",1,0,u),("plane1",1,.5,k),("standing1",1,.5,standing)]


def quantum_case(n,arm,g,state,return_matrices=False):
    name,num,js,psi=state
    adj=cycle(n); L=lap(adj); hp=L+np.eye(n)
    W=[]
    for i in range(n):
        T=np.zeros((n,n))
        if arm=="onsite": T[i,i]=1
        else:
            for j in np.flatnonzero(adj[i]):
                v=np.eye(n)[i]-np.eye(n)[j]
                T+=np.outer(v,v)
        W.append(T)
    W=np.array(W)
    if num==0:
        energy=0; residual=0; H=hp; Z=np.eye(n)
        density=np.zeros(n)
    elif num==2:
        # Static localized source; exact invariant block, n_s(0)=2.
        energy=2.; residual=0.; H=hp+2*g*W[0]+energy*np.eye(n)
        Z=np.eye(n); density=2*abs(psi)**2
    else:
        hs=np.eye(n)+js*L
        energy=float(np.real(psi.conj()@hs@psi))
        residual=float(np.linalg.norm(hs@psi-energy*psi))
        H=np.kron(hs,np.eye(n))+np.kron(np.eye(n),hp)
        for i in range(n):
            H[i*n:(i+1)*n,i*n:(i+1)*n]+=g*W[i]
        Z=np.kron(psi[:,None],np.eye(n))
        density=abs(psi)**2
    H=H-energy*np.eye(len(H))
    vals,U=eigh(H)
    F=Z.conj().T@U
    def Q(t):
        return (F*np.exp(-1j*vals*t))@F.conj().T
    def moments(h):
        q=Q(h)
        return 1j*(q-q.conj().T)/(2*h),-(q+q.conj().T-2*np.eye(n))/h**2
    T1,S1=moments(.01); T2,S2=moments(.005)
    T=(4*T2-T1)/3; S=(4*S2-S1)/3
    trueT=Z.conj().T@H@Z
    trueS=Z.conj().T@H@H@Z
    expected=hp+g*np.einsum("i,ijk->jk",density,W)
    E=eigh(hp)[0]; V=eigh(hp)[1]
    base=lambda t:(V*np.exp(-1j*E*t))@V.T
    edge=-T.real.copy(); np.fill_diagonal(edge,0)
    delta=edge-adj
    onsite_diag=np.diag(T).real-1-edge.sum(axis=1)
    distances=shortest_path(adj,directed=False,unweighted=True)
    # For localized test the nearest endpoint of each original edge.
    endpoints=np.minimum(distances[0,:,None],distances[0,None,:])
    far=(endpoints>=2)&(adj>0)
    result=dict(arm=arm,g=g,state=name,energy=energy,source_residual=residual,
                T_error=rel(T,trueT),S_error=rel(S,trueS),projected_error=rel(trueT,expected),
                response_t1=rel(Q(1),base(1)),response_t2=rel(Q(2),base(2)),
                moment_variance=float(np.linalg.norm(S-T@T)),
                exact_variance=float(np.linalg.norm(trueS-trueT@trueT)),
                edge_change_max=float(np.max(abs(delta))),
                edge_change_min_on_edge=float(np.min(delta[adj>0])),
                edge_change_norm=float(np.linalg.norm(delta)),
                edge_change_far=float(np.max(abs(delta[far]))),
                nonedge_max=float(np.max(abs(delta[(adj==0)]))),
                onsite_max=float(np.max(abs(onsite_diag))),
                onsite_norm=float(np.linalg.norm(onsite_diag)),
                edge_profile=np.round(delta[0],8).tolist(),
                onsite_profile=np.round(onsite_diag,8).tolist())
    return (result,T,Q(1),Q(2)) if return_matrices else result


def main():
    t0=time.process_time(); wall=time.monotonic()
    graphs=[(f"cycle{n}",cycle(n)) for n in (12,24,48)]
    graphs += [(f"torus{n}",torus(n)) for n in (4,6,8)]
    irregular=cycle(24)
    for i,j in [(0,5),(3,12),(7,19)]: irregular[i,j]=irregular[j,i]=1
    graphs.append(("irregular24",irregular))
    rows=[harmonic_case(tag,a,b,s) for tag,a in graphs for b in (0,.25,.5) for s in (17,61)]
    control=coherent_control()
    gate=all(r["edge_ok"] and r["distance_ok"] and r["Aerr"]<1e-5 and
             r["Berr"]<1e-5 and r["ratio_error"]<1e-4 and r["common"]==(r["b"]==0)
             for r in rows) and all(r["differential_error"]<1e-10 for r in control)
    dims={}
    for prefix in ["cycle","torus"]:
        group=[r for r in rows if r["graph"].startswith(prefix) and r["b"]==0 and r["seed"]==17]
        exponent=np.polyfit(np.log([r["n"] for r in group]),np.log([r["gap"] for r in group]),1)[0]
        dims[prefix]=float(-2/exponent)
    nonlinear=[]
    if gate:
        original=[quantum_case(12,arm,g,state,True) for arm in ["onsite","bond"]
                  for g in [0,.125,.25,.5,1] for state in source_states(12)]
        nonlinear=[r[0] for r in original]
    stage2_gate=bool(nonlinear) and all(r["T_error"]<1e-5 and r["S_error"]<1e-5 and
                     r["source_residual"]<1e-10 and
                     (r["g"]!=0 or r["response_t2"]<1e-10) for r in nonlinear)
    follow=[]
    if stage2_gate:
        # Registered after the original 60 cases, before these additional cases.
        for arm in ["onsite","bond"]:
            for g in [0,.125,.25,.5,1]:
                original_plane=next(r for r in original if r[0]["arm"]==arm and
                                    r[0]["g"]==g and r[0]["state"]=="plane1")
                vals=[quantum_case(12,arm,g,(f"plane_m{m}",1,.5,
                       np.exp(2j*np.pi*m*np.arange(12)/12)/np.sqrt(12)),True)
                      for m in [0,3,6]]
                ref=vals[0]
                for row in [vals[0],original_plane,*vals[1:]]:
                    r=dict(row[0])
                    r.update(first_difference_k0=rel(row[1],ref[1]),
                             response_difference_k0_t1=rel(row[2],ref[2]),
                             response_difference_k0_t2=rel(row[3],ref[3]))
                    follow.append(r)
    follow_gate=bool(follow) and all(r["source_residual"]<1e-10 and r["T_error"]<1e-5
                   and r["S_error"]<1e-5 and r["first_difference_k0"]<1e-5 and
                   (r["g"]!=0 or r["response_t2"]<1e-10) for r in follow)
    print(json.dumps(dict(stage1_gate=gate,stage2_gate=stage2_gate,follow_gate=follow_gate,harmonic=rows,
                          coherent=control,dimension_fit=dims,nonlinear=nonlinear,followup=follow,
                          cpu_seconds=time.process_time()-t0,wall_seconds=time.monotonic()-wall),indent=2))


if __name__=="__main__":
    main()
