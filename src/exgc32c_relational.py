"""GC32C fixed-budget vacuum energy and link-force screen.
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python exgc32c_relational.py
Prints JSON. No external geometry in the objective or swap proposals.
"""
import json
import time
import itertools
import numpy as np
from scipy.sparse.csgraph import connected_components, shortest_path

MASSES=(.1,1.,5.)

def lap(a):
    return np.diag(a.sum(axis=1))-a

def spec(a):
    return np.maximum(np.linalg.eigvalsh(lap(a)),0)

def energy(a,mu):
    return float(np.sqrt(mu*mu+spec(a)).sum()/2)

def edges(a):
    return np.argwhere(np.triu(a,1)>0)

def components(a):
    return int(connected_components(a,directed=False,return_labels=False))

def describe(a,with_edges=False):
    n=len(a)
    d=a.sum(axis=1)
    c=components(a)
    tr=float(np.trace(a@a@a)/6)
    denom=float(np.sum(d*(d-1))/2)
    r=dict(n=n,edges=int(d.sum()/2),degrees=sorted(map(int,d)),components=c,
           triangles=int(round(tr)),clustering=3*tr/denom if denom else 0.,
           gap=float(spec(a)[1]),lambda_max=float(spec(a)[-1]))
    dist=shortest_path(a,directed=False,unweighted=True)
    r.update(ball=[float(np.mean((dist<=k).sum(axis=1))) for k in (1,2,3)],
             diameter=float(dist.max()) if c==1 else None,
             mean_distance=float(dist.sum()/(n*(n-1))) if c==1 else None)
    if with_edges:r["edge_list"]=edges(a).tolist()
    return r

def exact():
    n=6; pairs=list(itertools.combinations(range(n),2))
    masks=np.arange(1<<len(pairs),dtype=np.uint32)
    bits=((masks[:,None]>>np.arange(len(pairs),dtype=np.uint32))&1).astype(float)
    a=np.zeros((len(masks),n,n))
    ii,jj=np.array(pairs).T
    a[:,ii,jj]=a[:,jj,ii]=bits
    degrees=a.sum(axis=2)
    l=-a.copy()
    ix=np.arange(n); l[:,ix,ix]=degrees
    eig=np.maximum(np.linalg.eigvalsh(l),0)
    connected=eig[:,1]>1e-8
    m=bits.sum(axis=1)
    rows=[]
    for mu in MASSES:
        u=np.sqrt(mu*mu+eig).sum(axis=1)/2
        for tag,choose in [(f"m{j}",m==j) for j in range(16)]+[(f"reg{j}",(degrees==j).all(axis=1)) for j in (2,3)]:
            for conn in (False,True):
                loc=np.flatnonzero(choose & (connected if conn else True))
                if len(loc)==0: continue
                opt=u[loc].min(); ties=loc[np.abs(u[loc]-opt)<=1e-10]
                # Record every tied spectral type by graph invariants; counts separate.
                types={}
                for k in ties:
                    tri=int(round(np.trace(a[k]@a[k]@a[k])/6))
                    key=(tuple(sorted(degrees[k].astype(int))),int((eig[k]<1e-8).sum()),tri,
                         tuple(np.round(eig[k],8)))
                    if key not in types:
                        types[key]=[int(k),0]
                    types[key][1]+=1
                rows.append(dict(mu=mu,group=tag,connected_only=conn,count=len(loc),
                                 minimum=float(opt),mean=float(u[loc].mean()),
                                 maximum=float(u[loc].max()),ties=len(ties),
                                 examples=[dict(count=v[1],**describe(a[v[0]],True)) for v in types.values()]))
        assert np.count_nonzero(abs(u-u.min())<=1e-10)==1 and np.argmin(u)==0
    return rows

def torus(shape):
    n=int(np.prod(shape)); a=np.zeros((n,n))
    for i in range(n):
        x=list(np.unravel_index(i,shape))
        for j,s in enumerate(shape):
            y=x.copy(); y[j]=(y[j]+1)%s
            k=np.ravel_multi_index(y,shape); a[i,k]=a[k,i]=1
    return a

def blocks(n,d):
    size=5 if d==4 else 8
    b=np.ones((size,size))-np.eye(size)
    if d==6:
        for i in range(0,size,2): b[i,i+1]=b[i+1,i]=0
    a=np.kron(np.eye(n//size),b)
    assert np.all(a.sum(axis=1)==d)
    return a

def necklace(a,d):
    size=5 if d==4 else 8; out=a.copy()
    i,j=edges(a[:size,:size])[0]
    for s in range(0,len(a),size):
        out[s+i,s+j]=out[s+j,s+i]=0
        k=(s+size)%len(a)+i
        out[s+j,k]=out[k,s+j]=1
    return out

def propose(a,rng,connected):
    e=edges(a)
    for _ in range(500):
        id1,id2=rng.choice(len(e),2,replace=False)
        i,j=e[id1]; k,l=e[id2]
        if len({i,j,k,l})!=4: continue
        if rng.integers(2): k,l=l,k
        if a[i,k] or a[j,l]:continue
        b=a.copy()
        b[i,j]=b[j,i]=b[k,l]=b[l,k]=0
        b[i,k]=b[k,i]=b[j,l]=b[l,j]=1
        if connected and components(b)!=1:continue
        return b
    raise RuntimeError("no valid two-edge swap found")

def randomized(a,seed):
    rng=np.random.default_rng(seed)
    for _ in range(20*len(a)):a=propose(a,rng,True)
    return a

def search(a,mu,connected,seed):
    rng=np.random.default_rng(seed)
    u=energy(a,mu); start=u; best=u; best_a=a.copy(); accepted=0
    for j in range(800):
        b=propose(a,rng,connected)
        v=energy(b,mu); t=.05-.045*j/199 if j<200 else 0
        if v<=u+1e-12 or (t>0 and rng.random()<np.exp(-(v-u)/t)):
            a,u=b,v; accepted+=1
            if u<best-1e-12:best,best_a=u,a.copy()
    return dict(start_energy=start,end_energy=u,best_energy=best,accepted=accepted,
                attempted=800,best=describe(best_a,True),end=describe(a,True))

def perturb(a,mu,seed):
    rng=np.random.default_rng(seed); e0=energy(a,mu)
    d=np.array([energy(propose(a,rng,True),mu)-e0 for _ in range(128)])
    return dict(lower=int((d< -1e-10).sum()),higher=int((d>1e-10).sum()),
                tied=int((abs(d)<=1e-10).sum()),min_delta=float(d.min()),
                max_delta=float(d.max()))

def weights(a):
    i,j=np.triu_indices(len(a),1)
    return i,j,a[i,j]>0

def mode_response(a,mu):
    n=len(a); K=lap(a)+mu*mu*np.eye(n)
    lam,V=np.linalg.eigh(K); w=np.sqrt(lam)
    i,j,active=weights(a)
    dx=V[i,:]-V[j,:]
    f0=(dx*dx)@(1/(4*w))
    v=np.eye(n)[i[active][0]]-np.eye(n)[j[active][0]]
    eps=1e-5; delta=np.outer(v,v)*eps
    fd0=(np.sqrt(np.linalg.eigvalsh(K+delta)).sum()-np.sqrt(np.linalg.eigvalsh(K-delta)).sum())/(4*eps)
    true0=f0[np.flatnonzero(active)[0]]
    result=[]
    for k in (0,1,n-1):
        force=dx[:,k]**2/(2*w[k])
        u=V[:,k]
        # Fixed one-quantum excess expectation, not adiabatic mode switching.
        def expec(mat):return w[k]/2+float(u@mat@u)/(2*w[k])
        fd=(expec(K+delta)-expec(K-delta))/(2*eps)
        analytic=force[np.flatnonzero(active)[0]]
        result.append(dict(index=k,source_E=float(w[k]),total_E=float(w.sum()/2+w[k]),
                           source_gradient_active_mean=float(force[active].mean()),
                           source_gradient_active_max=float(force[active].max()),
                           source_gradient_absent_mean=float(force[~active].mean()),
                           source_gradient_absent_max=float(force[~active].max()),
                           fixed_state_fd_abs_error=float(abs(fd-analytic))))
    return dict(vacuum_energy=float(w.sum()/2),vac_gradient_active_mean=float(f0[active].mean()),
                vac_gradient_absent_mean=float(f0[~active].mean()),
                vac_fd_abs_error=float(abs(fd0-true0)),vac_fd_rel_error=float(abs(fd0-true0)/abs(true0)),
                modes=result)

def coherent(a,mu):
    n=len(a); K=lap(a)+mu*mu*np.eye(n)
    lam,V=np.linalg.eigh(K); w=np.sqrt(lam)
    ts=np.linspace(0,40,401); i,j,active=weights(a); out=[]
    for arm in ("q","p"):
        q=np.zeros(n);p=np.zeros(n)
        if arm=="q":q[0]=np.sqrt(2/K[0,0])
        else:p[0]=np.sqrt(2)
        Q=V.T@q;P=V.T@p
        cs=np.cos(w[:,None]*ts);sn=np.sin(w[:,None]*ts)
        qt=Q[:,None]*cs+P[:,None]*sn/w[:,None]
        pt=P[:,None]*cs-Q[:,None]*sn*w[:,None]
        e=(pt**2+lam[:,None]*qt**2).sum(axis=0)/2
        state=V@qt
        f=(state[i,:]-state[j,:])**2/2
        out.append(dict(phase=arm,energy=float(e[0]),max_energy_err=float(abs(e-1).max()),
                        initial_active_max=float(f[active,0].max()),initial_absent_max=float(f[~active,0].max()),
                        timeavg_active_mean=float(f[active,:].mean()),
                        timeavg_absent_mean=float(f[~active,:].mean()),
                        timeavg_active_max=float(f[active,:].mean(axis=1).max()),
                        timeavg_absent_max=float(f[~active,:].mean(axis=1).max())))
    return out

def main():
    start=time.process_time(); wall=time.monotonic()
    small=exact()
    refs=[]; searches=[]; perts=[]; mode=[]; coh=[]
    for shape in [(4,5),(5,8),(6,10),(4,4,4)]:
        a=torus(shape); n=len(a); d=int(a.sum(axis=1)[0]); family=f"N{n}d{d}"
        b=blocks(n,d); c=necklace(b,d)
        graphs={"lattice":a,"blocks":b,"necklace":c,"random17":randomized(a,17),"random61":randomized(a,61)}
        assert all(np.all(x.sum(axis=1)==d) for x in graphs.values())
        for name,x in graphs.items():
            p=np.random.default_rng(701).permutation(n)
            assert np.max(abs(spec(x)-spec(x[np.ix_(p,p)])))<1e-10
            refs.append(dict(family=family,name=name,**describe(x,True),
                             energies={str(mu):energy(x,mu) for mu in MASSES}))
        for mu in (.1,1.):
            perts.append(dict(family=family,mu=mu,**perturb(a,mu,1234+100*len(perts))))
            for conn in (True,False):
                for name in ("random17","lattice","necklace"):
                    r=search(graphs[name],mu,conn,1234+100*len(searches))
                    searches.append(dict(family=family,mu=mu,connected_only=conn,start_graph=name,**r))
        for mu in MASSES:
            for name in ("lattice","necklace","random17"):
                x=graphs[name]
                mode.append(dict(family=family,graph=name,mu=mu,**mode_response(x,mu)))
                coh.append(dict(family=family,graph=name,mu=mu,controls=coherent(x,mu)))
    print(json.dumps(dict(exact=small,references=refs,searches=searches,perturbations=perts,
                          mode_response=mode,coherent=coh,cpu_seconds=time.process_time()-start,
                          wall_seconds=time.monotonic()-wall),separators=(",",":")))

if __name__=="__main__":main()
