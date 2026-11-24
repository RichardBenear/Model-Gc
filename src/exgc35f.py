"""GC35F: versioned local events for staggered linear Y1 and a scalar."""
import hashlib
import json
from itertools import product
from pathlib import Path
import random
import time

import exgc34p as P
import numpy as np

N=5
H=.5
AXES=np.eye(3,dtype=int)
CELLS=list(product(range(N),repeat=3))
SHAPE=(N,)*3
PATCH=list(product((-1,0,1),repeat=3))
ZERO=(0,0,0)
PAIRS=P.N.PAIRS
OFFSETS=np.asarray(P.N.OFFS)


def translate(x,dx):
    return tuple((a+b)%N for a,b in zip(x,dx))


class Host:
    def __init__(self,dt):
        self.dt=dt
        waves=np.meshgrid(*([2*np.pi*np.fft.fftfreq(N,d=H)]*3),indexing="ij")
        k=np.stack(waves,axis=-1).reshape(-1,3)
        self.K=2*np.sin(k*H/2)/H
        self.S=P.System(self.K)
        self.phase=np.exp(1j*k@(H*OFFSETS).T)
        spatial=self.S.AiC*self.phase[:,4:,None]*self.phase[:,None,4:].conj()
        stencil=P.sf.ifftn(spatial.reshape((*SHAPE,6,6)),axes=(0,1,2))
        self.ops=[(dx,stencil[tuple((-np.asarray(dx))%N)].real) for dx in PATCH]
        support={tuple((-np.asarray(dx))%N) for dx in PATCH}
        self.stencil_error=max(float(np.max(abs(stencil.imag))),
                               max(float(np.max(abs(stencil[x]))) for x in CELLS if x not in support))
        self.force=-self.S.Ainv[0]*(2*P.KAP)
        phase_force=self.S.Ainv*self.phase[:,4:,None]*self.phase[:,None,4:].conj()
        self.force_local_error=float(np.max(abs(phase_force-self.S.Ainv[0])))
        points=np.random.default_rng(3506).uniform(-np.pi,np.pi,size=(128,3))
        Ktest=2*np.sin(points/2)/H
        phases=np.exp(1j*points@OFFSETS[4:].T)
        expected=P.System(Ktest).AiC*phases[:,:,None]*phases[:,None,:].conj()
        predicted=sum(np.exp(1j*(points@dx))[:,None,None]*w for dx,w in self.ops)
        self.off_grid_error=float(np.max(abs(predicted-expected)))
        assert max(self.stencil_error,self.force_local_error,self.off_grid_error)<1e-11

    def flat_lap(self,phi):
        return P.N.lap(phi,H)-phi

    def fourier(self,field,part):
        values=P.sf.fftn(field,axes=(0,1,2)).reshape(N**3,-1)/N**3
        return values*self.phase[:,part].conj()

    def inverse(self,values,part):
        return P.sf.ifftn((values*self.phase[:,part]*N**3).reshape((*SHAPE,-1)),axes=(0,1,2)).real

    def source(self,phi,v):
        lapse,b,d,o=P.O.msources2(phi,v,H)
        return np.stack((lapse,*b,*d,*o),axis=-1)*2*P.KAP/H**3

    def initial(self):
        x=np.indices(SHAPE)
        phi=.3*np.cos(2*np.pi*(x[0]+x[1])/N)+.15*np.sin(2*np.pi*x[2]/N)
        vel=.11*np.cos(2*np.pi*(x[0]+x[2])/N)+.07*np.sin(2*np.pi*x[1]/N)
        f=self.fourier(self.source(phi,vel),slice(None))
        mh,mhv,_=P.System(self.K[1:]).init(f[1:])
        mh=np.vstack([np.zeros((1,6)),mh])
        mhv=np.vstack([np.zeros((1,6)),mhv])
        start=self.inverse(mh,slice(4,None))
        startv=self.inverse(mhv,slice(4,None))
        accel=self.reference_force(start,phi,vel)
        return {("s",0):phi,("s",-1):phi-self.dt*vel+.5*self.dt**2*self.flat_lap(phi),
                ("m",0):start,("m",-1):start-self.dt*startv+.5*self.dt**2*accel}

    def reference_force(self,m,s,v):
        coef=self.fourier(m,slice(4,None))
        f=self.fourier(self.source(s,v),slice(None))
        return self.inverse(self.S.acc(coef,coef*0,f[:,4:]),slice(4,None))

    def synchronous(self,initial,depth):
        fields=initial.copy()
        for n in range(depth):
            s=fields["s",n]
            m=fields["m",n]
            prev=fields["s",n-1]
            sn=2*s-prev+self.dt**2*self.flat_lap(s)
            v=(sn-prev)/(2*self.dt)
            mn=2*m-fields["m",n-1]+self.dt**2*self.reference_force(m,s,v)
            fields["s",n+1],fields["m",n+1]=sn,mn
        return fields

    def event(self,node,fields,stale=False):
        kind,n,x=node
        def read(key,depth,dx=ZERO):
            return fields[key,depth,translate(x,dx)]
        if kind=="s":
            center=read("s",n-1)
            grad=sum(read("s",n-1,a)+read("s",n-1,-a)-2*center for a in AXES)/H**2
            return 2*center-read("s",n-2)+self.dt**2*(grad-center)
        center=read("s",n-1)
        prev=read("s",n-2)
        future=(center if stale and ("s",n,x) not in fields else read("s",n))
        v=(future-prev)/(2*self.dt)
        gp=[(read("s",n-1,a)-center)/H for a in AXES]
        gm=[(center-read("s",n-1,-a))/H for a in AXES]
        diagonal=[.25*(v*v-center*center)-.125*sum(gp[b]**2+gm[b]**2 for b in range(3) if b!=a)
                  +.25*gp[a]*gm[a] for a in range(3)]
        shear=[]
        for a,b in PAIRS[3:]:
            c=read("s",n-1,AXES[a]+AXES[b])
            za,zb=read("s",n-1,AXES[a]),read("s",n-1,AXES[b])
            shear.append(.25*(gp[a]+(c-zb)/H)*(gp[b]+(c-za)/H))
        acc=sum(matrix@read("m",n-1,dx) for dx,matrix in self.ops)
        acc+=self.force@np.array(diagonal+shear)
        return 2*read("m",n-1)-read("m",n-2)+self.dt**2*acc

    def graph(self,depth,stale):
        counts,followers={},{}
        for n in range(1,depth+1):
            for x in CELLS:
                scalar=("s",n,x)
                gravity=("m",n,x)
                for node in (scalar,gravity):
                    kind=node[0]
                    shifts=[ZERO,*AXES,*(-AXES)] if kind=="s" else PATCH
                    parents={(kind,n-1,translate(x,v)) for v in shifts if n>1}
                    if n>2:
                        parents.add((kind,n-2,x))
                    if kind=="m":
                        if n>1:
                            parents.update(("s",n-1,translate(x,v)) for v in PATCH)
                        if n>2:
                            parents.add(("s",n-2,x))
                        if not stale:
                            parents.add(scalar)
                    counts[node]=len(parents)
                    for p in parents:
                        followers.setdefault(p,[]).append(node)
        return counts,followers

    def asynchronous(self,initial,graph,depth,seed,stale=False):
        counts,followers=graph
        counts=counts.copy()
        fields={(kind,n,x):v[x] for (kind,n),v in initial.items() for x in CELLS}
        ready=[key for key,val in counts.items() if val==0]
        progress=[0]*(depth+1)
        witness=None
        rng=random.Random(seed)
        while ready:
            index=rng.randrange(len(ready))
            node=ready[index]
            ready[index]=ready[-1]
            ready.pop()
            kind,n,x=node
            fields[node]=self.event(node,fields,stale)
            progress[n]+=1
            if witness is None and n>1 and progress[n-1]<2*N**3:
                witness=dict(event=[kind,n,list(x)],preceding_finished=progress[n-1],preceding_total=2*N**3)
            for child in followers.get(node,[]):
                counts[child]-=1
                if counts[child]==0:
                    ready.append(child)
        assert sum(progress)==len(counts)
        result=initial.copy()
        for n in range(1,depth+1):
            for kind,shape in (("s",SHAPE),("m",(*SHAPE,6))):
                result[kind,n]=np.array([fields[kind,n,x] for x in CELLS]).reshape(shape)
        return result,witness

    def constraints(self,fields,depth):
        rows=[]
        for n in range(depth):
            s=fields["s",n]
            v=(fields["s",n+1]-fields["s",n-1])/(2*self.dt)
            d=(fields["m",n+1]-fields["m",n-1])/(2*self.dt)
            f=self.fourier(self.source(s,v),slice(None))
            r,a=self.S.cres(self.fourier(fields["m",n],slice(4,None)),self.fourier(d,slice(4,None)),f)
            dd=(fields["m",n+1]-2*fields["m",n]+fields["m",n-1])/self.dt**2
            a_event=self.fourier(dd,slice(4,None))
            # Check the original constraint equations using the event acceleration.
            r-=P.mv(self.S.Ach,a_event-a)
            energy=float(np.linalg.norm(r[1:,0])/np.linalg.norm(f[1:,:4]))
            momentum=float(np.linalg.norm(r[1:,1:])/np.linalg.norm(f[1:,:4]))
            rows.append(dict(t=n*self.dt,energy=energy,momentum=momentum,
                             combined=float(np.hypot(energy,momentum)),
                             spatial_max=float(np.max(abs(a_event-a)))))
        return rows


def run(dt):
    host=Host(dt)
    depth=round(.24/dt)+1
    initial=host.initial()
    ref=host.synchronous(initial,depth)
    graph=host.graph(depth,False)
    first,witness=host.asynchronous(initial,graph,depth,0)

    def checksum(f):
        return hashlib.sha256(b"".join(f[key].tobytes() for key in sorted(f))).hexdigest()

    expected=checksum(first)
    matches=1
    overlaps=int(witness is not None)
    for seed in range(1,24):
        other,w=host.asynchronous(initial,graph,depth,seed)
        matches+=all(first[key].tobytes()==other[key].tobytes() for key in first)
        overlaps+=w is not None
    badgraph=host.graph(depth,True)
    bad,_=host.asynchronous(initial,badgraph,depth,0,stale=True)
    bad2,_=host.asynchronous(initial,badgraph,depth,1,stale=True)
    def difference(a,b,kind):
        return max(float(np.max(abs(a[key]-b[key]))) for key in a if key[0]==kind)
    return dict(dt=dt,events=len(graph[0]),matched=matches,schedules=24,sha256=expected,
                overlap_count=overlaps,witness=witness,scalar_error=difference(first,ref,"s"),
                metric_error=difference(first,ref,"m"),negative_error=difference(first,bad,"m"),
                negative_order_difference=difference(bad,bad2,"m"),stencil_error=host.stencil_error,
                force_local_error=host.force_local_error,off_grid_error=host.off_grid_error,constraints=host.constraints(first,depth))


if __name__=="__main__":
    start=time.process_time()
    result=dict(runs=[run(.04),run(.02)],cpu_seconds=time.process_time()-start)
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35F"
    path.mkdir(parents=True,exist_ok=True)
    (path/"results.json").write_text(json.dumps(result,indent=2)+"\n")
    print(result,flush=True)
