"""GC27A/B: static vacuum-scale sensitivity, no gravitational dynamics."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[k]='1'
import hashlib,json,time,resource
from pathlib import Path
import numpy as np
from scipy.ndimage import affine_transform
import exgc26e1 as e
import exgc26b2 as b

ROOT=e.ROOT;OUT=ROOT/'results/EX-GC27AB';RAW=ROOT/'local_raw/EX-GC27AB';NORM=e.NORM
START=time.process_time();CAP=5500
ACASES=((2,0,0),(0,2,0),(2,-2,0),(2,2,0));BCASES=((2,0,4),(2,0,0),(0,0,1))

def save(tag,r):
    OUT.mkdir(parents=True,exist_ok=True);p=OUT/(tag+'.json');q=p.with_suffix('.tmp');q.write_text(json.dumps(r,indent=2,allow_nan=False)+'\n');q.replace(p)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gate():
    if time.process_time()-START>CAP:raise RuntimeError('Budget launch stop; completed checkpoints retained')
def coef(eta,exps,C0):return np.exp(eta*np.array(exps))*[1,1,C0]
def weighted(G,p,co):
    A,B,C=co;parts,g=e.evaluate(G,p,np.full(G.r.shape,A/B));parts*=B;g*=B;fac=G.h**3/NORM;V=C*float(np.sum(1-p[0])*fac);g[0]-=C*fac
    return np.r_[parts,V],g
def optimize(G,p,co,it):
    gate();e.a2.evaluate=lambda GG,q,mode=0:(lambda z:(np.array([z[0][0]+z[0][2],z[0][1]]),z[1]))(weighted(GG,q,co));return e.a2.optimize(G,p,0,it)
def metrics(G,p,co,degree=False):
    v,g=weighted(G,p,co);w=(1-p[0])**2;Q=w.sum();radius=float(np.sqrt(np.sum(w*G.r**2)/Q)) if Q else 0.;E=float(v.sum());r=dict(parts=v.tolist(),E=E,R=radius,gradient=float(np.max(np.abs(e.ref.tangent(p,g)*G.inner))),virial=float((v[0]-v[1]+3*v[2])/E) if E else 0.)
    if degree:
        r['degree']=e.a2.geometric(p);H=np.zeros(G.r.shape)
        for s in __import__('itertools').product((-1,1),repeat=3):
            mask,d,D=b.orient(G,p,s);tr=sum(D[i,i] for i in range(3));H+=(co[0]*tr+co[1]*.5*(tr**2-np.sum(D**2,axis=(0,1))))*mask/8
        H+=co[2]*(1-p[0]);wall=np.max(np.abs(np.array(G.X)),axis=0)>G.ax[-1]-1.5;r['wall_energy_fraction']=float(H[wall].sum()/H.sum()) if E else 0.
    return r
def rescale(G,p,lam):
    center=(G.n-1)/2;out=np.stack([affine_transform(p[i],np.eye(3)/lam,offset=np.full(3,center*(1-1/lam)),order=3,mode='constant',cval=float(i==0)) for i in range(4)]);return G.fix(out)
def tag(h,B):return f'h{h:g}_B{B}'
def baseline(h,B):
    G=e.ref.Grid(h,7);t=tag(h,B);gate();prior=None
    if B==1:
        fname=b.restpath(h,7);meta=json.loads((b.OUT/f'static_h{h:g}_L7.json').read_text());assert sha(fname)==meta['sha256'];p=np.load(fname);info=dict(reused=True,source_sha256=sha(fname))
    else:p,info=optimize(G,e.rational_map(G,2,1.3),[1,1,0],500)
    path=RAW/f'A_{t}.npy';np.save(path,p);m=metrics(G,p,[1,1,0],True);m.update(h=h,B=B,sha256=sha(path),optimization=info,cpu_total=time.process_time()-START);save('A_base_'+t,m);assert m['degree']['degrees']==[-B]*3;print('A base',t,m['E'],m['parts'],m['gradient'],flush=True);return G,p,m
def sensitivities(G,p,base,C0,exp_set,phase,B):
    E=base['E'];v=np.array(base['parts']);out=[]
    for ex in exp_set:
        k=float(np.dot(v,ex)/E);rec=dict(exponents=ex,k=k,s=k*E,points=[])
        for eta in ((-.05,-.02,.02,.05) if phase=='A' else (-.02,.02)):
            gate();co=coef(eta,ex,C0);fixed=metrics(G,p,co);row=dict(eta=eta,fixed=fixed)
            if phase=='A':
                lam=np.exp((ex[1]-ex[0])*eta/2);q=rescale(G,p,lam) if lam!=1 else p;row['dilation']=metrics(G,q,co);row['lambda']=float(lam)
                if G.h==.2 and abs(eta)==.02 and ex in ((2,0,0),(2,-2,0)):
                    q,info=optimize(G,q,co,120);row['relaxed']=metrics(G,q,co,True);row['optimization']=info
            elif G.h==.2 and ex==(2,0,0):
                q,info=optimize(G,p,co,120);row['relaxed']=metrics(G,q,co,True);row['optimization']=info
            rec['points'].append(row)
        for kind in ('fixed','dilation','relaxed'):
            vals={r['eta']:r[kind]['E'] for r in rec['points'] if kind in r};
            if -.02 in vals and .02 in vals:rec['k_fd_'+kind]=(vals[.02]-vals[-.02])/(.04*E)
        out.append(rec);save(phase+'_sensitivity_'+tag(G.h,B),dict(h=G.h,B=B,base=base,cases=out,cpu_total=time.process_time()-START));print(phase,tag(G.h,B),ex,k,flush=True)
def phaseB(h,B,C0):
    G=e.ref.Grid(h,7);t=tag(h,B);p=np.load(RAW/f'A_{t}.npy');p,info=optimize(G,p,[1,1,C0],350);path=RAW/f'B_{t}.npy';np.save(path,p);m=metrics(G,p,[1,1,C0],True);m.update(h=h,B=B,C0=C0,sha256=sha(path),optimization=info,cpu_total=time.process_time()-START);save('B_base_'+t,m);assert m['degree']['degrees']==[-B]*3;print('B base',t,m['E'],m['parts'],m['gradient'],flush=True);sensitivities(G,p,m,C0,BCASES,'B',B)
def validation():
    G=e.ref.Grid(.25,1);rng=np.random.default_rng(27);p=G.fix(G.vacuum()+.2*rng.normal(size=(4,G.n,G.n,G.n)));q=e.ref.tangent(p,rng.normal(size=p.shape))*G.inner;q/=np.linalg.norm(q);co=[1.17,.89,.31];parts,g=weighted(G,p,co);ep=1e-5;fd=(weighted(G,G.fix(p+ep*q),co)[0].sum()-weighted(G,G.fix(p-ep*q),co)[0].sum())/(2*ep);an=float(np.sum(g*q));vac=metrics(G,G.vacuum(),co,True);err=abs(fd-an)/max(abs(an),1e-12);save('validation',dict(fd=fd,analytic=an,error=err,vacuum=vac));assert err<1e-6 and vac['E']==0 and vac['gradient']==0
def main():
    resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));RAW.mkdir(parents=True,exist_ok=True);e.load();validation()
    for h in (.2,.15):
        for B in (1,2):
            G,p,m=baseline(h,B)
            if B==1:
                refl=metrics(G,p[:,::-1,:,:].copy(),[1,1,0],True);save('A_parity_'+tag(h,B),refl)
            sensitivities(G,p,m,0,ACASES,'A',B)
    m=json.loads((OUT/'A_base_h0.2_B1.json').read_text());C0=2*(.3/m['R'])**2;save('B_coefficient',dict(C0=C0,R_B1=m['R'],mass=.3/m['R'],cpu_total=time.process_time()-START))
    for h in (.2,.15):
        for B in (1,2):phaseB(h,B,C0)
    save('completion',dict(cpu_s=time.process_time()-START,complete=True));print('DONE',time.process_time()-START,flush=True)
if __name__=='__main__':main()
