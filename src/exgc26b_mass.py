"""EX-GC26B: mass calibration of the inherited massless Skyrme lattice."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS'): os.environ[key]='1'
import argparse, hashlib, json, time, platform
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates, spline_filter
import scipy
import exgc26a_skyrmion as ref
import exgc26a2_qualify as a2

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/EX-GC26B'; RAW=ROOT/'local_raw/EX-GC26B'; NORM=ref.NORM
def save(name,rec):
    OUT.mkdir(parents=True,exist_ok=True);p=OUT/(name+'.json');tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(rec,indent=1,allow_nan=False)+'\n');tmp.replace(p)
def path(h,L):return RAW/f'rest_h{h:g}_L{L:g}.npy'
def dot(a,b):return np.sum(a*b,axis=0)
def fields(G,p,v):
    """Energy, canonical momentum, continuum spatial stress sampled on inherited stencil."""
    H=np.zeros(p.shape[1:]); P=np.zeros((3,)+H.shape); S=np.zeros((3,3)+H.shape)
    for s in (1,-1):
        d,D=G.strain(p,s);tr=sum(D[i,i] for i in range(3));vv=dot(v,v);vd=np.array([dot(v,x) for x in d]);K=vv*(1+tr)-np.sum(vd**2,axis=0);e=tr+.5*(tr**2-np.sum(D**2,axis=(0,1)));mask=G.mask[s]
        H+=.5*(K+e)*mask
        for i in range(3):
            P[i]+=-((1+tr)*vd[i]-sum(D[i,j]*vd[j] for j in range(3)))*mask
            for j in range(3):
                # T_ji = -dL/d(phi_,i) dot phi_,j + delta_ij L
                S[i,j]+=((1+tr-vv)*D[i,j]-sum(D[i,k]*D[j,k] for k in range(3))+vd[i]*vd[j]+.5*(K-e)*(i==j))*mask
    return H,P,S
def static_tensor(G,p):
    z=np.zeros_like(p);H,P,S=fields(G,p,z);fac=G.h**3/NORM;E=H.sum()*fac;stress=S.sum(axis=(2,3,4))*fac
    d=[ref.tangent(p,np.gradient(p,G.h,axis=i+1))*G.inner for i in range(3)];M=np.zeros((3,3));matched=np.zeros((3,3))
    for s in (1,-1):
        dp,D=G.strain(p,s);tr=sum(D[i,i] for i in range(3));dd=np.array([[dot(x,y) for y in dp] for x in d]);
        for i in range(3):
            for j in range(3):
                M[i,j]+=np.sum(((1+tr)*dot(d[i],d[j])-sum(dd[i,k]*dd[j,k] for k in range(3)))*G.mask[s])*fac
                matched[i,j]+=np.sum(((1+tr)*D[i,j]-sum(D[i,k]*D[j,k] for k in range(3)))*G.mask[s])*fac
    return dict(E=E,stress=stress.tolist(),M=M.tolist(),M_over_E=(M/E).tolist(),eigen_over_E=(np.linalg.eigvalsh(M)/E).tolist(),matched_M_over_E=(matched/E).tolist(),stress_identity_max=float(np.max(abs(matched-stress-E*np.eye(3)))))
def prepare(h,L,maxiter,source=None):
    t=time.process_time();G=ref.Grid(h,L);RAW.mkdir(parents=True,exist_ok=True)
    if source:
        p=np.load(source);assert p.shape==(4,G.n,G.n,G.n);e,g=a2.evaluate(G,p,0);info=dict(E=float(e.sum()),ratio=float(e[0]/e[1]),projected_gradient_max=float(abs(ref.tangent(p,g)*G.inner).max()),reused=True,source_sha256=hashlib.sha256(Path(source).read_bytes()).hexdigest())
    else:p,info=a2.optimize(G,G.hedgehog(ref.PROF),0,maxiter)
    np.save(path(h,L),p);info.update(h=h,L=L,actual_L=float(G.ax[-1]),geometric=a2.geometric(p),tensor=static_tensor(G,p),sha256=hashlib.sha256(path(h,L).read_bytes()).hexdigest(),cpu_total=time.process_time()-t);save(f'static_h{h:g}_L{L:g}',info);print('static',h,L,info['E'],info['tensor']['eigen_over_E'],info['cpu_total'],flush=True)
def init_boost(G,p,u,n,t=0.):
    gamma=1/np.sqrt(1-u*u);X=np.array(G.X);long=np.einsum('i,ixyz->xyz',n,X);coords=(X+n[:,None,None,None]*((gamma-1)*long-gamma*u*t))/G.h+(G.n-1)/2
    # Differentiate the same cubic interpolant in the rest-frame direction.
    eps=1e-4;delta=eps*n[:,None,None,None]/G.h;out=np.empty_like(p);der=np.empty_like(p)
    for k in range(4):
        coeff=spline_filter(p[k],order=3,mode='constant');sample=lambda x:map_coordinates(coeff,x,order=3,mode='constant',cval=float(k==0),prefilter=False)
        out[k]=sample(coords);der[k]=(sample(coords+delta)-sample(coords-delta))/(2*eps)
    norm=np.sqrt(dot(out,out));q=G.fix(out);v=ref.tangent(q,-gamma*u*der/norm)*G.inner
    return q,v
def region(G,H,P,S,R):
    # closed cube, trapezoid weights; faces snapped to nearest grid plane
    k=int(np.argmin(abs(G.ax+R)));j=G.n-1-k;sl=slice(k,j+1);w=np.ones(j-k+1);w[[0,-1]]=.5; W=w[:,None,None]*w[None,:,None]*w[None,None,:];face=w[:,None]*w[None,:];fac=G.h**3/NORM
    p=np.sum(P[:,sl,sl,sl]*W,axis=(1,2,3))*fac;flux=np.zeros(3)
    for i in range(3):
        a=np.moveaxis(S[i],i+1,1);flux+=np.sum((a[:,j,sl,sl]-a[:,k,sl,sl])*face,axis=(1,2))*G.h**2/NORM
    h=H[sl,sl,sl]*W;cen=np.array([np.sum(x[sl,sl,sl]*h)/max(h.sum(),1e-300) for x in G.X]);return dict(R=float(G.ax[j]),P=p.tolist(),flux=flux.tolist(),centroid=cen.tolist())
def measure(G,p,v,t):
    H,P,S=fields(G,p,v);fac=G.h**3/NORM;return dict(t=t,H=float(H.sum()*fac),P=(P.sum(axis=(1,2,3))*fac).tolist(),centroid=[float(np.sum(H*x)/H.sum()) for x in G.X],regions=[region(G,H,P,S,R) for R in (3,5)])
def step(G,p,v,dt):
    a=a2.acceleration;k1v=a(G,p,v);k1p=v;k2p=v+dt*k1v/2;k2v=a(G,G.fix(p+dt*k1p/2),k2p);k3p=v+dt*k2v/2;k3v=a(G,G.fix(p+dt*k2p/2),k3p);k4p=v+dt*k3v;k4v=a(G,G.fix(p+dt*k3p),k4p);q=G.fix(p+dt*(k1p+2*k2p+2*k3p+k4p)/6);return q,ref.tangent(q,v+dt*(k1v+2*k2v+2*k3v+k4v)/6)*G.inner
def summary(rows,E0,n,u):
    t=np.array([r['t'] for r in rows]);P=np.array([r['P'] for r in rows]);H=np.array([r['H'] for r in rows]);C=np.array([r['centroid'] for r in rows]);vel=np.linalg.lstsq(np.stack([np.ones_like(t),t],axis=1),C,rcond=None)[0][1];speed=float(vel@n);g=1/np.sqrt(1-speed**2);gamma=1/np.sqrt(1-u*u);out=dict(velocity=vel.tolist(),speed=speed,H_over_gamma_E0=float(H.mean()/(gamma*E0)),P_over_gamma_u_E0=float(P.mean(0)@n/(gamma*u*E0)) if u else None,M_over_E0=float(P.mean(0)@n/(g*speed*E0)) if u else None,H_over_measured_gamma_E0=float(H.mean()/(g*E0)),max_dP_E0=float(np.max(np.linalg.norm(P-P[0],axis=1))/E0),max_dH_H0=float(np.max(abs(H/H[0]-1))),max_transverse_P_E0=float(np.max(np.linalg.norm(P-(P@n)[:,None]*n,axis=1))/E0),regions=[])
    for k in range(2):
        ps=np.array([r['regions'][k]['P'] for r in rows]);fs=np.array([r['regions'][k]['flux'] for r in rows]);cs=np.array([r['regions'][k]['centroid'] for r in rows]);integ=np.r_[np.zeros((1,3)),np.cumsum((fs[1:]+fs[:-1])*np.diff(t)[:,None]/2,axis=0)];res=ps-ps[0]+integ;out['regions'].append(dict(R=rows[0]['regions'][k]['R'],raw_max=float(np.max(np.linalg.norm(ps-ps[0],axis=1))/E0),residual_max=float(np.max(np.linalg.norm(res,axis=1))/E0),end_flux_over_E0=(integ[-1]/E0).tolist(),end_residual_over_E0=(res[-1]/E0).tolist(),speed=float(np.linalg.lstsq(np.stack([np.ones_like(t),t],axis=1),cs,rcond=None)[0][1]@n)))
    return out
def boost(h,L,u,axis,dtfac,T):
    start=time.process_time();G=ref.Grid(h,L);p=np.load(path(h,L));E0=G.energy(p);n=np.array(axis,float);n/=np.linalg.norm(n);p,v=init_boost(G,p,u,n) if u else (p,np.zeros_like(p));degree0=a2.geometric(p);dt=h*dtfac;steps=round(T/dt);every=max(1,round(.08/dt));rows=[]
    for i in range(steps+1):
        if i%every==0 or i==steps:rows.append(measure(G,p,v,i*dt))
        if i<steps:p,v=step(G,p,v,dt)
    rec=dict(h=h,L=L,u=u,axis=axis,dt=dt,E0=E0,rows=rows,summary=summary(rows,E0,n,u),degree_initial=degree0,degree_final=a2.geometric(p),cpu_s=time.process_time()-start);tag=f'h{h:g}_L{L:g}_u{u:g}_axis{"".join(map(str,axis))}_dt{dtfac:g}';save(tag,rec);np.save(RAW/f'{tag}_end.npy',p);print(tag,rec['summary'],rec['cpu_s'],flush=True)
def validate():
    G=ref.Grid(.2,1);p=G.hedgehog(ref.PROF);v=ref.tangent(p,np.random.default_rng(26).normal(scale=.05,size=p.shape))*G.inner;H,P,S=fields(G,p,v);fac=G.h**3/NORM;hr,pr=ref.hamiltonian(G,p,v);e=abs(H.sum()*fac-hr);pe=float(np.max(abs(P.sum(axis=(1,2,3))*fac-pr)));st=static_tensor(G,p);vac=fields(G,G.vacuum(),np.zeros_like(p));n=np.ones(3)/np.sqrt(3);q,w=init_boost(G,p,.3,n);eps=1e-5;fd=(init_boost(G,p,.3,n,eps)[0]-init_boost(G,p,.3,n,-eps)[0])/(2*eps);interior=G.r<.5;de=float(abs(fd[:,interior]-w[:,interior]).max());out=dict(H_error=e,P_error=pe,stress_identity=st['stress_identity_max'],vacuum_max=max(float(abs(a).max()) for a in vac),boost_derivative_interior_error=de,python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__,platform=platform.platform());save('validation',out);print(out,flush=True);assert max(e,pe,st['stress_identity_max'],out['vacuum_max'])<1e-10 and de<1e-5
def main():
    ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['validate','static','boost']);ap.add_argument('--h',type=float,default=.2);ap.add_argument('--L',type=float,default=7);ap.add_argument('--maxiter',type=int,default=350);ap.add_argument('--source');ap.add_argument('--axis',nargs=3,type=int,default=[1,0,0]);ap.add_argument('--u',type=float,default=.3);ap.add_argument('--dtfac',type=float,default=.05);ap.add_argument('--T',type=float,default=.8);a=ap.parse_args();a2.load()
    if a.phase=='validate':validate()
    elif a.phase=='static':prepare(a.h,a.L,a.maxiter,a.source)
    else:boost(a.h,a.L,a.u,a.axis,a.dtfac,a.T)
if __name__=='__main__':main()
