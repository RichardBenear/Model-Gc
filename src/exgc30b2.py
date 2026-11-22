"""EX-GC30B2: band geometry versus winding in explicit periodic Weyl models."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,time,platform,resource,itertools
from pathlib import Path
import numpy as np

I=np.eye(2); TSH=np.array([[1.1,.15],[.05,.9]]); E=np.eye(3)
DIR=np.array([*E,*[(E[i]+s*E[j])/np.sqrt(2) for i in range(3) for j in range(i+1,3) for s in (1,-1)],np.ones(3)/np.sqrt(3)])
X=np.array([[x*x,y*y,z*z,2*x*y,2*x*z,2*y*z] for x,y,z in DIR])
def mat(a):
    x,y,z,xy,xz,yz=a;return np.array([[x,xy,xz],[xy,y,yz],[xz,yz,z]])

def cases():
    def S(name,kind='protected',m=0.,a=.5,t=1.,T=I,u=0.,v=0.):
        return dict(name=name,arm='S',kind=kind,blocks=[dict(m=m,a=a,t=t,T=T,u=u,v=v)])
    def F(name,kind='common',shift=0.,z=1.,T=I,u=0.,v=0.,relative=1.,tuned=True):
        return dict(name=name,arm='F',kind=kind,blocks=[dict(m=m+shift,a=None,t=z/(np.sqrt(1-m*m) if tuned else 1.),T=T*(relative if j else 1.),u=u,v=v) for j,m in enumerate((.2,.6))])
    yield S('baseline');yield S('transverse',T=1.1*I);yield S('axial',t=1.1);yield S('shear',T=TSH);yield S('combined',T=TSH,t=1.1)
    for a in (.4,.6):yield S(f'root_separation_{a}',a=a)
    for m in (-.1,-.05,-.025,.025,.05,.1):yield S(f'offset_{m}','common',m=m)
    yield S('sine_tilt','common_drift',u=.1);yield S('cosine_tilt','common',v=.1);yield S('offset_combined','common',m=.05,T=TSH,t=1.1)
    yield F('natural','untuned',tuned=False);yield F('tuned','calibration');yield F('shear',T=TSH);yield F('axial',z=1.1)
    for m in (-.1,-.05,-.025,.025,.05,.1):yield F(f'offset_{m}',shift=m)
    yield F('offset_combined',shift=.05,T=TSH,z=1.1);yield F('relative_transverse','relative',relative=1.1)
    yield F('sine_tilt',u=.1);yield F('cosine_tilt',v=.1)

def dvec(b,k):
    x,y,z=k; f=np.cos(z)-b['m']; f=f if b['a'] is None else f*f-b['a']**2
    return np.r_[b['T']@np.sin([x,y]),b['t']*f+3*(2-np.cos(x)-np.cos(y))]
def d0(b,k):return b['u']*np.sin(k[0])+b['v']*np.cos(k[2])
def energies(b,k):
    x,y,z=dvec(b,k);s=d0(b,k);return np.linalg.eigvalsh(np.array([[s+z,x-1j*y],[x+1j*y,s-z]],dtype=complex))
def exact(b,k):
    J=np.zeros((3,3));J[:2,:2]=b['T'];z=k[2];J[2,2]=-b['t']*np.sin(z)*(1 if b['a'] is None else 2*(np.cos(z)-b['m']))
    return J,np.array([b['u'],0.,-b['v']*np.sin(z)])
def angle(a,b,c):return 2*np.arctan2(np.dot(a,np.cross(b,c)),1+a@b+b@c+c@a)
def winding(b,k,r,N):
    total=0.;small=np.inf;grid=np.linspace(-r,r,N+1)
    for axis in range(3):
        i,j=(axis+1)%3,(axis+2)%3
        for sign in (-1,1):
            pts=np.zeros((N+1,N+1,3))
            for x,y in itertools.product(range(N+1),repeat=2):
                delta=np.zeros(3);delta[axis]=sign*r;delta[i]=grid[x];delta[j]=grid[y];d=dvec(b,k+delta);norm=np.linalg.norm(d);small=min(small,norm);pts[x,y]=d/norm
            for x,y in itertools.product(range(N),repeat=2):
                A,B,C,D=pts[x,y],pts[x+1,y],pts[x+1,y+1],pts[x,y+1];total+=sign*(angle(A,B,C)+angle(A,C,D))
    return float(total/(4*np.pi)),float(small)

def measure(b,k):
    J,w=exact(b,k); B=J.T@J; assert np.linalg.norm(dvec(b,k))<1e-12 and min(np.linalg.eigvalsh(B))>0
    raw=[]
    for h in (.02,.01,.005):
        values=[]
        for u in DIR:
            ep,em=energies(b,k+h*u),energies(b,k-h*u)
            values.append([((ep[1]-ep[0])**2+(em[1]-em[0])**2)/(8*h*h),(ep.mean()-em.mean())/(2*h)])
        raw.append(np.array(values))
    fits=[]
    for a,c in zip(raw[:-1],raw[1:]):
        v=(4*c-a)/3; co=np.linalg.lstsq(X[:-1],v[:-1,0],rcond=None)[0]; tilt=np.linalg.lstsq(DIR[:-1],v[:-1,1],rcond=None)[0]
        fits.append(dict(B=mat(co),w=tilt,holdout_B=abs(X[-1]@co-v[-1,0])/(np.linalg.norm(B)/np.sqrt(3)),holdout_w=abs(DIR[-1]@tilt-v[-1,1])/np.sqrt(np.trace(B)/3)))
    det=np.linalg.det(J); winds=[dict(r=r,N=N,winding=z[0],min_boundary_norm=z[1]) for r in (.02,.04) for N in (6,12) for z in [winding(b,k,r,N)]]
    assert all(abs(z['winding']-np.sign(det))<1e-8 and z['min_boundary_norm']>0 for z in winds)
    fine=fits[-1]; norm=np.linalg.norm(B);speed=np.sqrt(np.trace(B)/3);errB=np.linalg.norm(fine['B']-B)/norm;errw=np.linalg.norm(fine['w']-w)/speed
    uncB=np.linalg.norm(fine['B']-fits[0]['B'])/norm;uncw=np.linalg.norm(fine['w']-fits[0]['w'])/speed
    assert max(errB,errw,uncB,uncw,fine['holdout_B'],fine['holdout_w'])<.002
    tilt_type=w@np.linalg.solve(B,w);assert tilt_type<1
    return dict(k=k,energy=d0(b,k),exact_B=B,exact_w=w,chirality=int(np.sign(det)),winding=winds,raw=raw,fits=fits,errB=errB,errw=errw,uncB=uncB,uncw=uncw,type_I_ratio=tilt_type)

def run(case):
    nodes=[]
    for j,b in enumerate(case['blocks']):
        roots=[b['m']] if b['a'] is None else [b['m']-b['a'],b['m']+b['a']]
        assert all(abs(x)<1 for x in roots) and abs(np.linalg.det(b['T']))>0
        # If either transverse component is pi, Wilson term >=6; bound negative f contribution.
        lower=-1-abs(b['m']) if b['a'] is None else -b['a']**2
        assert b['t']>0 and 6+b['t']*lower>0
        for x in roots:
            for sign in (-1,1):
                node=measure(b,np.array([0.,0.,sign*np.arccos(x)]));node['block']=j;node['root']=x;nodes.append(node)
    assert sum(n['chirality'] for n in nodes)==0
    pairs=[]
    for i,j in itertools.combinations(range(len(nodes)),2):
        a,b=[nodes[x]['fits'][-1] for x in (i,j)]; Bs=[a['B'],b['B']];speed=np.sqrt((np.trace(Bs[0])+np.trace(Bs[1]))/6)
        pairs.append(dict(i=i,j=j,B_mismatch=np.linalg.norm(Bs[0]-Bs[1])/np.mean([np.linalg.norm(z) for z in Bs]),w_mismatch=np.linalg.norm(a['w']-b['w'])/speed,energy_difference=abs(nodes[i]['energy']-nodes[j]['energy'])))
    return dict(case=case,nodes=nodes,pairs=pairs,max_B=max(p['B_mismatch'] for p in pairs),max_w=max(p['w_mismatch'] for p in pairs),max_energy=max(p['energy_difference'] for p in pairs))

def main():
    start=time.process_time();resource.setrlimit(resource.RLIMIT_CPU,(6900,6901));rows=[]
    for c in cases():
        r=run(c);rows.append(r);print(c['arm'],c['name'],r['max_B'],r['max_w'],[n['chirality'] for n in r['nodes']],flush=True)
    out=Path(__file__).resolve().parents[1]/'results/EX-GC30B2';out.mkdir(parents=True,exist_ok=True)
    def conv(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        raise TypeError(type(x))
    cpu=time.process_time()-start
    for i in range(0,len(rows),5):
        (out/f'cases_{i//5}.json').write_text(json.dumps(dict(rows=rows[i:i+5],cpu_s=cpu,python=platform.python_version(),numpy=np.__version__,directions=DIR),default=conv,separators=(',',':'),allow_nan=False)+'\n')
    print('CPU seconds',cpu)

if __name__=='__main__':main()
