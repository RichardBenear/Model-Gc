"""Post-hoc analytic energy-current diagnostic; no tuned criteria."""
import itertools,json,time
import numpy as np
import exgc26b2 as m
import exgc26b_mass as b
import exgc26a_skyrmion as ref

def current(G,p,v,orientations=8):
    # H_i = v_i pi_i - L_i, dH_i/dt = sum_j(v_i dL_j/dphi_i - v_j dL_i/dphi_j).
    # For fixed walls, sum_i x_i dH_i/dt = sum_i,j (x_j-x_i) v_j dL_i/dphi_j.
    # Neighbor derivative (s/h)L_di cancels displacement s*h, leaving L_di dot v_neighbor.
    # The canonical continuum P estimator instead uses v at the originating node.
    J=np.zeros(3);P=np.zeros(3);vv=b.dot(v,v)
    signs=list(itertools.product((-1,1),repeat=3)) if orientations==8 else [(-1,-1,-1),(1,1,1)]
    for s in signs:
        mask,d,D=m.orient(G,p,s);tr=sum(D[i,i] for i in range(3));vd=[b.dot(v,x) for x in d]
        for i in range(3):
            Ld=2*((vv-1-tr)*d[i]-vd[i]*v+sum(D[i,j]*d[j] for j in range(3)))
            J[i]+=np.sum(b.dot(Ld,np.roll(v,-s[i],axis=i+1))*mask)/orientations
            P[i]+=np.sum(b.dot(Ld,v)*mask)/orientations
    fac=G.h**3/ref.NORM;return J*fac,P*fac
def main():
    start=time.process_time();m.load();G=ref.Grid(.25,1);p=G.hedgehog(ref.PROF);v=ref.tangent(p,np.random.default_rng(2622).normal(scale=.05,size=p.shape))*G.inner;a=m.accel(G,p,v);j,pc=current(G,p,v);eps=1e-6
    def moment(p,v):H=m.fields(G,p,v)[0];return np.array([np.sum(x*H) for x in G.X])*G.h**3/ref.NORM
    qp=G.fix(p+eps*v);qm=G.fix(p-eps*v);fd=(moment(qp,ref.tangent(qp,v+eps*a))-moment(qm,ref.tangent(qm,v-eps*a)))/(2*eps);err=float(abs(j-fd).max());assert err<1e-8
    out=[]
    for h,L,u,axis in ((.2,7,.1,[1,0,0]),(.2,7,.3,[1,0,0]),(.2,7,.3,[1,1,1]),(.15,7,.1,[1,0,0]),(.15,7,.3,[1,0,0])):
        G=ref.Grid(h,L);p=np.load(m.restpath(h,L));n=np.array(axis)/np.linalg.norm(axis);p,v=b.init_boost(G,p,u,n);H,P,_=m.fields(G,p,v);energy=H.sum()*G.h**3/ref.NORM;j,pc=current(G,p,v);assert abs(pc-P.sum(axis=(1,2,3))*G.h**3/ref.NORM).max()<1e-10;out.append(dict(h=h,L=L,u=u,axis=axis,P_over_H=float(pc@n/energy),J_over_H=float(j@n/energy),J_over_P=float(j@n/(pc@n)),J=j.tolist(),P=pc.tolist()))
    r=dict(derivative_validation_max=err,rows=out,cpu_s=time.process_time()-start);m.save('discrete_current',r);print(json.dumps(r,indent=1),flush=True)
def inherited(source):
    from pathlib import Path
    start=time.process_time();m.load();out=[]
    for h in (.2,.15):
        G=ref.Grid(h,7);rest=np.load(Path(source)/f'rest_h{h:g}_L7.npy')
        for u in (.1,.3,.5):
            for axis in ([1,0,0],[1,1,1]):
                n=np.array(axis)/np.linalg.norm(axis);p,v=b.init_boost(G,rest,u,n);j,pc=current(G,p,v,2);rec=json.loads((m.ROOT/'results/EX-GC26B'/f'h{h:g}_L7_u{u:g}_axis{"".join(map(str,axis))}_dt0.05.json').read_text());measured=rec['summary']['speed']/np.mean([np.array(row['P'])@n/row['H'] for row in rec['rows']]);out.append(dict(h=h,u=u,axis=axis,J_over_P=float(j@n/(pc@n)),measured_ratio=float(measured)))
    m.save('inherited_current',dict(rows=out,cpu_s=time.process_time()-start));print(out,flush=True)
if __name__=='__main__':
    import sys
    if len(sys.argv)>1:inherited(sys.argv[1])
    else:main()
