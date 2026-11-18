import sys,os,json;sys.path.insert(0,'/home/user/model-gc/src');import numpy as np
os.chdir('/home/user/model-gc')
import exgc04_boundpair as b
pr,pu=b._prof()
d0,d1,step=[float(x) for x in sys.argv[1:4]]
a=b.AxiN([-d0/2,d0/2]); u,v,res,ok,it=a.solve(a.guess(pr,pu,b.M_P1),maxit=30)
rows=[dict(d=d0,v_right=float(v[1]),ok=bool(ok),res=res,Mc=a.consumption(u,1),Wmax=float(u[2].max()))]
print(rows[-1],flush=True)
d=d0
while (d1-d)*np.sign(step)>1e-9:
    dn=round(d+step,4); an=b.AxiN([-dn/2,dn/2])
    zo=a.Z[0]; zn=an.Z[0]; src=np.where(zn>0,zn-(dn-d)/2,zn+(dn-d)/2)
    g=np.array([[np.interp(src,zo,u[s].reshape(a.nr,a.nz)[i]) for i in range(a.nr)] for s in range(4)])
    # rho grids are identical (same rho_core, h, far)
    un,v,res,ok,it=an.solve(g.reshape(4,-1),maxit=30)
    rows.append(dict(d=dn,v_right=float(v[1]),ok=bool(ok),res=res,it=it,Mc=an.consumption(un,1),Wmax=float(un[2].max())))
    print(rows[-1],flush=True)
    if not ok: break
    a,u,d=an,un,dn
json.dump(rows,open(f'results/EX-GC04/vd_cont_{d0:g}_{d1:g}.json','w'),indent=1)
