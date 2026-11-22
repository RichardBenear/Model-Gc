"""Prospectively registered EX-GC29A phase-dispersion diagnostic."""
from exgc29a import *

def main():
    start=time.process_time(); resource.setrlimit(resource.RLIMIT_CPU,(6900,6901)); rows=[]
    R=24.; h0=.0125; r0=np.arange(1,round(R/h0))*h0
    for C in (0.,C0):
        base=canonical(h0,R,C,0.)
        fs=CubicSpline(np.r_[0,r0,R],np.r_[np.pi,base['profile'],0],bc_type='natural')
        ds=CubicSpline(r0,base['delta_potential'],bc_type=((1,0.),(2,0.)))
        for h in (h0,.00625,.003125):
            r=np.arange(1,round(R/h))*h; f=fs(r); K=r*r+2*np.sin(f)**2
            dv=2/h**2+2/r**2+C/2; off=np.full(len(r)-1,-1/h**2)
            for system,d in [('defect',dv+ds(r)),('vacuum',dv)]:
                ev,V=eigh_tridiagonal(d,off,select='v',select_range=(0.,64**2),lapack_driver='stebz'); om=np.sqrt(ev)
                for op,p in [('P1',np.sin(f)),('P2',2*(2*np.sin(f/2)**2)*np.sin(f))]:
                    J=r*r*p/np.sqrt(6*np.pi*K)
                    for drive,k in [('continuum',om),('matched',2*np.arcsin(h*om/2)/h)]:
                        w=np.empty_like(om)
                        for i in range(0,len(om),32):
                            s=slice(i,i+32); w[s]=h*np.sum(V[:,s]*J[:,None]*np.sinc(r[:,None]*k[None,s]/np.pi),axis=0)**2
                        cuts=[dict(cut=c,moments=[float(np.sum(w[om<=c]*om[om<=c]**pwr)) for pwr in (0,2,4)]) for c in (8,16,32,64)]
                        rows.append(dict(C=C,h=h,system=system,operator=op,drive=drive,cuts=cuts))
                print(C,h,system,len(om),flush=True)
    (OUT/'dispersion.json').write_text(json.dumps(dict(cpu_s=time.process_time()-start,rows=rows),separators=(',',':'))+'\n')

if __name__=='__main__':main()
