"""Plot and independently recompute compact EX-GC26B boost records."""
from pathlib import Path
import json, hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC26B'
def main():
    cases=[json.loads(p.read_text()) for p in sorted(OUT.glob('h*_L*_u*_axis*_dt*.json'))];fig,axs=plt.subplots(2,2,figsize=(11,8));checks=[]
    for r in cases:
        if not r['u']:continue
        h,L,u,axis=r['h'],r['L'],r['u'],r['axis'];n=np.array(axis)/np.linalg.norm(axis);rows=r['rows'];t=np.array([z['t'] for z in rows]);P=np.array([z['P'] for z in rows]);C=np.array([z['centroid'] for z in rows]);H=np.array([z['H'] for z in rows]);E=r['E0'];v=np.dot(t-t.mean(),C@n)/(np.sum((t-t.mean())**2));g=1/np.sqrt(1-v*v);mass=(P@n).mean()/(g*v*E);rel=H.mean()*np.sqrt(1-u*u)/E;assert abs(mass-r['summary']['M_over_E0'])<1e-11;assert abs(rel-r['summary']['H_over_gamma_E0'])<1e-11
        checks.append(dict(h=h,L=L,u=u,axis=axis,dt=r['dt'],speed=v,M_over_E0=mass,H_over_gamma_E0=rel))
    for h in sorted(set(x['h'] for x in checks)):
        for axis in ([1,0,0],[1,1,1]):
            subset=sorted([x for x in checks if x['h']==h and x['axis']==axis and x['L']==7 and abs(x['dt']-.05*h)<1e-12],key=lambda x:x['u']);u=[x['u'] for x in subset];lab=f'h={h}, '+('axis' if axis[1]==0 else 'diagonal')
            for ax,key in zip(axs.flat[:3],['M_over_E0','H_over_gamma_E0','speed']):
                y=[100*(x[key]-1) if key!='speed' else 100*(x[key]/x['u']-1) for x in subset];ax.plot(u,y,'o-',label=lab)
    for ax,title in zip(axs.flat[:3],['Dynamic mass error (%)','Boost energy error (%)','Measured speed error (%)']):ax.axhspan(-1,1,color='green',alpha=.09);ax.axhline(0,color='k',lw=.7);ax.set(xlabel='Imposed speed / c',ylabel=title);ax.grid(alpha=.2)
    axs[0,0].legend(fontsize=8)
    ax=axs[1,1]
    for h in sorted(set(r['h'] for r in cases)):
        rs=[r for r in cases if r['h']==h and r['L']==7 and r['axis']==[1,0,0] and r['u'] and abs(r['dt']-.05*h)<1e-12];rs.sort(key=lambda r:r['u']);u=[r['u'] for r in rs];ax.semilogy(u,[max(r['summary']['regions'][1]['raw_max'],1e-14) for r in rs],'o--',label=f'h={h} raw');ax.semilogy(u,[max(r['summary']['regions'][1]['residual_max'],1e-14) for r in rs],'s-',label=f'h={h} flux corrected')
    ax.axhline(1e-3,color='k',lw=.7);ax.set(xlabel='Imposed speed / c',ylabel='Regional momentum residual / E₀');ax.grid(alpha=.2);ax.legend(fontsize=8);fig.suptitle('EX-GC26B — numerical Lorentz mass calibration');fig.tight_layout();fig.savefig(ROOT/'paper/figures/exgc26b_mass.png',dpi=160);plt.close(fig);(OUT/'recompute.json').write_text(json.dumps(checks,indent=1)+'\n')
    raw=ROOT/'local_raw/EX-GC26B';manifest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(raw.glob('*.npy'))];(OUT/'raw_manifest.json').write_text(json.dumps(manifest,indent=1)+'\n')
if __name__=='__main__':main()
