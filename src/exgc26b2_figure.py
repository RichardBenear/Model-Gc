"""Independent compact-record recomputation, plot, retention manifest."""
from pathlib import Path
import hashlib,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC26B2'
def main():
    cases=[json.loads(p.read_text()) for p in sorted(OUT.glob('h*.json'))];checks=[];fig,axs=plt.subplots(2,2,figsize=(11,8));cur=json.loads((OUT/'discrete_current.json').read_text())['rows'];ax=axs[0,0]
    for r in cases:
        t=np.array([s['t'] for s in r['rows']]);n=np.array(r['axis'])/np.linalg.norm(r['axis']);C=np.array([s['centroid'] for s in r['rows']])@n;P=np.array([s['P'] for s in r['rows']])@n;H=np.array([s['H'] for s in r['rows']]);F=np.array([s['F'] for s in r['rows']]);v=np.dot(t-t.mean(),C)/np.sum((t-t.mean())**2);record=dict(h=r['h'],u=r['u'],T=r['T'],direction=r['axis'],force_sign=r['force_sign'],force_power=r.get('force_power',1),speed=v)
        if r['u']:
            ratio=v/np.mean(P/H);mass=P.mean()*np.sqrt(1-v*v)/(v*r['E0']);assert abs(ratio-r['summary']['speed_over_P_H'])<1e-10 and abs(mass-r['summary']['M_over_E0'])<1e-10;record.update(ratio=ratio,mass=mass)
            if r['T']<1:
                q=next(x for x in cur if x['h']==r['h'] and x['u']==r['u'] and x['axis']==r['axis']);ax.scatter(q['J_over_P'],ratio,s=50,label=f"h={r['h']}, u={r['u']}, {''.join(map(str,r['axis']))}")
            else:
                ax2=axs[0,1];ax2.plot(t,P/H,label='canonical P/H');ax2.plot((t[1:]+t[:-1])/2,np.diff(C)/np.diff(t),label='centroid speed');ax2.axhline(r['u'],color='k',ls=':',label='launch');ax2.set(xlabel='Time',ylabel='Speed / c',title='Long free coast');ax2.legend(fontsize=8)
        elif r['force_sign']:
            tc=t-t.mean();X=np.array([np.ones_like(tc),tc,tc**2]).T;coef=np.linalg.lstsq(X,C,rcond=None)[0];a=2*coef[2];mass=F.mean()/(a*r['E0']);assert abs(mass-r['summary']['mass_over_E0_force'])<1e-10;imp=np.r_[0.,np.cumsum((F[:-1]+F[1:])*np.diff(t)/2)];record.update(acceleration=a,mass=mass,impulse_error=float(max(abs(P-P[0]-imp))/abs(imp[-1])));lab=f"force {r['force_sign']:+d}, power {r.get('force_power',1)}";axs[1,0].plot(t,C-C[0],label=lab);axs[1,1].plot(t,(P-P[0])/r['E0'],label=lab);axs[1,1].plot(t,imp/r['E0'],'--',alpha=.5)
        else:axs[1,0].plot(t,C-C[0],label='rest control')
        checks.append(record)
    ax.plot([.98,1.002],[.98,1.002],':k');ax.axhspan(.995,1.005,color='green',alpha=.1);ax.set(xlabel='Initial discrete J/P',ylabel='Measured v / (P/H)',title='Exact-current diagnostic');ax.legend(fontsize=7)
    axs[1,0].set(xlabel='Time',ylabel='Centroid displacement',title='Signed force response');axs[1,0].legend(fontsize=8);axs[1,1].set(xlabel='Time',ylabel='Momentum change / E₀',title='Momentum (solid) vs applied impulse (dash)');axs[1,1].legend(fontsize=8)
    for a in axs.flat:a.grid(alpha=.25)
    fig.suptitle('EX-GC26B2 — mirror symmetry is not continuum Lorentz symmetry');fig.tight_layout();fig.savefig(ROOT/'paper/figures/exgc26b2_qualification.png',dpi=160);plt.close(fig);(OUT/'recompute.json').write_text(json.dumps(checks,indent=1)+'\n');raw=ROOT/'local_raw/EX-GC26B2';manifest=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(raw.glob('*.npy'))];(OUT/'raw_manifest.json').write_text(json.dumps(manifest,indent=1)+'\n')
if __name__=='__main__':main()
