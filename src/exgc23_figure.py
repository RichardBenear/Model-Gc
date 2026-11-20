"""Generate the EX-GC23 trajectory/velocity/ledger comparison from compact JSON."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).resolve().parents[1]; out=root/'results/EX-GC23'
fig,axs=plt.subplots(2,2,figsize=(11,8),layout='constrained')
styles={'strong_off':('tab:red','--','Strong, transport off'),'strong_on':('tab:blue','-','Strong, transport on'),'weak_off':('tab:orange','--','Weak, transport off'),'weak_on':('tab:green','-','Weak, transport on'),'wall_on':('tab:purple',':','Strong, larger box'),'strong_cons':('tab:brown','-.','Strong, conservative variant'),'flux_K02':('tab:blue','-','K = 0.2'),'flux_K04':('tab:orange','--','K = 0.4')}
for name,(c,ls,label) in styles.items():
    p=out/(name+'.json')
    if not p.exists():continue
    r=json.loads(p.read_text());a=r['rows'];t=np.array([x['t'] for x in a]);v=lambda k:np.array([x[k] for x in a]); kwargs=dict(color=c,ls=ls,label=label)
    if r['kick']:
        axs[0,0].plot(t,v('Wcen_z')-v('Wcen_z')[0],**kwargs);axs[0,1].plot(t,v('Vcore_z'),**kwargs)
        axs[1,0].plot(t,v('M')/v('M')[0],**kwargs)
    else:axs[1,1].plot(t,v('Wcen_z')-v('Wcen_z')[0],**kwargs)
t=np.linspace(0,20,101);axs[0,0].plot(t,.05*t,'k:',label='ballistic, v=0.05');axs[0,1].axhline(.05,color='k',ls=':');axs[1,1].plot(t,.5*.002*t*t,'k:',label='nominal 0.002 acceleration')
for ax,title,y in zip(axs.flat,['Kick: physical W displacement','Kick: internal flow','Kick: core survival','Uniform flux: W displacement'],['W centroid shift','s-weighted Vz','M/M(0)','W centroid shift']):
    ax.set(title=title,xlabel='Time (model units)',ylabel=y);ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle('EX-GC23 · Material transport of velocity (exploratory)',fontsize=15);fig.savefig(root/'paper/figures/exgc23_transport.png',dpi=160)
