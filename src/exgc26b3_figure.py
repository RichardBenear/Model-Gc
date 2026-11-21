"""Compact comparative error-budget figure from B3 saved records."""
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from exgc26b3_summary import OUT,ROOT,unpack

fig,ax=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
for h,color in ((.2,'tab:blue'),(.15,'tab:orange')):
    f,n=[json.loads((OUT/f'h{h:g}_force{s}.json').read_text()) for s in (1,0)]
    t,P,C,H,V,I,W,R=unpack(f);_,Pn,Cn,*_=unpack(n);E=f['E0'];p=P[:,0]-Pn[:,0];p-=p[0];x=C[:,0]-Cn[:,0];x-=x[0];label=f'h={h}'
    ax[0].plot(t,I/E,'--',color=color,alpha=.65);ax[0].plot(t,p/E,color=color,label=label)
    ax[1].plot(t,(p-I)/E,color=color,label=label)
    ax[2].plot(t,x,color=color,label=label)
for a in ax:
    a.axvspan(0,1,color='gray',alpha=.13);a.set_xlabel('time');a.grid(alpha=.2);a.legend()
ax[0].set_title('Applied impulse (dashed) vs momentum');ax[0].set_ylabel('normalized by E₀')
ax[1].set_title('Impulse transfer residual');ax[1].set_ylabel('(ΔP − I)/E₀')
ax[2].set_title('Control-corrected centroid trajectory');ax[2].set_ylabel('x displacement')
fig.savefig(ROOT/'paper/figures/exgc26b3_calibration.png',dpi=180)
fig.savefig(ROOT/'paper/figures/exgc26b3_calibration.svg')
