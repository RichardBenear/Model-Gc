"""Replot EX-GC22 from compact committed records; no simulation dependency."""
import json
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

root=Path(__file__).resolve().parents[1]
p=root/'results/EX-GC22'
fixed=json.loads((p/'fixed.json').read_text())
full=[json.loads(f.read_text()) for f in sorted(p.glob('full_*.json'))]
fig,axs=plt.subplots(1,3,figsize=(13,4.1),layout='constrained')
for D,col in [(1,'#b9552d'),(8,'#657b9e'),(80,'#147a76')]:
    rows=sorted([r for r in fixed['rows'] if r['D']==D and r['v']==.05],key=lambda r:r['K'])
    axs[0].loglog([r['K'] for r in rows],[r['gamma'] for r in rows],'o-',c=col,label=f'D={D}')
    axs[0].loglog([r['K'] for r in rows],[r['gamma_G1'] for r in rows],':',c=col)
    axs[1].semilogx([r['K'] for r in rows],[r['external_a']/.002 for r in rows],'o-',c=col,label=f'D={D}')
axs[0].set(xlabel='Sink strength K',ylabel='Velocity decay rate',title='A. Prescribed sink (steady)')
axs[0].legend(title='Dotted: G=1 in force'); axs[0].grid(alpha=.2,which='both')
axs[1].axhline(1,color='k',ls=':',lw=1); axs[1].set(xlabel='Sink strength K',ylabel='External acceleration / 0.002',title='B. Flux response (assigned inertia)'); axs[1].grid(alpha=.2)
for D,col in [(1,'#b9552d'),(80,'#147a76')]:
    cases=sorted([r for r in full if r['D']==D and r['h']==.5],key=lambda r:r['K'])
    if cases:
        axs[2].semilogx([r['K'] for r in cases],[r['gamma_late'] for r in cases],'o-',c=col,label=f'D={D}')
axs[2].axhline(0,color='k',ls=':',lw=1)
axs[2].set(xlabel='Chemical sink strength K',ylabel='Late fitted decay rate',title='C. Free chemical particle (h=0.5)')
if full: axs[2].legend()
axs[2].grid(alpha=.2)
fig.suptitle('EX-GC22 — self-wake drag and weak sourcing (exploratory)')
fig.savefig(root/'paper/figures/exgc22_drag.png',dpi=180)
