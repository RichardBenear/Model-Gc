"""Plot compact decision exports, not raw field data."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];r=ROOT/'results/EX-GC26A2'
def read(p):return json.loads(p.read_text())
st=[read(p) for p in r.glob('static*hedgehog.json')];boost=[read(p) for p in r.glob('boost*.json')];pin=[read(p) for p in r.glob('compact_pin*.json')]
fig,ax=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
for mode,lab in ((0,'original'),(1,'higher order')):
    s=sorted([a for a in st if a['mode']==mode and a['L']==5],key=lambda x:x['h'])
    if s:
        ax[0,0].plot([x['h'] for x in s],[abs(x['B_fd']) for x in s],'o-',label=lab+' finite difference');ax[0,1].plot([x['h'] for x in s],[x['ratio'] for x in s],'o-',label=lab)
        ax[0,0].plot([x['h'] for x in s],[abs(x['geometric']['degrees'][0]) for x in s],'x--',label=lab+' geometric')
for a in ax[0]:a.axhline(1,color='k',ls=':',lw=1);a.set_xlabel('Grid spacing h');a.legend(fontsize=8)
ax[0,0].set_ylabel('|Degree|');ax[0,1].set_ylabel('Derrick ratio E2/E4');ax[0,1].axhspan(.98,1.02,color='green',alpha=.1)
for m,l in ((0,'original'),(1,'higher order')):
    p=sorted([x for x in pin if x['mode']==m and x['L']==5],key=lambda x:x['h'])
    if p:
        ax[1,0].semilogy([x['h'] for x in p],[x['relative_range_fourier'] for x in p],'o-',label=l+' Fourier');ax[1,0].semilogy([x['h'] for x in p],[x['relative_range_cubic'] for x in p],'s--',label=l+' cubic')
ax[1,0].axhline(1e-4,color='k',ls=':');ax[1,0].set_ylabel('Compact translation energy range / E');ax[1,0].set_xlabel('Grid spacing h');ax[1,0].legend(fontsize=8)
for q,lab in (([1,0,0],'axis'),([1,1,1],'diagonal')):
    b=sorted([x for x in boost if x['axis']==q and x.get('L',5)==5 and abs(x['dt']-.05*x['h'])<1e-9],key=lambda x:x['h'])
    if b:ax[1,1].plot([x['h'] for x in b],[x['max_dP_over_E0'] for x in b],'o-',label=lab+' L=5')
for x in boost:
    if x.get('L',5)!=5:ax[1,1].scatter(x['h'],x['max_dP_over_E0'],marker='s',label='axis L=7')
ax[1,1].axhline(1e-3,color='k',ls=':',label='A2 threshold');ax[1,1].set_xlabel('Grid spacing h');ax[1,1].set_ylabel('Max ΔP / E0 (T=1.5)');ax[1,1].legend(fontsize=8)
fig.suptitle('EX-GC26A2 — Numerical qualification (exploratory; no inertia measurement)');fig.savefig(ROOT/'paper/figures/exgc26a2_qualification.png',dpi=160)
