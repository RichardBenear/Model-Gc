import sys,os,json;sys.path.insert(0,'/home/user/model-gc/src');import numpy as np
import exgc05_rotation as r5, exgc02_inertia as g2
M,u=g2.load('ref'); rows=[]
for tW in [float(x) for x in sys.argv[1:]]:
    for l in (0,2,3):
        sp_=r5.spectrum(M,u,l,{"W":tW}); loc=[m for m in sp_ if m['core']>0.5]
        rows.append(dict(tauW=tW,l=l,modes=sp_[:6],core_modes=loc[:6]))
        print('tauW',tW,'l',l,' | '.join(f"{m['re']:+.4f}{m['im']:+.3f}i c{m['core']:.2f}" for m in loc[:4]),flush=True)
json.dump(rows,open(f'results/EX-GC05/scan_tauW_{sys.argv[1]}.json','w'),indent=1)
