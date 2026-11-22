"""Recompute GC27 source differences and virial separation from checkpoints."""
import json
from pathlib import Path
import numpy as np

OUT=Path(__file__).resolve().parents[1]/'results/EX-GC27AB'
def load(name):return json.loads((OUT/(name+'.json')).read_text())
def main():
    result={}
    for phase in ('A','B'):
        result[phase]={}
        for h in (.2,.15):
            states=[]
            for B in (1,2):
                path=OUT/f'{phase}_sensitivity_h{h:g}_B{B}.json'
                if not path.exists():break
                states.append(json.loads(path.read_text()))
            if len(states)!=2:continue
            E=[x['base']['E'] for x in states];v=[x['base']['virial'] for x in states];p=[x['base']['parts'][2]/x['base']['E'] for x in states];bind=2*E[0]-E[1];rows=[]
            for a,b in zip(states[0]['cases'],states[1]['cases']):
                assert a['exponents']==b['exponents'];x,y,z=a['exponents'];k=[a['k'],b['k']];sb=2*a['s']-b['s'];vir=[(x-y)*q/2 for q in v];k0=[k[i]-vir[i] for i in range(2)];expected=[(x+y)/2+(z-2*x+y)*q for q in p]
                assert max(abs(k0[i]-expected[i]) for i in range(2))<1e-12
                rows.append(dict(exponents=[x,y,z],k=k,spread=k[1]-k[0],virial_component=vir,k_without_measured_virial=k0,spread_without_virial=k0[1]-k0[0],binding_source=sb,binding_sensitivity=sb/bind,binding_residual_vs_B1=sb-k[0]*bind,fixed_fd_error=[a['k_fd_fixed']-k[0],b['k_fd_fixed']-k[1]],relaxed_fd=[a.get('k_fd_relaxed'),b.get('k_fd_relaxed')],dilation_fd=[a.get('k_fd_dilation'),b.get('k_fd_dilation')]))
            result[phase][str(h)]=dict(E=E,binding=bind,binding_fraction=bind/(2*E[0]),virial=v,potential_fraction=p,radius=[x['base']['R'] for x in states],max_gradient=[x['base']['gradient'] for x in states],cases=rows)
    if (OUT/'completion.json').exists():result['completion']=load('completion')
    (OUT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
