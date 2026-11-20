"""Recompute EX-GC23 diagnostics from compact records (no solver dependency)."""
import json, sys
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'results/EX-GC23'

def summarize(r):
    rows=r['rows'];t=np.array([x['t'] for x in rows]);a={key:np.array([x[key] for x in rows]) for key in ('M','Wmax','Wcen_z','scen_z','Vcore_z','Vspread','depletion')}
    out={k:r[k] for k in ('name','K','L','mode','kick','flux','status','cpu_s')};out.update(T=float(t[-1]),M_ratio=float(a['M'][-1]/a['M'][0]),Wmax_ratio=float(a['Wmax'][-1]/a['Wmax'][0]),
        M_max_change=float(abs(a['M']/a['M'][0]-1).max()),Wmax_max_change=float(abs(a['Wmax']/a['Wmax'][0]-1).max()),displacement=float(a['Wcen_z'][-1]-a['Wcen_z'][0]),s_displacement=float(a['scen_z'][-1]-a['scen_z'][0]),V_initial=float(a['Vcore_z'][0]),V_final=float(a['Vcore_z'][-1]),
        max_ledger_residual=float(max(abs(np.array(x['residual'])).max() for x in rows)),ledger_z=np.array(rows[-1]['ledger'])[:,2].tolist(),P_z_change=float(rows[-1]['P'][2]-rows[0]['P'][2]))
    for lo,hi,label in ((2,10,'early'),(10,20,'late')):
        ix=(t>=lo-.0002)&(t<=hi+.0002)
        if ix.sum()<4 or t[-1]<hi-.001:continue
        for key in ('Wcen_z','scen_z','Vcore_z'):
            out[f'{key}_slope_{label}']=float(np.polyfit(t[ix],a[key][ix],1)[0])
        out[f'V_mean_{label}']=float(a['Vcore_z'][ix].mean());out[f'W_speed_over_V_{label}']=float(out[f'Wcen_z_slope_{label}']/out[f'V_mean_{label}']) if abs(out[f'V_mean_{label}'])>1e-10 else None
        if min(a['Vcore_z'][ix])>0:out[f'gamma_{label}']=float(-np.polyfit(t[ix],np.log(a['Vcore_z'][ix]),1)[0])
    if r['kick']:
        out['V_retention']=float(a['Vcore_z'][-1]/a['Vcore_z'][0]);out['displacement_over_ballistic']=float(out['displacement']/(r['kick']*t[-1])) if t[-1]>0 else None
    if r['flux']:
        for end in (5,10):
            ix=t<=end+.0002
            if t[-1]<end-.001:continue
            out[f'flow_accel_{end}']=float(np.polyfit(t[ix],a['Vcore_z'][ix],1)[0]);out[f'centroid_accel_{end}']=float(2*np.polyfit(t[ix],a['Wcen_z'][ix],2)[0]);out[f's_centroid_accel_{end}']=float(2*np.polyfit(t[ix],a['scen_z'][ix],2)[0])
    return out

def main():
    rows=[];settle=[]
    for f in sorted(OUT.glob('*.json')):
        r=json.loads(f.read_text())
        if 'name' in r:rows.append(summarize(r))
        elif f.stem.startswith('settle_'):
            a=r['rows'];settle.append(dict(K=r['K'],L=r['L'],cpu_s=r['cpu_s'],M=a[-1]['M'],Wmax=a[-1]['Wmax'],M_change_last5=a[-1]['M']/a[-6]['M']-1,sha256=r['sha256']))
    result=dict(cases=rows,settle=settle,total_cpu_s=sum(r['cpu_s'] for r in rows+settle))
    (OUT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    for r in rows:print(r['name'],r['status'],'W displacement',round(r['displacement'],5),'V final',round(r['V_final'],6),'M ratio',round(r['M_ratio'],5),'late W speed',r.get('Wcen_z_slope_late'),'ledger',r['max_ledger_residual'])
    print('CPU seconds',result['total_cpu_s'])

if __name__=='__main__':main()
