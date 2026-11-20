"""Compact independent analysis from saved EX-GC22 time-series (no solver calls)."""
import json
from pathlib import Path
import numpy as np

root=Path(__file__).resolve().parents[1]; p=root/'results/EX-GC22'
f=json.loads((p/'fixed.json').read_text()); rows=[]
for file in sorted(p.glob('full_*.json')):
    j=json.loads(file.read_text()); s=j['settle']; a=j['kick']; last=a[-1]; t=np.array([r['t'] for r in a]); v=np.array([r['Vcore_z'] for r in a]); z=np.array([r['Wcen_z'] for r in a]); late=(t>=10)&(t<=20)
    gamma=-np.polyfit(t[late],np.log(v[late]),1)[0] if np.all(v[late]>0) else None
    row=dict(K=j['K'],D=j['D'],h=j['h'],M0=a[0]['M'],depletion0=a[0]['depletion'],settle_M_rel=s[-1]['M']/s[-6]['M']-1,
             gamma_late=gamma,gamma_early=j['gamma_early'],velocity_retention=last['Vcore_z']/a[0]['Vcore_z'],M_retention=last['M']/a[0]['M'],Wmax_retention=last['Wmax']/a[0]['Wmax'],
             centroid_z=last['Wcen_z'],centroid_speed_early=np.polyfit(t[t<=5],z[t<=5],1)[0],centroid_speed_late=np.polyfit(t[late],z[late],1)[0],flow_speed_late=v[late].mean(),cpu_s=j['cpu_s'])
    row['centroid_flow_ratio']=row['centroid_speed_late']/row['flow_speed_late']
    if 'flux' in j:
        b=j['flux']; tt=np.array([r['t'] for r in b]); vv=np.array([r['Vcore_z'] for r in b]); zz=np.array([r['Wcen_z'] for r in b]); early=tt<=5
        row['flux_accel_early']=np.polyfit(tt[early],vv[early],1)[0]
        row['flux_centroid_accel_early']=2*np.polyfit(tt[early],zz[early],2)[0]
        row['null_max_abs_velocity']=max(abs(r['Vcore_z']) for r in j['null'])
        row['null_max_abs_centroid']=max(abs(r['Wcen_z']) for r in j['null'])
    rows.append(row)
summary=dict(full=rows,fixed_cpu_s=f['cpu_s'],full_cpu_s=sum(r['cpu_s'] for r in rows),fixed_cases=len(f['rows']),fixed_controls=len(f['controls']),
             max_fixed_residual=max(r['residual'] for r in f['rows']+f['controls']))
(p/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
for r in rows: print(' '.join(f'{k}={v:.6g}' for k,v in r.items() if v is not None))
