"""Recompute B3 paired-control force error budget from recorded trajectories."""
import json
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results/EX-GC26B3'

def unpack(r):
    rows=r['rows'];t=np.array([q['t'] for q in rows]);P=np.array([q['P'] for q in rows]);C=np.array([q['centroid'] for q in rows]);H=np.array([q['H'] for q in rows]);V=np.array([q['V'] for q in rows]);F=np.array([q['F'] for q in rows]);D=np.array([q['D'] for q in rows]);dr=np.where(t<1,.5*np.pi*np.sin(np.pi*t),0.)*r['f'];W=cumulative_trapezoid(dr*D,t,initial=0);I=cumulative_trapezoid(F,t,initial=0)
    R=np.array([[x['centroid'][0] for x in q['regions']] for q in rows]);return t,P,C,H,V,I,W,R

def pairs(f,n):
    t,P,C,H,V,I,W,R=unpack(f);tn,Pn,Cn,Hn,Vn,In,Wn,Rn=unpack(n);assert np.array_equal(t,tn);E=f['E0'];dP=P[:,0]-P[0,0];dPn=Pn[:,0]-Pn[0,0];rawC=C[:,0];c=rawC-Cn[:,0];p=dP-dPn;late=t>=4;speed=lambda y:float(np.polyfit(t[late],y[late],1)[0]);v=speed(c);vr=speed(rawC);pm=float(p[late].mean());pr=float(dP[late].mean());fits={}
    for start in (1,2):
        w=t>=start;acc=lambda y:float(2*np.polyfit(t[w],y[w],2)[0]);accel=acc(c);fits[str(start)]=dict(acceleration=accel,momentum_slope=float(np.polyfit(t[w],p[w],1)[0]),force_over_accel_E0=float(np.polyfit(t[w],I[w],1)[0]/accel/E),regional_acceleration=[acc(R[:,k]-Rn[:,k]) for k in range(2)])
    return dict(h=f['h'],E0=E,T=f['T'],impulse_E0=float(I[-1]/E),momentum_gain_E0=float(p[-1]/E),impulse_transfer=float(p[-1]/I[-1]),raw_impulse_transfer=float(dP[-1]/I[-1]),max_impulse_residual_E0=float(max(abs(p-I))/E),speed=v,mean_P_E0=pm/E,mass_P_over_v_E0=pm/v/E,mass_P_over_gamma_v_E0=pm/v/E*np.sqrt(1-v*v),raw_mass_P_over_v_E0=pr/vr/E,regional_late_speed=[speed(R[:,k]-Rn[:,k]) for k in range(2)],fits=fits,control_speed=speed(Cn[:,0]),control_max_P_E0=float(np.max(np.linalg.norm(Pn-Pn[0],axis=1))/E),max_transverse_P_E0=float(np.max(np.linalg.norm(P[:,1:],axis=1))/E),bare_E_change=float(H[-1]/H[0]-1),ledger_max_rel=float(max(abs(H+V+W-H[0]))/E),null_energy_drift=float(max(abs(Hn/Hn[0]-1))),degree=[f['degree_initial'],f['degree_final'],n['degree_initial'],n['degree_final']],cpu_s=f['cpu_s']+n['cpu_s'])

def main():
    out=[]
    for h in (.2,.15):
        f,n=[json.loads((OUT/f'h{h:g}_force{s}.json').read_text()) for s in (1,0)];out.append(pairs(f,n))
    result=dict(cases=out,core_hours=sum(r['cpu_s'] for r in out)/3600);(OUT/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
