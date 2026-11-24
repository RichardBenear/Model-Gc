"""GC35I: independent local proper-time events at a prescribed lapse step."""
import hashlib
import json
from pathlib import Path
import time
import numpy as np
from numba import njit

KY=.60
KX=1.55
TEND=36.
INTERFACE=40.
LENGTH=80.


@njit
def permitted(i,done,period,limit):
    if done[i]>=limit[i]:
        return False
    if i==0 or i==len(period)-1:
        return True
    tick=done[i]*period[i]
    return done[i-1]>=tick//period[i-1] and done[i+1]>=tick//period[i+1]


@njit
def predict(history,i,tick,period):
    step=tick//period[i]
    frac=(tick%period[i])/period[i]
    new=history[i,step+1]
    if frac==0.:
        return new
    return new+frac*(new-history[i,step])


@njit
def update(history,i,done,period,lapse,face,h,proper,ky2):
    if i==0 or i==len(period)-1:
        return 0j
    version=done[i]
    tick=version*period[i]
    u=history[i,version+1]
    left=predict(history,i-1,tick,period)
    right=predict(history,i+1,tick,period)
    flux=(face[i]*(right-u)-face[i-1]*(u-left))/(h*h*lapse[i])
    return 2*u-history[i,version]+proper*proper*(flux-ky2*u)


@njit
def simulate(initial,period,lapse,face,h,proper,ky2,end,seed):
    n=len(period)
    limit=end//period
    history=np.zeros((n,np.max(limit)+2),np.complex128)
    history[:,:2]=initial
    done=np.zeros(n,np.int64)
    total=np.sum(limit)
    witness=np.full(4,-1,np.int64)
    if seed<0:
        for tick in range(1,end+1):
            for i in range(n):
                if tick%period[i]==0:
                    history[i,done[i]+2]=update(history,i,done,period,lapse,face,h,proper,ky2)
                    done[i]+=1
        return history,witness,total
    np.random.seed(seed)
    queued=np.ones(n,np.uint8)
    ready=np.arange(n)
    count=n
    produced=0
    while count>0:
        j=np.random.randint(count)
        i=ready[j]
        count-=1
        ready[j]=ready[count]
        queued[i]=0
        history[i,done[i]+2]=update(history,i,done,period,lapse,face,h,proper,ky2)
        done[i]+=1
        produced+=1
        if witness[0]<0 and 0<i<n-1 and produced%128==0:
            now=done[i]*period[i]
            other=np.argmin(done[1:-1]*period[1:-1])+1
            prior=done[other]*period[other]
            if now-prior>=12:
                witness[:]=i,now,other,prior
        for p in range(max(i-1,0),min(i+2,n)):
            if queued[p]==0 and permitted(p,done,period,limit):
                ready[count]=p
                count+=1
                queued[p]=1
    if produced!=total:
        raise RuntimeError("incomplete schedule")
    return history,witness,produced


def phase_slope(t,z):
    mask=abs(z)>.45*np.max(abs(z))
    tp=t[mask]
    weight=abs(z[mask])**2
    phase=np.unwrap(np.angle(z))[mask]
    center=np.average(tp,weights=weight)
    return float(-np.sum(weight*(tp-center)*phase)/np.sum(weight*(tp-center)**2))


def wave_angle(x,u,sitecut,h):
    z=u[sitecut]
    k=float(np.angle(np.sum(z[:-1].conj()*z[1:]))/h)
    return float(np.degrees(np.arctan2(KY,k))),k


def run(h,proper,step):
    x=np.linspace(0,LENGTH,round(LENGTH/h)+1)
    lapse=np.where((x>=INTERFACE)&step,1.5,1.)
    period=np.where(lapse>1.,2,3).astype(np.int64)
    quantum=proper/3
    end=round(TEND/quantum)
    face=2/(1/lapse[:-1]+1/lapse[1:])
    kick=4*(np.sin(KY*h/2)/h)**2
    init=np.exp(-.5*((x-21)/4.2)**2+1j*KX*x)
    k=2*np.pi*np.fft.fftfreq(len(x),h)
    w=np.sqrt(4*(np.sin(k*h/2)/h)**2+kick)
    initial=np.column_stack([np.empty(len(x),complex),init])
    for p in np.unique(period):
        prev=np.fft.ifft(np.fft.fft(init)*np.exp(1j*w*(p*quantum)))
        initial[period==p,0]=prev[period==p]
    initial[[0,-1]]=0.
    ref,_,events=simulate(initial,period,lapse,face,h,proper,kick,end,-1)
    raw=ref.tobytes()
    matches=0
    witnesses=[]
    for seed in range(6):
        hist,witness,produced=simulate(initial,period,lapse,face,h,proper,kick,end,350900+seed)
        matches+=hist.tobytes()==raw
        witnesses.append(witness.tolist())
        assert produced==events
    active=float(np.max(abs(ref)))
    def phase_probe(location,t0,t1):
        i=round(location/h)
        t=np.arange(end//period[i]+1)*period[i]*quantum
        y=ref[i,1:1+len(t)]
        mask=(t>=t0)&(t<=t1)
        return phase_slope(t[mask],y[mask])
    def field(t):
        tick=round(t/quantum)
        return ref[np.arange(len(x)),tick//period+1]
    before=phase_probe(24,0,14)
    after=phase_probe(53,21,36)
    proper_ratio=(after/lapse[round(53/h)])/before
    incident,k0=wave_angle(x,field(6),(x>=14)&(x<=35),h)
    transmitted,k1=wave_angle(x,field(36),(x>=46)&(x<=72),h)
    expected_incident=float(np.degrees(np.arctan2(KY,KX)))
    at_right=lapse[round(53/h)]
    expected_out=float(np.degrees(np.arcsin(at_right*KY/np.sqrt(KX*KX+KY*KY))))
    snell_from_measured_incident=float(np.degrees(np.arcsin(at_right*np.sin(np.deg2rad(incident)))))
    return dict(h=h,phase_samples=len(x),proper_step=proper,stepped=step,events=int(events),schedules=6,
                matches=matches,sha256=hashlib.sha256(raw).hexdigest(),witnesses=witnesses,
                slow_versions=int(end//3),fast_versions=int(end//period[round(53/h)]),
                coordinate_frequencies=[before,after],proper_frequencies=[before,after/at_right],
                frequency_ratio=proper_ratio,expected_ratio=1/at_right,
                incident_angle=incident,transmitted_angle=transmitted,expected_incident=expected_incident,
                expected_transmitted=expected_out,snell_from_measured_incident=snell_from_measured_incident,
                longitudinal_wavenumbers=[k0,k1],finite=bool(np.isfinite(ref).all()),peak_field=active)


if __name__=="__main__":
    start=time.process_time()
    rows=[run(.2,.06,True),run(.1,.03,True),run(.1,.03,False)]
    fine=rows[1]
    gates=dict(schedule=all(x["matches"]==x["schedules"] and all(w[0]>=0 for w in x["witnesses"]) for x in rows),
               finite=all(x["finite"] for x in rows),
               frequency=abs(fine["frequency_ratio"]/fine["expected_ratio"]-1)<=.02,
               angle=abs(fine["transmitted_angle"]-fine["expected_transmitted"])<=2.)
    result=dict(rows=rows,gates=gates,cpu_seconds=time.process_time()-start)
    path=Path(__file__).resolve().parents[1]/"results"/"EX-GC35I"
    path.mkdir(parents=True,exist_ok=True)
    (path/"results.json").write_text(json.dumps(result,indent=2,default=lambda x:x.item())+"\n")
    print(result,flush=True)
