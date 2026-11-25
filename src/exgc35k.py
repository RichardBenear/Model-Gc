"""GC35K: static Y1 mass solve with a proper-time event packet."""
import argparse
import json
from pathlib import Path
import time
import exgc34p as P
import numpy as np
import scipy.fft as fft
from scipy.ndimage import map_coordinates
from numba import njit

GM=.095
IMPACT=7.5
CARRIER=4.
RADIUS=1.4
WAIST=1.65
STOP=50.
CLOCK_TICKS=2048


def metric():
    n=160
    h=1.
    x=(np.arange(n)-n//2)*h
    r2=x[:,None,None]**2+x[None,:,None]**2+x[None,None,:]**2
    rho=np.exp(-r2/(2*RADIUS**2))
    rho*=8*np.pi*GM/(P.KAP*rho.sum()*h**3)
    f=fft.rfftn(rho)
    k=2*np.sin(np.pi*np.fft.fftfreq(n))/h
    kr=k[:n//2+1].copy()
    kr[-1]=abs(kr[-1])
    K2=k[:,None,None]**2+k[None,:,None]**2+kr[None,None,:]**2
    K2[0,0,0]=1
    lapse=-P.KAP*f/(2*K2)
    lapse[0,0,0]=0
    modes=np.array([[1,0,0],[1,2,3],[7,-4,5],[20,20,17],[4,5,0]])
    i,j,l=modes.T
    Kv=np.column_stack([k[i],k[j],kr[l]])
    probe=lapse[i,j,l]
    state=np.zeros((len(i),10),complex)
    state[:,0]=probe
    state[:,4:7]=-2*probe[:,None]
    source=np.zeros_like(state)
    source[:,0]=-2*P.KAP*f[i,j,l]
    audit=float(np.linalg.norm(P.mv(P.kernel_batch(Kv,0.),state)-source)/np.linalg.norm(source))
    assert audit<1e-12, audit
    pot=fft.irfftn(lapse,s=rho.shape,axes=(0,1,2))[:,:,n//2].copy()
    return pot,h,audit


@njit
def read(history,period,done,j,tick,deadline,order):
    v=tick//period[j]
    if v>done[j] or done[j]-v>=3:
        raise ValueError("required earlier versions unavailable")
    u=history[j,v%5]
    x=(tick%period[j])/period[j]
    if x==0:
        return u
    before=history[j,(v-1)%5]
    if order==1:
        return u+x*(u-before)
    if (v+1)*period[j]<deadline:
        if done[j]<v+1:
            raise ValueError("missing bracketing version")
        upper=history[j,(v+1)%5]
        return u+.5*x*(upper-before)+.5*x*x*(upper-2*u+before)
    if v==0:
        return u+x*(u-before)
    past=history[j,(v-2)%5]
    return u+x*(u-before)+.5*x*(x+1)*(u-2*before+past)


@njit
def events(u,old,period,lapse,A,nx,ny,step,h,ticks,order=2,reverse=False):
    n=nx*ny
    history=np.zeros((n,5),np.complex128)
    history[:,0]=u
    history[:,4]=old
    size=np.max(period)+1
    head=np.full(size,-1,np.int64)
    links=np.empty(n,np.int64)
    done=np.zeros(n,np.int64)
    for k in range(n):
        i=n-1-k if reverse else k
        at=period[i]%size
        links[i]=head[at]
        head[at]=i
    count=0
    errors=0.
    max_gap=0
    for when in range(1,ticks+1):
        at=when%size
        i=head[at]
        head[at]=-1
        while i>=0:
            next_i=links[i]
            pi=period[i]
            earlier=when-pi
            v=done[i]
            if when!=(v+1)*pi:
                raise ValueError("event queued at the wrong time")
            a,b=i//ny,i%ny
            z=history[i,v%5]
            new=0j
            if 0<a<nx-1 and 0<b<ny-1:
                flux=0j
                for j in (i-ny,i+ny,i-1,i+1):
                    neighbor=read(history,period,done,j,earlier,when,order)
                    face=2*lapse[i]*lapse[j]/(lapse[i]+lapse[j])
                    flux+=face*(neighbor-z)
                    gap=earlier-(earlier//period[j])*period[j]
                    max_gap=max(max_gap,gap)
                new=2*z-history[i,(v-1)%5]+step[i]*step[i]*flux/(lapse[i]*A[i]*A[i]*h*h)
            history[i,(v+1)%5]=new
            done[i]+=1
            count+=1
            if when+pi<=ticks:
                nxt=(when+pi)%size
                links[i]=head[nxt]
                head[nxt]=i
            i=next_i
    output=np.empty(n,np.complex128)
    for i in range(n):
        output[i]=read(history,period,done,i,ticks,ticks+1,order)
    return output,count,int(np.min(done)),int(np.max(done)),max_gap


def phase_data(h):
    nx,ny=round(72/h)+1,round(46/h)+1
    x=(np.arange(nx)-nx//2)*h
    y=(np.arange(ny)-ny//2)*h
    grid=np.meshgrid(x,y,indexing="ij")
    kx,ky=np.meshgrid(2*np.pi*fft.fftfreq(nx,h),2*np.pi*fft.fftfreq(ny,h),indexing="ij")
    waist=np.exp(-.5*((grid[0]/2.55)**2+((grid[1]-IMPACT)/WAIST)**2)+1j*CARRIER*grid[0])
    spatial=2/h*np.sqrt(np.sin(h*kx/2)**2+np.sin(h*ky/2)**2)
    focused=fft.fft2(waist)*np.exp(1j*spatial*25.)
    init=fft.ifft2(focused)
    vel=fft.ifft2(-1j*spatial*focused)
    acc=fft.ifft2(-spatial**2*focused)
    for arr in (init,vel,acc):
        arr[[0,-1],:]=0
        arr[:,[0,-1]]=0
    return grid,kx,ky,init,vel,acc


def center(grid,u):
    w=abs(u)**2
    return [float(np.sum(z*w)/np.sum(w)) for z in grid]


def directions(h,kx,ky,u):
    p=abs(fft.fft2(u))**2
    forward=kx>CARRIER/2
    p*=forward
    p/=p.sum()
    angle=lambda x,y:float(np.degrees(np.arctan2(np.sum(p*y),np.sum(p*x))))
    return angle(kx,ky),angle(np.sin(kx*h),np.sin(ky*h))


def ray_prediction(grid,n,A,edge):
    x=grid[0][:,0]
    y=grid[1][0,:]
    signed=np.gradient(np.log(A/n),y,axis=1)
    lower,upper=sorted(edge)
    mask=(x>=lower)&(x<=upper)
    profile=np.trapezoid(signed[mask],x[mask],axis=0)
    weight=np.exp(-((y-IMPACT)/WAIST)**2)
    return float(np.degrees(np.sum(weight*profile)/np.sum(weight)))


def run(pot,hg,h,step,mode,duplicate=False,order=2,uniform=False):
    grid,kx,ky,u,v,a=phase_data(h)
    nx,ny=u.shape
    units=step/CLOCK_TICKS
    ticks=round(STOP/units)
    p=map_coordinates(pot,[z/hg+pot.shape[0]//2 for z in grid],order=3,mode="wrap")
    target=1+p if mode!="flat" else np.ones_like(p)
    period=np.rint(CLOCK_TICKS/target).astype(np.int64)
    if uniform:
        period=np.full(period.shape,CLOCK_TICKS,dtype=np.int64)
        lapse=target
    else:
        lapse=CLOCK_TICKS/period
    A=np.sqrt(1-2*p) if mode=="full" else np.ones_like(p)
    dt=units*period
    velocity=v*lapse/A
    accel=a*(lapse/A)**2
    old=u-dt*velocity+.5*dt**2*accel
    args=(u.ravel(),old.ravel(),period.ravel(),lapse.ravel(),A.ravel(),nx,ny,(dt*lapse).ravel(),h,ticks,order)
    out,count,low,high,gap=events(*args)
    same=None
    if duplicate:
        again=events(*args,True)[0]
        same=again.tobytes()==out.tobytes()
    out=out.reshape(u.shape)
    angles=directions(h,kx,ky,out)
    start,end=center(grid,u),center(grid,out)
    prediction=ray_prediction(grid,target,A,[start[0],end[0]])
    return dict(mode=mode,h=h,predictor=order,uniform_clock=uniform,proper_step=step,nx=nx,ny=ny,events=count,site_event_range=[low,high],
                maximum_predictor_offset=float(gap*units),lapse_rounding=float(np.max(abs(lapse-target))),
                phase_angle=angles[0],group_angle=angles[1],initial_centroid=start,final_centroid=end,
                eikonal_angle=prediction,finite=bool(np.isfinite(out).all()),peak=float(np.max(abs(out))),
                norm_ratio=float(np.linalg.norm(out)/np.linalg.norm(u)),tie_order_match=same)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("tag")
    parser.add_argument("--h",type=float,default=.16)
    parser.add_argument("--step",type=float,default=.04)
    parser.add_argument("--mode",choices=("flat","lapse","full"),default="full")
    parser.add_argument("--order",type=int,choices=(1,2),default=2)
    parser.add_argument("--uniform",action="store_true")
    parser.add_argument("--duplicate",action="store_true")
    args=parser.parse_args()
    start=time.process_time()
    pot,hg,audit=metric()
    row=run(pot,hg,args.h,args.step,args.mode,args.duplicate,args.order,args.uniform)
    row.update(static_kernel_relative=audit,GM=GM,impact=IMPACT,
                ideal_full_deg=float(np.degrees(4*GM/IMPACT)),cpu_seconds=time.process_time()-start)
    out=Path(__file__).resolve().parents[1]/"results"/"EX-GC35K"
    out.mkdir(parents=True,exist_ok=True)
    (out/f"{args.tag}.json").write_text(json.dumps(row,indent=2)+"\n")
    print(row,flush=True)


if __name__=="__main__":
    main()
