// EX-GC23: the inherited EX-GC22 chemistry and a fused explicit wave step.
#include "exgc22_kernel.cpp"

extern "C" void wave_step(int n, double h, double dt, double K, double tauJ,
                          double E, const double* u, double* jx, double* jy,
                          double* jz, double* gn) {
    const int N=n*n*n;
    const double* G=u+3*N; const double* W=u+2*N;
    for(int x=0;x<n;++x) for(int y=0;y<n;++y) for(int z=0;z<n;++z){
        int p=(x*n+y)*n+z, py=(x*(n+1)+y)*n+z, pz=(x*n+y)*(n+1)+z;
        gn[p]=G[p]+dt*(-(jx[p+n*n]-jx[p]+jy[py+n]-jy[py]+jz[pz+1]-jz[pz])/h-K*G[p]*W[p]);
    }
    // Face velocities use the old G as in EX-GC22. Ghost boundary is 1+Ez.
    for(int x=0;x<=n;++x) for(int y=0;y<n;++y) for(int z=0;z<n;++z){
        int p=(x*n+y)*n+z; double bc=1+E*(z-(n-1)*.5)*h;
        double lo=x?G[p-n*n]:bc, hi=x<n?G[p]:bc;
        jx[p]+=dt*(-4*(hi-lo)/h-jx[p]/tauJ);
    }
    for(int x=0;x<n;++x) for(int y=0;y<=n;++y) for(int z=0;z<n;++z){
        int p=(x*(n+1)+y)*n+z, q=(x*n+y)*n+z;
        double bc=1+E*(z-(n-1)*.5)*h;
        double lo=y?G[q-n]:bc, hi=y<n?G[q]:bc;
        jy[p]+=dt*(-4*(hi-lo)/h-jy[p]/tauJ);
    }
    for(int x=0;x<n;++x) for(int y=0;y<n;++y) for(int z=0;z<=n;++z){
        int p=(x*n+y)*(n+1)+z,q=(x*n+y)*n+z;
        double lo=z?G[q-1]:1+E*(-(n+1)*.5)*h;
        double hi=z<n?G[q]:1+E*((n+1)*.5)*h;
        jz[p]+=dt*(-4*(hi-lo)/h-jz[p]/tauJ);
    }
}
