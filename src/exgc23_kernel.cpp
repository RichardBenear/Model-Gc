// EX-GC23: the inherited EX-GC22 chemistry and a fused explicit wave step.
#include "exgc22_kernel.cpp"

// Same centered chemistry, with the z-interior loop exposed for SIMD.
extern "C" void chem_fast(int n, double h, double K, double A, double B, double DY,
                         double DW, double C, const double* __restrict u,
                         const double* __restrict v, double* __restrict out) {
    const int N=n*n*n; const double h2=h*h; const double y0=B/A;
    const double* __restrict X=u; const double* __restrict Y=u+N; const double* __restrict W=u+2*N;
    const double* __restrict G=u+3*N; const double* __restrict vx=v; const double* __restrict vy=v+N; const double* __restrict vz=v+2*N;
    double* __restrict dx=out;double* __restrict dy=out+N;double* __restrict dw=out+2*N;
    for(int x=0;x<n;++x) for(int y=0;y<n;++y) {
        const int q=(x*n+y)*n, xl=x? -n*n:0, xh=x<n-1?n*n:0, yl=y?-n:0, yh=y<n-1?n:0;
        const double wx=(xl&&xh)?2*h:h, wy=(yl&&yh)?2*h:h;
        // Boundary points use the same one-sided derivative convention.
        auto point=[&](int z,int zl,int zh,double wz){
            int p=q+z; double s=K*G[p]*W[p],autoXY=X[p]*X[p]*Y[p];
            auto rhs=[&](const double* f,double f0,double D){
                double lx=(f[p+xl]+f[p+xh]-2*f[p])/h2,ly=(f[p+yl]+f[p+yh]-2*f[p])/h2,lz=(f[p+zl]+f[p+zh]-2*f[p])/h2;
                double ax=-((f[p+xh]-f0)*vx[p+xh]-(f[p+xl]-f0)*vx[p+xl])/wx;
                double ay=-((f[p+yh]-f0)*vy[p+yh]-(f[p+yl]-f0)*vy[p+yl])/wy;
                double az=-((f[p+zh]-f0)*vz[p+zh]-(f[p+zl]-f0)*vz[p+zl])/wz;
                return D*(lx+ly+lz)+(ax+ay+az);
            };
            dx[p]=rhs(X,A,1)+A-(1+B)*X[p]+autoXY+s;
            dy[p]=rhs(Y,y0,DY)+B*X[p]-autoXY;
            dw[p]=rhs(W,0,DW)-W[p]+C*Y[p]*W[p]*W[p]-W[p]*W[p]*W[p];
        };
        point(0,0,1,h);
        #pragma GCC ivdep
        for(int z=1;z<n-1;++z) point(z,-1,1,2*h);
        point(n-1,-1,0,h);
    }
}

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
