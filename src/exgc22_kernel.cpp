// Fused EX-GC22 chemistry stencil. Matches numpy.gradient edge_order=1 and
// the reflecting seven-point Laplacian of exgc09_lattice.py.
#include <cstddef>
extern "C" void chem(int n, double h, double K, double A, double B, double DY,
                    double DW, double C, const double* u, const double* V, double* out) {
    const int N=n*n*n; const int stride[3]={n*n,n,1}; const double h2=h*h;
    const double base[3]={A,B/A,0.};
    for(int x=0;x<n;++x) for(int y=0;y<n;++y) for(int z=0;z<n;++z){
        int p=(x*n+y)*n+z, pos[3]={x,y,z}; double lap[3]={},adv[3]={};
        for(int a=0;a<3;++a){
            int lo=pos[a]>0?p-stride[a]:p, hi=pos[a]<n-1?p+stride[a]:p;
            double width=(lo==p||hi==p)?h:2*h;
            for(int j=0;j<3;++j){
                const double* f=u+j*N; const double* v=V+a*N;
                lap[j]+=(f[lo]+f[hi]-2*f[p])/h2;
                adv[j]-=((f[hi]-base[j])*v[hi]-(f[lo]-base[j])*v[lo])/width;
            }
        }
        double X=u[p],Y=u[N+p],W=u[2*N+p],G=u[3*N+p],s=K*G*W,autoXY=X*X*Y;
        out[p]=lap[0]+A-(1+B)*X+autoXY+s+adv[0];
        out[N+p]=DY*lap[1]+B*X-autoXY+adv[1];
        out[2*N+p]=DW*lap[2]-W+C*Y*W*W-W*W*W+adv[2];
    }
}
