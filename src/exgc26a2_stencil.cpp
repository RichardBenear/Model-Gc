// Independently implemented Skyrme static functional and exact gradient.
#include <cmath>
#include <algorithm>
#include <vector>
extern "C" void evaluate(int n, double h, int mode, const double* p, double* g, double* parts) {
  const int N=n*n*n, pad=mode?3:0; const double norm=h*h*h/(12*std::acos(-1.)*std::acos(-1.));
  std::fill(g,g+4*N,0.); parts[0]=parts[1]=0.;
  const double coeff[4]={-11./6,3.,-1.5,1./3};
  auto idx=[&](int x,int y,int z){return (x<0||y<0||z<0||x>=n||y>=n||z>=n)?-1:((x*n+y)*n+z);};
  const int count=mode?8:2, terms=mode?4:2;
  for(int x=-pad;x<n+pad;++x)for(int y=-pad;y<n+pad;++y)for(int z=-pad;z<n+pad;++z){
    int base=idx(x,y,z);
    for(int s=0;s<count;++s){
      int dir[3]; for(int a=0;a<3;++a)dir[a]=mode?((s>>a)&1?1:-1):(s?1:-1);
      if(!mode && ((s && (x==n-1||y==n-1||z==n-1))||(!s && (x==0||y==0||z==0))))continue;
      double d[3][4]={{0}}, D[3][3]={{0}}, w[3][4]; int ids[3][4];
      for(int a=0;a<3;++a)for(int t=0;t<terms;++t){
        int q[3]={x,y,z};q[a]+=t*dir[a];int j=idx(q[0],q[1],q[2]);ids[a][t]=j;
        w[a][t]=dir[a]*(mode?coeff[t]:(t?1.:-1.))/h;
        for(int b=0;b<4;++b)d[a][b]+=w[a][t]*(j<0?(b==0?1.:0.):p[4*j+b]);
      }
      for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int k=0;k<4;++k)D[a][b]+=d[a][k]*d[b][k];
      double tr=D[0][0]+D[1][1]+D[2][2], e4=0.;
      for(int a=0;a<3;++a)for(int b=a+1;b<3;++b)e4+=D[a][a]*D[b][b]-D[a][b]*D[a][b];
      const double f=norm/count;parts[0]+=f*tr;parts[1]+=f*e4;
      for(int a=0;a<3;++a)for(int k=0;k<4;++k){double v=2*d[a][k]*(1+tr);for(int b=0;b<3;++b)v-=2*D[a][b]*d[b][k];
        for(int t=0;t<terms;++t)if(ids[a][t]>=0)g[4*ids[a][t]+k]+=f*v*w[a][t];}
    }
  }
}

// Exact same Euler-Lagrange dynamics as inherited two-orientation stencil,
// independently accumulated into per-node mass matrices and RHS.
extern "C" void acceleration(int n,double h,const double* p,const double* v,double* out){
  const int N=n*n*n, stride[3]={n*n,n,1};
  std::vector<double> M(16*N,0.),R(4*N,0.);
  for(int i=0;i<N;++i)for(int k=0;k<4;++k)M[16*i+5*k]=1.;
  for(int x=0;x<n;++x)for(int y=0;y<n;++y)for(int z=0;z<n;++z){
    int i=(x*n+y)*n+z;
    double vv=0.;for(int k=0;k<4;++k)vv+=v[4*i+k]*v[4*i+k];
    for(int s=-1;s<=1;s+=2){
      if((s==1&&(x==n-1||y==n-1||z==n-1))||(s==-1&&(x==0||y==0||z==0)))continue;
      double d[3][4],dv[3][4],D[3][3]={{0}},vd[3]={0},dvv[3]={0},sdd=0.;
      for(int a=0;a<3;++a)for(int k=0;k<4;++k){int j=i+s*stride[a];d[a][k]=s*(p[4*j+k]-p[4*i+k])/h;dv[a][k]=s*(v[4*j+k]-v[4*i+k])/h;vd[a]+=v[4*i+k]*d[a][k];dvv[a]+=v[4*i+k]*dv[a][k];sdd+=d[a][k]*dv[a][k];}
      for(int a=0;a<3;++a)for(int b=0;b<3;++b)for(int k=0;k<4;++k)D[a][b]+=d[a][k]*d[b][k];
      double tr=D[0][0]+D[1][1]+D[2][2];
      for(int k=0;k<4;++k){
        for(int l=0;l<4;++l){double m=(k==l?tr:0.);for(int a=0;a<3;++a)m-=d[a][k]*d[a][l];M[16*i+4*k+l]+=.5*m;}
        double mdv=2*sdd*v[4*i+k];for(int a=0;a<3;++a)mdv-=dv[a][k]*vd[a]+d[a][k]*dvv[a];R[4*i+k]-=.5*mdv;
        for(int a=0;a<3;++a){double flux=(vv-1-tr)*d[a][k]-vd[a]*v[4*i+k];for(int b=0;b<3;++b)flux+=D[a][b]*d[b][k];flux*=.5*s/h;R[4*i+k]-=flux;R[4*(i+s*stride[a])+k]+=flux;}
      }
    }
  }
  std::fill(out,out+4*N,0.);
  for(int x=1;x<n-1;++x)for(int y=1;y<n-1;++y)for(int z=1;z<n-1;++z){
    int i=(x*n+y)*n+z;double L[4][4]={{0}},a[4],b[4],vv=0.;
    for(int k=0;k<4;++k){for(int l=0;l<=k;++l){double t=M[16*i+4*k+l];for(int j=0;j<l;++j)t-=L[k][j]*L[l][j];L[k][l]=k==l?std::sqrt(t):t/L[l][l];}vv+=v[4*i+k]*v[4*i+k];}
    for(int k=0;k<4;++k){double x1=R[4*i+k],x2=p[4*i+k];for(int j=0;j<k;++j){x1-=L[k][j]*a[j];x2-=L[k][j]*b[j];}a[k]=x1/L[k][k];b[k]=x2/L[k][k];}
    for(int k=3;k>=0;--k){for(int j=k+1;j<4;++j){a[k]-=L[j][k]*a[j];b[k]-=L[j][k]*b[j];}a[k]/=L[k][k];b[k]/=L[k][k];}
    double pa=0.,pb=0.;for(int k=0;k<4;++k){pa+=p[4*i+k]*a[k];pb+=p[4*i+k]*b[k];}
    double mu=(-vv-pa)/pb;for(int k=0;k<4;++k)out[4*i+k]=a[k]+mu*b[k];
  }
}
