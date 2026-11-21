// Eight nearest-neighbor orientations, same discrete static and kinetic action.
#include <cmath>
#include <algorithm>
#include <vector>
static constexpr double w=0.125;
template<bool dyn> void core(int n,double h,const double* p,const double* v,double f,int power,double* out,double* parts){
  const int N=n*n*n, stride[3]={n*n,n,1};const double norm=h*h*h/(12*std::acos(-1.)*std::acos(-1.));
  std::vector<double> M(dyn?16*N:0,0.),R(dyn?4*N:0,0.);std::fill(out,out+4*N,0.);if(parts)parts[0]=parts[1]=0.;
  for(int x=0;x<n;++x)for(int y=0;y<n;++y)for(int z=0;z<n;++z){
    int i=(x*n+y)*n+z, coord[3]={x,y,z};double vv=0.;if(dyn){for(int k=0;k<4;++k){vv+=v[4*i+k]*v[4*i+k];M[16*i+5*k]=1.;}R[4*i]-=.5*f*(x-(n-1)*.5)*h*(power==2?2*(1-p[4*i]):1.);}
    // Prepare the six derivatives once per site.
    double ds[3][2][4]={},dvs[3][2][4]={};bool valid[3][2];
    for(int a=0;a<3;++a)for(int b=0;b<2;++b){int s=2*b-1;valid[a][b]=(coord[a]+s>=0&&coord[a]+s<n);if(!valid[a][b])continue;int j=i+s*stride[a];for(int k=0;k<4;++k){ds[a][b][k]=s*(p[4*j+k]-p[4*i+k])/h;if(dyn)dvs[a][b][k]=s*(v[4*j+k]-v[4*i+k])/h;}}
    for(int o=0;o<8;++o){int b[3]={o&1,(o>>1)&1,(o>>2)&1};if(!valid[0][b[0]]||!valid[1][b[1]]||!valid[2][b[2]])continue;
      const double* d[3]={ds[0][b[0]],ds[1][b[1]],ds[2][b[2]]};const double* dv[3]={dvs[0][b[0]],dvs[1][b[1]],dvs[2][b[2]]};double D[3][3]={},vd[3]={},dvv[3]={},sdd=0.;
      for(int a=0;a<3;++a){for(int c=a;c<3;++c){double t=0.;for(int k=0;k<4;++k)t+=d[a][k]*d[c][k];D[a][c]=D[c][a]=t;}if(dyn)for(int k=0;k<4;++k){vd[a]+=v[4*i+k]*d[a][k];dvv[a]+=v[4*i+k]*dv[a][k];sdd+=d[a][k]*dv[a][k];}}
      double tr=D[0][0]+D[1][1]+D[2][2];
      if(!dyn){double e4=0.;for(int a=0;a<3;++a)for(int c=a+1;c<3;++c)e4+=D[a][a]*D[c][c]-D[a][c]*D[a][c];parts[0]+=norm*w*tr;parts[1]+=norm*w*e4;}
      for(int k=0;k<4;++k){
        if(dyn){for(int l=0;l<4;++l){double m=(k==l?tr:0.);for(int a=0;a<3;++a)m-=d[a][k]*d[a][l];M[16*i+4*k+l]+=w*m;}double mdv=2*sdd*v[4*i+k];for(int a=0;a<3;++a)mdv-=dv[a][k]*vd[a]+d[a][k]*dvv[a];R[4*i+k]-=w*mdv;}
        for(int a=0;a<3;++a){double q=(1+tr)*d[a][k];for(int c=0;c<3;++c)q-=D[a][c]*d[c][k];int s=2*b[a]-1,j=i+s*stride[a];if(dyn){double flux=w*s/h*(-q+vv*d[a][k]-vd[a]*v[4*i+k]);R[4*i+k]-=flux;R[4*j+k]+=flux;}else{double flux=2*norm*w*s/h*q;out[4*i+k]-=flux;out[4*j+k]+=flux;}}
      }
    }
  }
  if(!dyn)return;
  for(int x=1;x<n-1;++x)for(int y=1;y<n-1;++y)for(int z=1;z<n-1;++z){int i=(x*n+y)*n+z;double L[4][4]={},a[4],b[4],vv=0.;
    for(int k=0;k<4;++k){for(int l=0;l<=k;++l){double t=M[16*i+4*k+l];for(int j=0;j<l;++j)t-=L[k][j]*L[l][j];L[k][l]=k==l?std::sqrt(t):t/L[l][l];}vv+=v[4*i+k]*v[4*i+k];}
    for(int k=0;k<4;++k){double t=R[4*i+k],u=p[4*i+k];for(int j=0;j<k;++j){t-=L[k][j]*a[j];u-=L[k][j]*b[j];}a[k]=t/L[k][k];b[k]=u/L[k][k];}
    for(int k=3;k>=0;--k){for(int j=k+1;j<4;++j){a[k]-=L[j][k]*a[j];b[k]-=L[j][k]*b[j];}a[k]/=L[k][k];b[k]/=L[k][k];}double pa=0.,pb=0.;for(int k=0;k<4;++k){pa+=p[4*i+k]*a[k];pb+=p[4*i+k]*b[k];}double mu=(-vv-pa)/pb;for(int k=0;k<4;++k)out[4*i+k]=a[k]+mu*b[k];
  }
}
extern "C" void evaluate8(int n,double h,const double* p,double* g,double* parts){core<false>(n,h,p,nullptr,0.,1,g,parts);}
extern "C" void acceleration8(int n,double h,const double* p,const double* v,double f,double* out){core<true>(n,h,p,v,f,1,out,nullptr);}
extern "C" void acceleration8_coreforce(int n,double h,const double* p,const double* v,double f,double* out){core<true>(n,h,p,v,f,2,out,nullptr);}
