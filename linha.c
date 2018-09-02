#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#define NMAX2 5000

double distancia(double *p){
  int i;
  double b=0;
  for(i=0;i<3;++i){
    b+=pow(p[i],2);
  }
  return sqrt(b);
}

double *V(double *p){
  double r2=pow(p[0],2)+pow(p[1],2)+pow(p[2],2);
  double v[3];
  double d=sqrt(3*pow(p[2],2)*r2+pow(r2,2));
  if (r2<0.0001){
    v[0]=0;
    v[1]=0;
    v[2]=1;
  }
  else{
    v[0]=3*p[2]*p[0]/d;
    v[1]=3*p[2]*p[1]/d;
    v[2]=(3*pow(p[2],2)-r2)/d;
  }  
  return v;
}

double Kurv(double *p){
  double r2=pow(p[0],2)+pow(p[1],2)+pow(p[2],2);
  double r4=pow(r2,2);
  double r6=pow(r2,3);
  double z2=pow(p[2],2);
  double z4=pow(z2,2);
  double d=3*z2+r2;
  return 3*sqrt(8.5*z2*r4-9.5*z4*r2+r6)/pow(d,2);
}

void RungK4(double *p, double dt){
  double k1[3],k2[3],k3[3],k4[3],r[3];
  int i;
  for(i=0;i<3;++i){
    r[i]=p[i];
  }
  for(i=0;i<3;++i){
    k1[i]=V(r)[i]*dt;
  }
  for(i=0;i<3;++i){
    r[i]=p[i]+k1[i]/2.0;
  }
  for(i=0;i<3;++i){
    k2[i]=V(r)[i]*dt;
  }
  for(i=0;i<3;++i){
    r[i]=p[i]+k2[i]/2.0;
  }
  for(i=0;i<3;++i){
    k3[i]=V(r)[i]*dt;
  }
  for(i=0;i<3;++i){
    r[i]=p[i]+k3[i];
  }
  for(i=0;i<3;++i){
    k4[i]=V(r)[i]*dt;
  }
  for(i=0;i<3;++i){
    p[i]+=k1[i]/6.0+k2[i]/3.0+k3[i]/3.0+k4[i]/6.0;
  }
}

void linhadecampo(double *p){
  FILE *linha;
  double v[3],c,ds;
  int N=0,i,j;
  linha=fopen("linha.dat","w");
  for(j=0;j<2;++j){
    for(i=0;i<3;++i){
      v[i]=p[i];
    }
    N=0;
    while(N<NMAX2){
      if (distancia(v)<1.0)
	break;
      fprintf(linha,"%lf %lf %lf\n",v[0],v[1],v[2]);
      c=Kurv(v);
      if(c!=0)
	ds=0.05/c;
      else
	ds=0.01;
      if(j==1)
	ds=-ds;
      RungK4(v,ds);
      N+=1;
    }
    fprintf(linha,"\n\n");
  }
  fclose(linha);
}






