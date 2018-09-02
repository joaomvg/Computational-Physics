#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>

#define NMAX 10000

/* � neces�rio linkar este programa com linha.o e gnuplot_i.o */

extern void linhadecampo(double *); 
extern double distancia(double *);
extern double *V(double*);
extern void *RungK4(double *,double);
extern double Kurv(double *);

static double *Tang(double *);

static double *B(double *,double ); 

double *f(double *,double ); 

double curv(double *, double ); 

void RK4(double *,double , double); 

int main(){
  
  /* Nmax: n�mero m�ximo de pontos a calcular */
  
  /* Fun��es:
   * linhadecampo: constroi trajectoria linha de campo
   * distancia: calcula distancia euclidiana entre duas posi��es
   * V:devolve o vector unit�rio tangente � linha de campo
   * RungK4: integra��o Runge-Kutta para a linha de campo
   * Kurv: calcula a curvatura no programa linha de campo
   * Tang: devolve o vector unit�rio tangente � linha de campo neste programa
   * B: fun��o campo magn�tico, devolve vector campo magn�tico
   * f: devolve o vector 6-dimensional velocidade generalizada, ou seja, devolve dp/dt
   * curv: devolve a curvatura para a traject�ria da part�cula
   * RK4: Runge-Kutta de 4� ordem para o problema da integra��o da traject�ria da part�cula
   */
  
  /* Variaveis:
   * t=tempo;
   * T=tempo total requerido
   * dt=passo da integra��o
   * p[6]= vector 6-dimensional posi��o generalizada
   * v[6],r[6],w[6],q[6], vectores auxiliares
   * erro[6]= array que cont�m os erros relativos correspondentes a cada coordenada
   * c= par�metro auxiliar
   * K=curvatura da traject�ria num determinado ponto
   * alpha=constante adimensional da equa��o
   * modV=norma da velocidade para cada ponto
   * modulo= norma da velocidade inicial
   * n=conta o n�mero de pontos calculados
   * u=conta o n�mero de pontos que n�o passaram no teste
   * j=conta o n�mero de subitera��es efectuadas durante o teste
   * s[0],f[0]=strings utilizadas para guardar op��es do utilizador
   */
  
  double t=0.0,T,dt,v[6],p[6],r[6],w[6],q[6],erro[6],c,K,alpha,modV=0,modulo=0;
  int i,j=1,u=0,n=0,e=0;
  char s[0],f[0];/* string utilizada para averiguar se o utilizador continua o prog ou nao*/
  FILE *gnu;

  gnu=popen("gnuplot","w");
  
  printf("alpha=");
  scanf("%lf",&alpha);
  printf("Posi��o inicial (unidades R):\n");
  printf("X(>1)=");
  scanf("%lf",&p[0]);
  printf("Y(>1)=");
  scanf("%lf",&p[1]);
  printf("Z(>1)=");
  scanf("%lf",&p[2]);
  printf("Vector linha de campo nesse ponto:\n");
  for(i=0;i<3;++i){
    printf("vecT%c=%lf\n",'x'+i,Tang(p)[i]);
  }
  printf("Velocidade inicial: \n");
  printf("Vx=");
  scanf("%lf",&p[3]);
  printf("Vy=");
  scanf("%lf",&p[4]);
  printf("Vz=");
  scanf("%lf",&p[5]);
  printf("tempo total=");
  scanf("%lf", &T);
  printf("Pretende visualizar linha de campo?(y or n):");
  scanf("%s",f);
  if(f[0]=='y')
    linhadecampo(p);
  alpha=pow(alpha,3)*2*M_PI;
  
  for(i=0;i<3;++i){
    v[i]=p[i];
    v[i+3]=p[i]+p[i+3];
  }
  
  /* modulo da velocidade inicial*/
  for(i=0;i<3;++i){
    modulo+=pow(p[i+3],2);
  }
  modulo=sqrt(modulo);
  
  FILE *fi;
  fi=fopen("trajelect.dat","w");
  fprintf(fi,"#alpha=%lf\n#X=%lf,Y=%lf,Z=%lf\n#Vx=%lf,Vy=%lf,Vz=%lf\n",alpha,p[0],p[1],p[2],p[3],p[4],p[5]); /* imprime no  ficheiro as condi��es iniciais*/
  
  while(t<T && n<NMAX){
    fprintf(fi,"%lf %lf %lf #modulo(velocidade)=%lf\n",p[0],p[1],p[2],sqrt(modV));
   
    /* primeiro controlo do passo, dt � ajustado de acordo com a curvatura */
    K=curv(p,alpha);
    if (K!=0.0){ 
      dt=0.01/K;
      if(dt<0.0008)
	dt=0.0008;
    }
    else
      dt=0.01;
    
    for(i=0;i<6;++i){
      r[i]=p[i];
    }
    
    RK4(r,alpha,dt);
    
    while(j<=10){ /* segundo controlo do passo, dt � variado de modo a diminuir o erro relativo entre itera��es */
      c=dt/pow(2,j);
      for(i=0;i<6;++i){
	w[i]=p[i];
      }
      
      for(i=0;i<pow(2,j);++i){
	RK4(w,alpha,c);
      }
    
      /* calculo do erro relativo para cada coordenada */
      erro[0]=0;
      for(i=0;i<6;++i){
	if (w[i]!=0.0)
	  erro[i]=fabs((w[i]-r[i])/w[i]);
	else 
	  erro[i]=erro[i-1];
      }
      
      for(i=0;i<6;++i){
	q[i]=r[i];
	r[i]=w[i];
      }
      
      e=0;
      for(i=0;i<6;++i){
	if (erro[i]>0.001 ){ /* se o erro relativo de uma qualquer coordenada for superior a 0.1% o ciclo termina */
	  break;
	}
	else
	  ++e;
      }
      if(e==6)
	break;
      j+=1;
    }
    
    for(i=0;i<6;++i){              /* m�todo utilizado para diminuir o erro associado � posi��o */
      p[i]=(pow(2,5)*w[i]-q[i])/(pow(2,5)-1.0);
    }
    
    modV=0; 
    for(i=0;i<3;++i){
      modV+=pow(p[i+3],2);
    }
    /* ajuste da velocidade */
    for(i=0;i<3;++i){
      p[i+3]*=(modulo/sqrt(modV));
    }
    /* se j>10, ou seja, se n�o foram calculados com sucesso pontos ap�s 10 subinterac��es, u=u+1. u final= n�mero de pontos que foram mal calculados*/
    if (j>10)
      u+=1;
    j=1;
    n+=1;
    t+=dt;
  }
  
  if (u!=0)
    printf("N�o foram calculados com sucesso %d pontos\n",u);
  fclose(fi);
  fflush(fi);
  /* utiliza��o do gnuplot para visualiza��o das traject�rias*/
  
  fprintf(gnu,"\nset title 'Ri=(%.2lf,%.2lf,%.2lf) ; Vi=(%.2lf,%.2lf,%.2lf)'\n",v[0],v[1],v[2],v[3]-v[0],v[4]-v[1],v[5]-v[2]);
  fprintf(gnu,"\n set parametric\n set urange [0:%lf]\n set vrange [0:%lf]",2*M_PI,M_PI);
  fprintf(gnu,"\n set arr 1 from 0,0,-1.3 to 0,0,1.3 lt 6 lw 3\n set style arrow 1 head filled");/* vector momento magnetico*/
  fprintf(gnu,"\n set arr 2 from %lf,%lf,%lf to %lf,%lf,%lf lt 4 lw 1\n",v[0],v[1],v[2],v[3],v[4],v[5]);
  fprintf(gnu,"\n set style arrow 1 head filled\n set style arrow 2 head filled\n set isosamples 10\n");
  if(f[0]=='y')
    fprintf(gnu,"\n spl cos(u)*sin(v),sin(u)*sin(v),cos(v) title 'dipolo magnetico' lw 1 , 'trajelect.dat' title 'traject�ria do electr�o' w l 3, 'linha.dat' title 'linha de campo' w l 6\n");
  else 
    fprintf(gnu,"\n set mouse\n spl cos(u)*sin(v),sin(u)*sin(v),cos(v) title 'dipolo magnetico' lw 1 , 'trajelect.dat' title 'traject�ria do electr�o' w l 3\n pause mouse \n");
  fflush(gnu);
  printf("Deseja continuar? (y or n): ");
  scanf("%s",s);
  if(s[0]!='n')
    main();

  fclose(gnu);
  return 0;
}

static double *Tang(double *p){
  double r2=pow(p[0],2)+pow(p[1],2)+pow(p[2],2);
  static double v[3];
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

static double *B(double *p,double a){
  double R=sqrt(pow(p[0],2)+pow(p[1],2)+pow(p[2],2));
  double m[3]={0.0,0.0,1.0};
  double b=pow(R,5);
  double c=p[0]*m[0]+p[1]*m[1]+p[2]*m[2];
  static double B[3];
  int i;
  for(i=0;i<3;++i){
    B[i]=(3*p[i]*c-pow(R,2)*m[i])*a/b;
    }
  return B;
}

double *f(double *p,double a){  /* f devolve um array cujas componentes sao as derivadas do vector posic�o generalizada p */
  double g[3]={p[4]*B(p,a)[2]-p[5]*B(p,a)[1],
	       p[5]*B(p,a)[0]-p[3]*B(p,a)[2],
	       p[3]*B(p,a)[1]-p[4]*B(p,a)[0]};
  
  double h[6]={p[3],p[4],p[5],g[0],g[1],g[2]};/* p[] � um vector generalizado, p[0..2]=posicao  e p[3..5]=velocidade*/
  return h;
}

double curv(double *p, double a){ /*calculo da curvatura para determina��o do passo dt mais adequado */
  double NORMAacel=sqrt(pow(f(p,a)[3],2)+pow(f(p,a)[4],2)+pow(f(p,a)[5],2)); /* norma da acelera��o*/
  double NORMAvel=sqrt(pow(p[3],2)+pow(p[4],2)+pow(p[5],2)); /* norma da velocidade */
  return NORMAacel/NORMAvel;
}

void RK4(double *p,double a, double dt){
  double k1[6],k2[6],k3[6],k4[6],r[6];
  int i;
  for(i=0;i<6;++i){
    r[i]=p[i];
  }
  for(i=0;i<6;++i){
    k1[i]=f(r,a)[i]*dt;
  }
  for(i=0;i<6;++i){
    r[i]=p[i]+k1[i]/2.0;
  }
  for(i=0;i<6;++i){
    k2[i]=f(r,a)[i]*dt;
  }
  for(i=0;i<6;++i){
    r[i]=p[i]+k2[i]/2.0;
  }
  for(i=0;i<6;++i){
    k3[i]=f(r,a)[i]*dt;
  }
  for(i=0;i<6;++i){
    r[i]=p[i]+k3[i];
  }
  for(i=0;i<6;++i){
    k4[i]=f(r,a)[i]*dt;
  }
  for(i=0;i<6;++i){
    p[i]+=k1[i]/6.0+k2[i]/3.0+k3[i]/3.0+k4[i]/6.0;
  }
}






