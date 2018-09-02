#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

#define N 32 //grelha NxN
#define Flmax 10000

// Hamiltoniano: H=-JS(i)S(i+1)

static gsl_rng *rnduniform(int j){ //j:seed
  int seed;
  gsl_rng_type *T;
  static gsl_rng *a;
  T=gsl_rng_default;
  a=gsl_rng_alloc(T);
  srand(time(NULL));
  seed=rand()%(100*j);
  gsl_rng_set(a,seed);
  return a;
}

int viz(int n,int j){ //coordenada do vizinho; j=+-1
  if(n+j<0 || n+j>N-1)
    return N-1-n;
  else
    return n+j;
}

int main(){
  int i,j,flips=0; //flips: conta o numero de flips 
  double **spin,sum,T,M=0,E=0; // M: magnetização, Mmed: mag media; E=energia, Emed: energia media 
  gsl_rng *x,*y;
  FILE *gnu,*EnMag;
  
  gnu=popen("gnuplot","w");
  EnMag=fopen("isingEMtempo.dat","w");
  
  printf("T=");
  scanf("%lf",&T);
  
  spin=calloc(N,sizeof(double));
  for(i=0;i<N;++i){
    spin[i]=calloc(N,sizeof(double));
  }
  
  //algoritmo de Metropolis

  //inicialização
  
   for(i=0;i<N;++i){
      for(j=0;j<N;++j){
        spin[i][j]=1;
      }
    }
   
   x=rnduniform(10);
   y=rnduniform(100);
       
   
   while(flips<Flmax){
     
     for(i=0;i<N;++i){
       for(j=0;j<N;++j){
         E+=-spin[i][j]/(2.0*N*N)*(spin[viz(i,1)][j]+spin[viz(i,-1)][j]+spin[i][viz(j,1)]+spin[i][viz(j,-1)]);// energia por spin
         M+=spin[i][j]/(N*N);// magnetização por spin         
       }
     }
       
     //escolher aleatoriamente um spin
     i=gsl_rng_uniform(x)*N;
     j=gsl_rng_uniform(x)*N;

     //somar sobre os vizinhos
     sum=spin[viz(i,1)][j]+spin[viz(i,-1)][j]+spin[i][viz(j,1)]+spin[i][viz(j,-1)];
     
     if(sum*spin[i][j]<0)
       spin[i][j]*=-1;
     else{
       if(exp(-2*sum*spin[i][j]/T)>gsl_rng_uniform(y))
         spin[i][j]*=-1;
     }
     
     fprintf(EnMag,"%lf %lf %lf\n",(double)flips,E,M);
     flips+=1;
     E=0;
     M=0;
   }

   free(spin);
   fclose(gnu);
   fclose(EnMag);
   
   return 0;
}


