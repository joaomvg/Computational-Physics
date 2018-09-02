#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

/*
 * É necessário compilar com a gsl: gcc Ising2DMp2.c -o "..." -lgsl -lgslcblas -lm
 */

#define N 6 //grelha NxN
#define Flmax 10000 // tempo máximo
#define dT 0.01 //intervalo entre temperaturas

// Hamiltoniano: H=-JS(i)S(i+1)

static gsl_rng *rnduniform(int j){ //função geradora de numeros aleatórios; j:seed
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

int viz(int n,int j){ //esta função devolve a coordenada do vizinho n á esquerda\baixo (j=-1) ou á direita\cima (j=1)
  if(n+j<0 || n+j>N-1)
    return N-1-n;
  else
    return n+j;
}

int main(){
  /*
   * //flips: conta o numero de flips 
   */
  int i,j,n=1,flips=0;
  /*
   *  M: magnetização
   *  Mmed: mag media
   * E=energia
   *  Emed: energia media
   * pE: desvio para E
   * pM: desvio para M
   * T: temperatura (começamos com T=1 uma vez que para temperaturas inferiores quer a magnetização quer a energia são aproximadamente constantes)
   * Tmax: temperatura máxima
   */
    
  double **spin,sum=0,T=1,Tmax,M[4],E[4],Mmed,Emed,pE=1,pM=1; 
  gsl_rng *x,*y;
  FILE *gnu,*EnMag;
  
  gnu=popen("gnuplot","w");
  EnMag=fopen("isingEnMag.dat","w");// ficheiro onde é guardada a informação sobre 

  printf("Tmax=");
  scanf("%lf",&Tmax);
  
  spin=calloc(N,sizeof(double));
  for(i=0;i<N;++i){
    spin[i]=calloc(N,sizeof(double));
  }

  //configuração completamente ordenada
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      spin[i][j]=1;
    }
  }
  
  //algoritmo de Metropolis

  //inicialização
  
   for(i=0;i<N;++i){
      for(j=0;j<N;++j){
        spin[i][j]=1;
      }
    }
   
  while(T<Tmax){
    x=rnduniform(10*n);
    y=rnduniform(100*n);
       
    
    while(flips<Flmax){
              
        for(i=0;i<N;++i){
          for(j=0;j<N;++j){
            E[2]+=-spin[i][j]/(2.0*N*N)*(spin[viz(i,1)][j]+spin[viz(i,-1)][j]+spin[i][viz(j,1)]+spin[i][viz(j,-1)]);// energia por spin
            M[2]+=spin[i][j]/(N*N);// magnetização por spin         
          }
        }

        E[0]+=0.5*E[1]+0.5*E[2];
        M[0]+=0.5*M[1]+0.5*M[2];
        E[3]+=0.5*pow(E[1],2)+0.5*pow(E[2],2); //calculo da media de E^2
        M[3]+=0.5*pow(M[1],2)+0.5*pow(M[2],2); //calculo da media de M^2
        E[1]=E[2];
        M[1]=M[2];
        
            
      //escolher aleatoriamente um spin
      i=gsl_rng_uniform(x)*N;
      j=gsl_rng_uniform(x)*N;

      //somar sobre os vizinhos
      sum+=spin[viz(i,1)][j]+spin[viz(i,-1)][j]+spin[i][viz(j,1)]+spin[i][viz(j,-1)];

      if(sum*spin[i][j]<0)
        spin[i][j]*=-1;
      else{
        if(exp(-2*sum*spin[i][j]/T)>gsl_rng_uniform(y))
          spin[i][j]*=-1;
      }

      sum=0;
      
      flips+=1;
      E[2]=0;
      M[2]=0;
    }

    Emed=E[0]/((double)flips); //energia media
    Mmed=fabs(M[0])/((double)flips);//magnetização media
    pE=sqrt(E[3]/((double)flips)-pow(Emed,2));//desvio medio
    pM=sqrt(M[3]/((double)flips)-pow(Mmed,2));//desvio medio
      
    fprintf(EnMag,"%lf %lf %lf %lf %lf\n",T,Emed,pE,Mmed,pM);
    T+=dT;
    flips=1;

    //inicialização a zero
    for(i=0;i<4;++i){ 
      E[i]=0;
      M[i]=0;
    }

    n+=1;
  }
  
  free(spin);
  fclose(gnu);
  fclose(EnMag);
  return 0;
}


