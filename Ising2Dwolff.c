#include <stdio.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <unistd.h>
#include <time.h>
#include <stdlib.h>

/*
 * compilação: gcc Ising2Dwolff.c -o "..." -lgsl -lgslcblas -lm
 * gera ficheiro "Wolffenmag.dat"
 */

#define N 64 //grelha NxN
#define Flmax 5000 //numero maximo de passos Wolff
#define DT 0.1

// Hamiltoniano: H=-JS(i)S(i+1)

typedef struct{
  double w[2];//w[0]=spin {-1,1}; w[1]=agre(1) ou n agreg(0)
}sdat;

typedef struct{
  int c[2]; // c[0]=1ª coordenada  (i); c[1]=2ª coordenada (j)
}nb; // neighbour coordinates;

static gsl_rng *rnduniform(int); //gera numeros aleatorios
static nb *viz(int ,int ); // determina os quatro vizinhos de um determinado spin
void cluster(sdat *,double ,gsl_rng *,gsl_rng *); //construção do cluster
double Mag(sdat *); // determina a magnetização
double En(sdat *); // determina a energia

int main(){
  /*
   * clusterflip=número de clusters invertidos
   * n= num de agregados fronteira
   * spin=array[i+N*j]; toma em (i,j) sdat[0,1]
   * T: temperatura;
   * Tmax: temperatura máxima
   * E:energia
   * M:magnetização
   * DE: desvio para a energia
   * DM: desvio para a magnetização
   * Eq: E²
   * Mq: M²
   */
  
  int i,j,clusterflip=0; 
  double T,Tmax,m,e,E=0,Eq=0,DE,Mq=0,M=0,DM;
  sdat *spin;
  gsl_rng *x,*y;
  FILE *mag;

  printf("Ti=");// temperatura inicial
  scanf("%lf",&T);
  printf("Tf=");//temperatura final
  scanf("%lf",&Tmax);
  
  mag=fopen("Wolffenmag.dat","w");//ficheiro pnde é guardada a informação sobre T,E,DE,DM
  spin=malloc(N*N*sizeof(sdat));
  
  // inicialização dos spins
  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      spin[i+N*j].w[0]=1; //spin +1
      spin[i+N*j].w[1]=0; // n agregado
    }
  }

  //gerar numeros aleatorios
  x=rnduniform(100);
  y=rnduniform(1000);
  
  while(T<Tmax){
    
    while(clusterflip<Flmax){
      m=Mag(spin);
      M+=fabs(m); //magnetização
      Mq+=pow(m,2);
      e=En(spin);
      E+=e;
      Eq+=pow(e,2);
      
      //construção do cluster
      cluster(spin,T,y,x);
      clusterflip+=1;
      
      for(i=0;i<N;++i){
        for(j=0;j<N;++j){
          spin[i+N*j].w[1]=0; // iniciar a 0 a agregação dos spins
        }
      }
      
    }

    E=E/(double)clusterflip;//energia média
    Eq=Eq/(double)clusterflip; // média de E^2
    DE=sqrt(Eq-pow(E,2));
    M=M/(double)clusterflip; // magnetização média
    Mq=Mq/(double)clusterflip; // média de M^2
    DM=sqrt(Mq-pow(M,2));
    
    fprintf(mag,"%lf %lf %lf %lf %lf\n",T,E,DE,M,DM);
    clusterflip=0;
    T+=DT;
    M=0;
    E=0;
    Mq=0;
    Eq=0;
  }
  
  free(spin);
  fclose(mag);
  
  return 0;
}

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

static nb *viz(int i,int j){ //coordenada do vizinho (i,j)
  static nb v[4]; // 4 vizinhos
  
  if(j==N-1 || j==0){ //condições periodicas
    v[0].c[0]=i;
    v[0].c[1]=N-j-1;
    v[1].c[0]=i;
    v[1].c[1]=abs(j-1);
  }
  else{
    v[0].c[0]=i;
    v[0].c[1]=j+1;
    v[1].c[0]=i;
    v[1].c[1]=j-1;
  }

  if(i==N-1 || i==0){ //condições periodicas
    v[2].c[0]=N-i-1;
    v[2].c[1]=j;
    v[3].c[0]=abs(i-1);
    v[3].c[1]=j;
  }
  else{
    v[2].c[0]=i+1;
    v[2].c[1]=j;
    v[3].c[0]=i-1;
    v[3].c[1]=j;
  }
 
  return v;
}

void cluster(sdat *spin,double T,gsl_rng *y,gsl_rng *x){

  /* Começamos com um array q contem o spin semente e visitamos os spins vizinhos. Agregamo-los de acordo com as regras de Wolff. Os que são agragados passam a
   * ter numero de agregação (spin[].w[1]=>1)igual a 1.Procedemos da mesma forma para os novos spins agregados
   * fg: array q contem os spins fronteira
   * spin: array spin
   * s: orientação do spin semente
   * T: temperatura
   */
  
  int i,j,a,b,n,m,M=1,Nagreg=0,s;
  nb *tmp,*fg;

  tmp=malloc(1*sizeof(nb)); //array q vai armazenar os spins agregados
  fg=malloc(1*sizeof(nb));
  
  fg[0].c[0]=i=gsl_rng_uniform(x)*N;
  fg[0].c[1]=j=gsl_rng_uniform(x)*N;
  s=spin[i+N*j].w[0];
  spin[i+N*j].w[1]=1; // obviamente este spin pertence ao agregado
  
  while(M!=0){
    
    for(m=0;m<M;++m){
      i=fg[m].c[0];
      j=fg[m].c[1];
      for(n=0;n<4;++n){
        a=viz(i,j)[n].c[0];
        b=viz(i,j)[n].c[1];
        if(spin[a+N*b].w[0]==s && spin[a+N*b].w[1]==0 && gsl_rng_uniform(y)<(1-exp(-2.0/T))){
          spin[a+N*b].w[0]*=-1; // flip spin
          spin[a+N*b].w[1]=1; // spin é agregado
          
          tmp=realloc(tmp,(Nagreg+1)*sizeof(nb));
          tmp[Nagreg].c[0]=a;
          tmp[Nagreg].c[1]=b;
          Nagreg+=1;//aumenta +1 o numero de spins ja agregados
        }
      }
      
    }
    
    fg=realloc(fg,Nagreg*sizeof(nb));
    for(m=0;m<Nagreg;++m){
      fg[m]=tmp[m];
    }
    
    M=Nagreg;
    Nagreg=0;
  }
  free(tmp);
  free(fg);
  
}

double Mag(sdat *spin){
  int i,j;
  double m=0;

  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      m+=spin[i+N*j].w[0]/(N*N);
    }
  }
  return m;
}

double En(sdat *spin){
  int i,j,a,b,n;
  double m=0;

  for(i=0;i<N;++i){
    for(j=0;j<N;++j){
      for(n=0;n<4;++n){
        a=viz(i,j)[n].c[0];
        b=viz(i,j)[n].c[1];
        m+=-spin[a+N*b].w[0]*spin[i+N*j].w[0]/(N*N);
      }
    }
  }
  
  return m/2;
}


