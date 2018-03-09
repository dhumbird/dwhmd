#include <fstream.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "vector.cpp"
#include "particle.cpp"

int main(int argc, char *argv[]){

  //***********give default args****************
  double T=0;
  int I=3;
  double rho=0.01166;
  double mass=39.948/9648.531; //reduced mass (argon)
  char outfile[256]="temp.cfg";

  //***********test input***********************

  int i;

  for (i=1; i<argc; i++){
    if (!strcmp(argv[i],"-T")) T=atof(argv[++i]);
    else if (!strcmp(argv[i], "-I")) I=atoi(argv[++i]);
    else if (!strcmp(argv[i], "-rho")) rho=atof(argv[++i]);
    else if (!strcmp(argv[i], "-of")) strcpy(outfile, argv[++i]);
  }
  
  ofstream fout;

  fout.open(outfile, ios::out);
  fout.fill(' ');
  fout.precision(24);
  fout.setf(ios::fixed);
  
  const int N=4*(int)pow(I,3); //number of real particles
  const double a=pow(4.0/rho, 1.0/3.0); //lattice constant
  const double L=I*a; //box length

  //************make some atoms*******************  

  particle atom[N]; 
  for (i=0; i<N; i++){
    atom[i].F.clear();  
    atom[i].R.clear();
    atom[i].V.clear();
    atom[i].listinit(N);
    atom[i].m=mass;
  }

  //*****************position init**************

  i=0;
  for (int x=0; x<I; x++){
    for (int y=0; y<I; y++){
      for (int z=0; z<I; z++){	
	atom[i++].R.set(x*a, y*a, z*a);     
	atom[i++].R.set((x+0.5)*a, y*a, (z+0.5)*a);
	atom[i++].R.set(x*a, (y+0.5)*a, (z+0.5)*a);
	atom[i++].R.set((x+0.5)*a, (y+0.5)*a, z*a);
      }    
    }
  }

  for (i=0; i<N; i++) atom[i].R-=L/2;

  //**************velocity init*****************

  int j;
  double Ta, Tscale;
  double sumvx, sumvy, sumvz, sumv;
  vector sumV;
  sumvx=0.0;
  sumvy=0.0;
  sumvz=0.0;
  sumv=0.0;
  sumV.clear();
  const double kb=8.6173857e-5; //boltzmann
  double X,Y,Z;

  srand((unsigned)clock()+(unsigned)time(NULL));

  for (j=0; j<N; j++){
    X=0; Y=0; Z=0;
    for (i=0; i<12; i++) X+=(double)rand()/(double)RAND_MAX;
    X-=6.0;
    for (i=0; i<12; i++) Y+=(double)rand()/(double)RAND_MAX;
    Y-=6.0;
    for (i=0; i<12; i++) Z+=(double)rand()/(double)RAND_MAX;
    Z-=6.0;
    
    atom[j].V.set(X,Y,Z);
  }
  
  //scale to temperature
  for(i=0; i<N; i++) sumv+=atom[i].m*atom[i].V.sqmag();
  Ta=sumv/3.0/(double)N/kb;
  Tscale=sqrt(T/Ta);
  for (i=0; i<N; i++){
    atom[i].V*=Tscale;
  } 

  for (i=0; i<N; i++){
    sumvx+=atom[i].V.X();
    sumvy+=atom[i].V.Y();
    sumvz+=atom[i].V.Z();
  }
  
  sumV.set(sumvx/(double)N, sumvy/(double)N, sumvz/(double)N);
  
  //remove net momentum
  for (i=0; i<N; i++) atom[i].V-=sumV;

  //****************write file*******************
  fout<<N<<" "<<T<<" "<<L<<endl;
  for (i=0; i<N; i++){
    fout<<atom[i]<<endl;
  }

  fout.close();

  return(0);
};
