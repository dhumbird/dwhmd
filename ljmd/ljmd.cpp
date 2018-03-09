//ljmd.cpp
//molecular dynamics simulation code for a Lennard-Jones fluid
//default settings are for argon.

#include <fstream.h>
#include <iostream.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "vector.cpp"
#include "particle.cpp"
#include "toolbox.cpp"

/**********************************MAIN BODY******************************/
int main(int argc, char* argv[]){

  char cfgfile[256];

  if (argc!=3){
    cout<<"Usage: ljmd [input file] [time steps]"<<endl;
    exit(1);
  }

  strcpy(cfgfile, argv[1]);
  const int steps=atoi(argv[2]);

  /********variable declaration***************/
  const double kb=8.6173857e-5; //boltzmann
  const double pi=3.1415926559;
  int i,j,t;
  int N;
  const double sigma=3.405; //HS radius, ang
  const double sig6=pow(sigma,6);
  const double sig12=sig6*sig6;
  const double sig12_2=2*sig12;
  const double eps=kb*119.8; //L-J param, eV
  const double eps_4=4.0*eps;
  const double eps_24=24.0*eps;
  double L=0.0;
  double r=0;
  const double rc=2.5*sigma; //cutoff radius
  const double rlist=rc+0.5*sigma; //list cutoff
  const double rlistsq=rlist*rlist;
  const double rcsq=rc*rc;
  const double urc=(sig12/pow(rc,12)-sig6/pow(rc,6));
  double u=0.0; //system potential energy
  double k=0.0; //total kinetic energy
  double T=0.0; //reduced temperature
  const double h=0.001; //time increment
  const double hsq2=h*h/2.0;
  const double h2=h/2.0;  
  vector Rij, Fij;
  double rij=0.0;
  double rijsq=0.0;
  int atom_count;
  double rlower=0;
  double rupper=0;
  double nideal=0;

  /***********end variable declaration**********/

  /************read in atom cfg******************/
  
  ifstream fin;
  fin.open(cfgfile, ios::in);
  fin>>N>>T>>L;
  particle atom[N]; 
  for (i=0; i<N; i++){
    fin>>atom[i];
    //atom[i].listinit(N);
  }
  fin.close();
  const double rho=(double)N/L/L/L;
  //ReNeighbor(N, atom, rlistsq);
     
  double dr=0.1;
  int bin=0;
  int num_bins=(int)(L/dr/2)+1;
  int hist[num_bins];
  double g[num_bins];
  for (i=0; i<num_bins; i++){
    hist[i]=0;
    g[i]=0.0;
  }

  /***********begin time loop*****************/
  
  for (t=0; t<steps; t++){
    
    /******begin first part of velocity Verlet****/
    if (t>0){
      for (i=0; i<N; i++){
	atom[i].R+=h*atom[i].V+hsq2*atom[i].F/atom[i].m;
	atom[i].R.minimg(L);
	atom[i].V+=h2*atom[i].F/atom[i].m;
      }
    }
    /********end first part of VV***********/

    /********begin force routine**********/ 
  
    //re-init forces, potential
    for (i=0; i<N; i++) atom[i].F.clear();
    u=0.0;
    atom_count=0;
    for (i=0; i<N-1; i++){
      for (j=i+1; j<N; j++){
		
	//if(atom[i].nlist[j]=1){
	  //calculate interparticle distances
	  Rij=atom[i].R-atom[j].R;
	
	  //min image
	  Rij.minimg(L);
	  rij=Rij.mag();
	  bin=(int)(rij/dr);
	  if (bin<num_bins) hist[bin]+=2;
	  if (rij<rc){
	    atom_count++;
	    //rij=sqrt(rijsq);
	    
	    //evaluate forces
	    Fij=(sig12_2/pow(rij,13.0)-sig6/pow(rij,7.0))*Rij/rij;
	    
	    atom[i].F+=Fij;
	    atom[j].F-=Fij;
	    
	    //compute potential energy
	    u+=(sig12/pow(rij,12)-sig6/pow(rij,6));
	  }
	  //}
      }
    }
    for(i=0; i<N; i++) atom[i].F*=eps_24;
    u-=(double)atom_count*urc;
    u*=eps_4;
    /**************end force routine*********/
    
    /******second part of velocity verlet******/

    //reinit k
    k=0.0;
    for (i=0; i<N; i++) atom[i].V+=h2*atom[i].F/atom[i].m;
    
    /*********end 2nd VV***********************/

    //compute kinetic energy    
    for (i=0; i<N; i++) k+=atom[i].m*atom[i].V.sqmag(); 
    
    k/=2.0;
    
    //printf("\n%d\t%e\t%e\t%e", t, u, k, u+k);
    //if (t % 20 == 0) ReNeighbor(N, atom, rlistsq);
  
  } 
  
  /************end time loop***************/	       
  //cout<<"\n";
  double c=4.0*pi*rho/3.0;
  
  for (i=0; i<num_bins; i++){
    rlower=(double)i*dr;
    rupper=rlower+dr;
    nideal=c*(pow(rupper,3)-pow(rlower,3));
    g[i]=(double)hist[i]/(double)steps/(double)N/nideal;
    cout<<(double)i*dr<<"\t"<<g[i]<<endl;
  }
  return(0);
};








