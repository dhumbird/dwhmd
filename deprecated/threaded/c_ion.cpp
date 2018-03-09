//c_ion.cpp
//config's member functions for adding/deleting ions

#include "config.h"

//*************************************************************************
void config::AddIon(float e, short id, float To, bool choose_inc, float inc,
		    bool choose_azi, float azi, int seed, bool quiet,
		    bool no_noise, bool choose_pos, bool confine, bool round,
		    float conf){
  double speed=0;
  double vx=0;
  double vy=0;
  double vz=0;
  double zmax=0;

  if (seed==0) seed=(unsigned)clock()+(unsigned)time(NULL);
  srand(seed);

  //make new atom
  resize(++N);
  j=end-1; //j points to the ion
  j->id=id;
  
  j->SetProps();
  
  //assign ion position
  if (choose_pos){
    if (!confine)
      j->R.set(rand01()*fabs(Lx), rand01()*fabs(Ly), MaxZ()+sqrt(j->rc)+0.001);
    else{
      double rndx=rand01()*conf*2-conf;
      double rndy=rand01()*conf*2-conf;
      if (round){
	while (sqrt(pow(rndx,2)+pow(rndy,2))>conf){
	  rndx=rand01()*conf*2-conf;
	  rndy=rand01()*conf*2-conf;
	}
      }
      j->R.set(rndx, rndy, MaxZ()+sqrt(j->rc)+0.001);
    }
    j->R.minimg(Lx, Ly, Lz);
  }
  else j->R.set(0, 0, MaxZ()+sqrt(j->rc)+0.001);
  
  //pick incident & azimuthal angles
  if (choose_azi) azi=rand01()*2*PI;
  if (choose_inc) inc=rand01()*0.9444*PI;

  //assign random velocity
  if (!no_noise){
    int q;
    for (q=0; q<12; q++) vx+=rand01();
    vx-=6.0;
    for (q=0; q<12; q++) vy+=rand01();
    vy-=6.0;
    for (q=0; q<12; q++) vz+=rand01();
    vz-=6.0;

    //scale to temperature
    double Tscale=sqrt(To/(j->m*(vx*vx+vy*vy+vz*vz)/3.0/KB));
    vx*=Tscale; vy*=Tscale; vz*=Tscale;
  }

  //add sheath effects
  speed=sqrt(2/j->m*e);
  vx+=speed*sin(inc)*cos(azi);
  vy+=speed*sin(inc)*sin(azi);
  vz-=speed*cos(inc);

  j->V.set(vx, vy, vz);
  j->is_fixed=0; 
  
  if (!quiet){
    cerr<<"* Type "<<id<<" ion added at "<<j->R<<". Total atoms: "<<N<<endl;
    cerr<<"  Seed: "<<seed<<". Energy: "<<j->Ek()<<" eV"<<endl;
    cerr<<"  Incidence: "<<inc*PI/180<<" deg; Azimuth: "<<azi*PI/180<<" deg"<<endl;
  }
  RePartition();
}
//*************************************************************************
void config::DelIon(bool quiet){  
  resize(--N);
  if (!quiet) cerr<<"* Ion deleted. Total number of atoms in cfg: "<<N<<endl;
}
