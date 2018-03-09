//c_ion.cpp
//config's member functions for adding/deleting ions

#include "config.h"

//*************************************************************************
VPI config::AddIon(float e, int id, float To, bool choose_inc, float inc,
		    bool choose_azi, float azi, int seed, bool quiet){
  double speed=0;
  double vx=0;
  double vy=0;
  double vz=0;

  if (seed!=-1){
    if (seed==0) seed=TimeSeed();
    srand(seed);
  }
  cerr<<"* Adding type "<<id<<" ion. Seed: "<<seed<<endl;
  
  VPI ion=append(id); //ion points to the ion
  if (id==18 || id==2 || id==10 || id==36) inert=ion->ix;

  //assign ion x, y position
  ion->R.set(rand01()*fabs(Lx), rand01()*fabs(Ly), 0);
  ion->R.minimg(Lx, Ly);
  
  //assign ion's minimum z
  ion->R.z=MaxZ()+sqrt(ion->rc)+0.001;
  Partition();
  bool got_neighbors=0;
  while (!got_neighbors){
    subcell* sc=WhichCell(ion->R);
    for (c=sc->nbegin(); c!=sc->nend(); c++)
      for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	if (!got_neighbors)
	  if (&(*ion)!=(*j0))
	    if ((ion->R - (*j0)->R).minsqmag(Lx, Ly) < (ion->rc >? (*j0)->rc))
	      got_neighbors=1;
    if (!got_neighbors) ion->R.z-=0.01;
  }
  ion->R.z+=0.01;
  Partition();
  //pick incident & azimuthal angle
  if (inc!=0) inc=inc/180*PI;
  if (choose_azi) azi=rand01()*2*PI;
  if (choose_inc) inc=0.9444*PI;
  
  //assign random velocity
  int q;
  for (q=0; q<12; q++) vx+=rand01();
  vx-=6.0;
  for (q=0; q<12; q++) vy+=rand01();
  vy-=6.0;
  for (q=0; q<12; q++) vz+=rand01();
  vz-=6.0;
  
  //scale to temperature
  double Tscale=sqrt(To/(ion->m*(vx*vx+vy*vy+vz*vz)/3.0/KB));
  vx*=Tscale; vy*=Tscale; vz*=Tscale;
  
  //add sheath effects
  speed=sqrt(2/ion->m*e);
  vx+=speed*sin(inc)*cos(azi);
  vy+=speed*sin(inc)*sin(azi);
  vz-=speed*cos(inc);
  
  ion->V.set(vx, vy, vz);
  ion->is_fixed=0; 
  
  if (!quiet){
    cerr<<"* Ion added at "<<ion->R<<". Total atoms: "<<N<<endl;
    cerr<<"  Seed: "<<seed<<". Energy: "<<ion->Ek()<<" eV"<<endl;
    cerr<<"  Incidence: "<<inc/PI*180<<" deg; Azimuth: "<<azi/PI*180<<" deg"<<endl;
  }
  return ion;
}
//*************************************************************************
void config::DelIon(bool quiet){  
  if (inert!=-1){
    erase(atomix(inert));
    if (!quiet) cerr<<"* Ion deleted. Total number of atoms in cfg: "<<N<<endl;
    Partition();
    ReNeighbor();
    inert=-1;
  }
  else
    if (!quiet) cerr<<"* No ion to delete. Total number of atoms in cfg: "
		    <<N<<endl;
}
