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
  
  
  VPI ion=append(id); //ion points to the ion
  if (id==18 || id==2 || id==10 || id==36) inert=ion->ix;
  cerr<<"* Adding type "<<id<<" ion. (# "<<ion->ix<<") Seed: "<<seed<<endl;

  //assign ion x, y position
  ion->R.set(rand01()*fabs(Lx), rand01()*fabs(Ly), 0);
  ion->R.minimg(Lx, Ly);
  
  //assign ion's minimum z
  ion->R.z=MaxZ()+sqrt(ion->rc)+0.001;
  Partition();
  //pick incident & azimuthal angle
  if (inc!=0) inc=inc/180*PI;
  if (choose_azi) azi=rand01()*2*PI;
  if (choose_inc) inc=rand01()*0.5*PI;
  //assign random velocity
  do{
    vx=vy=vz=0;
    double kT=sqrt(KB*To/ion->m);
    for (int q=0; q<12; q++){
      vx+=rand01();
      vy+=rand01();
      vz+=rand01();
    }
    vx-=6.0; vx*=kT; 
    vy-=6.0; vy*=kT;
    vz-=6.0; vz*=kT;
    svector Vth(vx,vy,vz);
    svector Vsh;
    //add sheath effects or if e==0 use half-maxwellian directed down
    if (e==0){
      Vth.z=-fabs(Vth.z);
      ion->V=Vth;
    }
    else{
      speed=sqrt(2/ion->m*e);
      Vsh.set(speed*sin(inc)*cos(azi),
	      speed*sin(inc)*sin(azi),
	      -speed*cos(inc));
      ion->V.x=pow(Vth.x,2)*signof(Vth.x)+pow(Vsh.x,2)*signof(Vsh.x);
      ion->V.x=ion->V.x > 0 ? sqrt(ion->V.x) : -sqrt(-1*ion->V.x);
      ion->V.y=pow(Vth.y,2)*signof(Vth.y)+pow(Vsh.y,2)*signof(Vsh.y);
      ion->V.y=ion->V.y > 0 ? sqrt(ion->V.y) : -sqrt(-1*ion->V.y);
      ion->V.z=pow(Vth.z,2)*signof(Vth.z)+pow(Vsh.z,2)*signof(Vsh.z);
      ion->V.z=ion->V.z > 0 ? sqrt(ion->V.z) : -sqrt(-1*ion->V.z);
    }
    ion->is_fixed=0;
    //ionek=ion->Ek();
  }while((ion->R.z - the_top)/vz > 10);
  svector AV=ion->V;
  //take down to first neighbor
  svector dir=0.01*AV/AV.mag();
  bool got_neighbors=0;
  while (!got_neighbors){
    subcell* sc=WhichCell(ion->R);
    for (c=sc->nbegin(); c!=sc->nend(); c++)
      for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	if (!got_neighbors)
	  if (&(*ion)!=(*j0))
	    if ((ion->R - (*j0)->R).minsqmag(Lx, Ly) < (ion->rc >? (*j0)->rc))
	      got_neighbors=1;
    if (!got_neighbors){
      ion->R+=dir;
      ion->R.minimg(Lx,Ly);
      if ((this_cell=WhichCell(ion->R)) != ion->my_cell){
	((subcell*) ion->my_cell)->erase(&(*ion));
	((subcell*) this_cell)->insert(&(*ion));
      }
    }
  }
  ion->R-=dir; ion->R.minimg(Lx,Ly);
  Partition();
  ReNeighbor();
  
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
