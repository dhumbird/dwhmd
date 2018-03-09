//c_ion.cpp
//config's member functions for adding/deleting ions

#include "config.h"
extern map<short,double> RC_SQ;
//*************************************************************************
VAI config::AddIon(float e, string species, float To, bool choose_inc, 
		   float inc, bool choose_azi, float azi, int seed, 
		   bool quiet, float confine){
  string type;
  double speed=0;
  double vx=0;
  double vy=0;
  double vz=0;
  double ionek=0;
  float the_top=MaxZ();
 
  if (seed==0) seed=TimeSeed();
  srand(seed);
  svector AV;

  cerr<<"* Adding species "<<species<<". (# "<<Nmax<<") Seed: "<<seed<<endl;
  VAI ion;
  svector Vcom;
  double rcsq;
  if (species=="F" || species=="Ar" || species=="C" || species=="Cl"){
    if (species=="Ar"){
      rcsq=36;
      ion=append(18);
      inert=ion->ix;
    }
    else if (species=="F") ion=append(9);
    else if (species=="Cl") ion=append(17);
    else if (species=="C") ion=append(6);
    //assign ion x, y position
    if (confine){
      double rndx=rand01()*confine*2-confine;
      double rndy=rand01()*Ly;
      //double rndy=rand01()*confine*2-confine;
      //while (sqrt(pow(rndx,2)+pow(rndy,2))>confine){
      //	rndx=rand01()*confine*2-confine;
      //	rndy=rand01()*confine*2-confine;
      //}
      ion->R.set(rndx, rndy, 0); 
      ion->R.minimg(Lx, Ly);
    }
    else{
      ion->R.set(rand01()*Lx, rand01()*Ly, 0);
      ion->R.minimg(Lx, Ly);
    }
    //assign initially high z
    ion->R.z=the_top+6.01;
    Partition();
    //pick incident & azimuthal angle
    if (inc!=0) inc=inc/180*PI;
    if (azi!=0) azi=azi/180*PI;
    if (choose_azi) azi=rand01()*2*PI;
    if (choose_inc) inc=rand01()*PI/2;
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
      ionek=ion->Ek();
    }while((ion->R.z - the_top)/vz > 10);
    AV=ion->V;
    //take down to first neighbor
    svector dir=0.01*AV/AV.mag();
    bool got_neighbors=0;
    while (!got_neighbors){
      for (VAI i=begin; i<ion; i++){
	if (inert==ion->ix) rcsq=36;
	else rcsq=RC_SQ[i->id+ion->id];
	if ((ion->R - i->R).minsqmag(Lx, Ly) < rcsq){
	  got_neighbors=1;
	}
      }
      if (!got_neighbors){
	ion->R+=dir;
	//cerr<<ion->R<<endl;
	ion->R.minimg(Lx,Ly);
	if ((this_cell=WhichCell(ion->R)) != ion->my_cell){
	  ((subcell*) ion->my_cell)->erase(&(*ion));
	  ((subcell*) this_cell)->insert(&(*ion));
	}
      }
    } 
//     ion->R+=50*dir;   ion->R.minimg(Lx,Ly);
//     Partition();
//     ReNeighbor(); TimeStepInit();
//     ForceEval();
//     cerr<<"* immediate force "<<ion->F<<endl;
//     ion->R-=50*dir;
    ion->R-=dir;
    ion->R.minimg(Lx,Ly);
    Partition();
    ReNeighbor(); 
  }
  else{ 
    //VAI ion will point to the carbon
    //assign ion x, y center of mass position
    double ionmass=0;
    if (species=="CF3"){
      ion=append(6);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(0, 1.267, 0);
      VAI H2=append(9);
      H2->R.set(-1.097, -0.6336, 0);
      VAI H3=append(9);
      H3->R.set(1.097, -0.6336, 0);
      ion=end-4;
    }
    else if (species=="CF2"){
      ion=append(6);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(-1.267, 0, 0);
      VAI H2=append(9);
      H2->R.set(1.267, 0, 0);
      ion=end-3;
    }
    else if (species=="CF"){
      ion=append(6);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(-1.267, 0, 0);
      ion=end-2;
    }
    else if (species=="SiF3"){
      ion=append(14);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(0, 1.59, 0);
      VAI H2=append(9);
      H2->R.set(-1.377, -0.795, 0);
      VAI H3=append(9);
      H3->R.set(1.377, -0.795, 0);
      ion=end-4;
    }
    else if (species=="SiF2"){
      ion=append(14);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(0, 1.59, 0);
      VAI H2=append(9);
      H2->R.set(0, -1.59, 0);
      ion=end-3;
    }
    else if (species=="SiF"){
      ion=append(14);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(0, 1.59, 0);
      ion=end-2;
    }
    else if (species=="F2"){
      ion=append(9);
      ion->R.clear();
      VAI H1=append(9);
      H1->R.set(0, 1.435, 0);
      ion=end-2;
    }
    else if (species=="Cl2"){
      ion=append(17);
      ion->R.clear();
      VAI H1=append(17);
      H1->R.set(0, 2.02, 0);
      ion=end-2;
    }
    //angle of orientation
    double phi=rand01()*2*PI;
    double theta=rand01()*2*PI;
    double psi=rand01()*2*PI;
    svector R;
    if (confine){
      double rndx=rand01()*confine*2-confine;
      double rndy=rand01()*confine*2-confine;
      while (sqrt(pow(rndx,2)+pow(rndy,2))>confine){
	rndx=rand01()*confine*2-confine;
	rndy=rand01()*confine*2-confine;
      }
      R.set(rndx, rndy, 0);
    }
    else R=svector(rand01()*Lx, rand01()*Ly, 0);
    for (VAI i=ion; i<end; i++){
      ionmass+=i->m;
      i->R.EulerTrans(phi, theta, psi);
      i->R+=R;
      i->R.minimg(Lx,Ly);
      i->R.z+=(the_top+5);
    }
    Partition();
    //pick incident & azimuthal angles
    if (inc!=0) inc=inc/180*PI;
    if (choose_azi) azi=rand01()*2*PI;
    if (choose_inc) inc=rand01()*PI/2;
    
    double netvz;
    double kT=sqrt(KB*To/ionmass);
    
    //assign random thermal velocity
    do{
      netvz=0;
      Vcom.clear();
      ionek=0;
      for (VAI i=ion; i<end; i++){
	vx=vy=vz=0;
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
	//add sheath effects
	speed=sqrt(2/ionmass*e);
	Vsh.set(speed*sin(inc)*cos(azi),
		speed*sin(inc)*sin(azi),
		-speed*cos(inc));
	i->V.x=pow(Vth.x,2)*signof(Vth.x)+pow(Vsh.x,2)*signof(Vsh.x);
	i->V.x=i->V.x > 0 ? sqrt(i->V.x) : -sqrt(-1*i->V.x);
	i->V.y=pow(Vth.y,2)*signof(Vth.y)+pow(Vsh.y,2)*signof(Vsh.y);
	i->V.y=i->V.y > 0 ? sqrt(i->V.y) : -sqrt(-1*i->V.y);
	i->V.z=pow(Vth.z,2)*signof(Vth.z)+pow(Vsh.z,2)*signof(Vsh.z);
	i->V.z=i->V.z > 0 ? sqrt(i->V.z) : -sqrt(-1*i->V.z);
	i->is_fixed=0;
	ionek+=i->Ek();
	//netvz+=i->m*vz;
	Vcom+=i->m*i->V;
      }
      Vcom/=ionmass;
    }while(Vcom.z > 0);
    AV=Vcom;
    //take down to first neighbor
    svector dir=0.01*AV/AV.mag();
    bool got_neighbors=0;
    while (!got_neighbors){
      for (VAI i=ion; i<end; i++){
	for (VAI j=begin; j<ion; j++)
	  if (!got_neighbors){
	    if ((i->R - j->R).minsqmag(Lx, Ly) < RC_SQ[i->id+j->id])
	      got_neighbors=1;
	  }
      }
      if (!got_neighbors)
	for (VAI i=ion; i<end; i++){
	  i->R+=dir; i->R.minimg(Lx,Ly);
	  if ((this_cell=WhichCell(i->R)) != i->my_cell){
	    ((subcell*) i->my_cell)->erase(&(*i));
	    ((subcell*) this_cell)->insert(&(*i));
	  }
	}
    }
    for (VAI i=ion; i<end; i++){
      i->R-=dir;
      i->R.minimg(Lx,Ly);
    }
    Partition();
    ReNeighbor();
  }
  if (!quiet){
    AV.Rec2Sph();
    cerr<<"* Species added at "<<ion->R<<". Total atoms: "<<N<<endl;
    cerr<<"  Energy: "<<ionek<<" eV; ";
    cerr<<"Actual Incidence: "<<180-AV.z/PI*180<<
      " deg; Azimuth: "<<AV.y/PI*180<<" deg"<<endl;
    if (confine) cerr<<"  Confined to "<<confine<<" Ang around center\n";
  }
  return ion;
}
//*************************************************************************
void config::DelIon(bool quiet){  
  if (inert!=-1){
    erase(atomix(inert));
    inert=-1;
    if (!quiet) cerr<<"* Ion deleted. Total number of atoms in cfg: "<<N<<endl;
    Partition();
    ReNeighbor();
  }
  else
    if (!quiet) cerr<<"* No ion to delete. Total number of atoms in cfg: "
		    <<N<<endl;
}
