//c_cryst.cpp
//config's member functions for making crystals

#include "config.h"

void config::MakeCryst(float To, short Ix, short Iy, short Iz, float rho,
		       short id, int seed, string cfgname){
  
  if (rho==0) rho=0.049852723;
  name=cfgname;
  NFixed=0;
  double Ta, Tscale;
  double sumvx=0; double sumvy=0; double sumvz=0;
  svector sumV;
  double X,Y,Z;
  int q=0;
  T=To;
  N=8*Ix*Iy*Iz; //number of real particles
  resize(N);
  mat=id;
  double a=pow(8.0/rho, 1.0/3.0); //lattice constant, typ 2.717 A.
  Lx=Ix*a; //box length
  Ly=Iy*a;
  Lz=Iz*a;
  
  for(i=begin; i<end; i++) i->id=mat;

  SetProps();
  
  //position init
  i=begin;
  
  if (mat==14)
    for (short x=0; x<Ix; x++)
      for (short y=0; y<Iy; y++)
	for (short z=0; z<Iz; z++){	
	  i->R.set((x+0.25)*a, (y+0.75)*a, (z+0.75)*a);
	  (i+1)->R.set((x+0.75)*a, (y+0.25)*a, (z+0.75)*a);
	  (i+2)->R.set(x*a, (y+0.5)*a, (z+0.5)*a);
	  (i+3)->R.set((x+0.5)*a, (y+1)*a, (z+0.5)*a);
	  (i+4)->R.set((x+0.25)*a, (y+0.25)*a, (z+0.25)*a);
	  (i+5)->R.set((x+0.75)*a, (y+0.75)*a, (z+0.25)*a);
	  (i+6)->R.set(x*a, y*a, z*a);     
	  (i+7)->R.set((x+0.5)*a, (y+0.5)*a, z*a);
	  if(z==0){
	    (i+4)->is_fixed=1; NFixed++;
	    (i+5)->is_fixed=1; NFixed++;
	    (i+6)->is_fixed=1; NFixed++;
	    (i+7)->is_fixed=1; NFixed++;
	  }
	  i+=8;
	}

  //origin is currently at bottom corner of crystal.
  //this puts it in the middle.
  //also use this loop to set the atoms' ix.
  
  for (i=begin; i<end; i++){
    i->R.x-=Lx/2;
    i->R.y-=Ly/2;
    i->R.z-=Lz/2;
    i->ix=q++;
    i->O=i->R;
  }
  
  //velocity init
  if (seed==0)
    srand((unsigned)clock()+(unsigned)time(NULL));
  else srand(seed);

  for (i=begin; i<end; i++){
    X=0; Y=0; Z=0;
    if (!i->is_fixed){
      for (q=0; q<12; q++) X+=(double)rand()/(double)RAND_MAX;
      X-=6.0;
      for (q=0; q<12; q++) Y+=(double)rand()/(double)RAND_MAX;
      Y-=6.0;
      for (q=0; q<12; q++) Z+=(double)rand()/(double)RAND_MAX;
      Z-=6.0;
    }
    i->V.set(X,Y,Z);
  }

  //scale to temperature
  for(i=begin; i<end; i++) if (!i->is_fixed) sumvx+=i->m*i->V.sqmag();
  Ta=sumvx/3.0/(double)(N-NFixed)/KB;
  Tscale=sqrt(T/Ta);
  for (i=begin; i<end; i++) i->V*=Tscale;
  sumvx=0;

  for (i=begin; i<end; i++){
    sumvx+=i->V.x; sumvy+=i->V.y; sumvz+=i->V.z;
  }

  sumV.set(sumvx/(double)N, sumvy/(double)N, sumvz/(double)N);

  //remove net momentum
  for (i=begin; i<end; i++) if(!i->is_fixed) i->V-=sumV;

  //finally, initialize time
  t=0;

  Nmax=N;
}
//**************************************************************************
void config::AddFixed(float rho){
  if (rho==0) rho=0.049852723;
  double a=pow(8.0/rho, 1.0/3.0);
  Lz+=a;
  particle dummy;
  float temp=T;
  set<int> fixed_list;
  int q; double X,Y,Z;
  for (i=begin; i<end; i++)
    if (i->is_fixed){
      i->is_fixed=0;
      X=Y=Z=0;
      for (q=0; q<12; q++) X+=(double)rand()/(double)RAND_MAX;
      X-=6.0;
      for (q=0; q<12; q++) Y+=(double)rand()/(double)RAND_MAX;
      Y-=6.0;
      for (q=0; q<12; q++) Z+=(double)rand()/(double)RAND_MAX;
      Z-=6.0;
      i->V.set(X,Y,Z);
      i->V*=sqrt(T/(i->Ek()/3/KB));
      fixed_list.insert(i->ix);
    }
  for (set<int>::iterator si=fixed_list.begin(); si!=fixed_list.end(); si++){  
    i=atomix(*si);
    dummy.R=i->R;
    dummy.R.z-=a/2;
    dummy.O=dummy.R;
    dummy.id=i->id;
    dummy.ix=Nmax++;
    dummy.is_fixed=1;
    dummy.SetProps();
    atom->push_back(dummy);
    begin=atom->begin();
    end=atom->end();
  }
  Partition();
  N+=NFixed;
  a=t;
  //Run(0.05,0,0,0,0,-1,1,0);
  //Thermo(temp, bhb_dt, 0, dt, 0, 0, 0);
  t=a;
}
	
    
  
