//c_cryst.cpp
//config's member functions for making crystals

#include "config.h"

void config::MakeCryst(float To, short Ix, short Iy, short Iz, float rho,
		       short id, int seed, string cfgname){
  
  if (rho==0){
    //if (id==6) rho=0.1778; //REBO density
    if (id==6) rho = 0.18637; //from CFA code
    else if (id==14) rho = 0.049853; //from CFA code
  }
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
  double a=pow(8.0/rho, 1.0/3.0); //lattice constant, typ 2.717 A.
  Lx=Ix*a; //box length
  Ly=Iy*a;
  Lz=Iz*a;
  
  for(VAI i=begin; i<end; i++) i->id=id;

  SetProps();
  
  //position init
  VAI s=begin;
  
  for (short x=0; x<Ix; x++)
    for (short y=0; y<Iy; y++)
      for (short z=0; z<Iz; z++){	
	s->R.set((x+0.25)*a, (y+0.75)*a, (z+0.75)*a);
	(s+1)->R.set((x+0.75)*a, (y+0.25)*a, (z+0.75)*a);
	(s+2)->R.set(x*a, (y+0.5)*a, (z+0.5)*a);
	(s+3)->R.set((x+0.5)*a, (y+1)*a, (z+0.5)*a);
	(s+4)->R.set((x+0.25)*a, (y+0.25)*a, (z+0.25)*a);
	(s+5)->R.set((x+0.75)*a, (y+0.75)*a, (z+0.25)*a);
	(s+6)->R.set(x*a, y*a, z*a);     
	(s+7)->R.set((x+0.5)*a, (y+0.5)*a, z*a);
	if(z==0){
	  //i->is_fixed=1; NFixed++;
	  //(i+1)->is_fixed=1; NFixed++;
	  //(i+2)->is_fixed=1; NFixed++;
	  //(i+3)->is_fixed=1; NFixed++;
	  (s+4)->is_fixed=1; NFixed++;
	  (s+5)->is_fixed=1; NFixed++;
	  (s+6)->is_fixed=1; NFixed++;
	  (s+7)->is_fixed=1; NFixed++;
	}
	s+=8;
      }
  
  //origin is currently at bottom corner of crystal.
  //this puts it in the middle.
  //also use this loop to set the atoms' ix.
  
  for (VAI i=begin; i<end; i++){
    i->R.x-=Lx/2;
    i->R.y-=Ly/2;
    i->R.z-=Lz/2;
    i->ix=q++;
  }
  
  //velocity init
  if (seed==0)
    srand((unsigned)clock()+(unsigned)time(NULL));
  else srand(seed);

  for (VAI i=begin; i<end; i++){
    X=Y=Z=0;
    if (!i->is_fixed){
      double kT=sqrt(KB*To/i->m);
      for (q=0; q<12; q++){
	X+=rand01();
	Y+=rand01();
	Z+=rand01();
      }
      X-=6.0; X*=kT;
      Y-=6.0; Y*=kT;
      Z-=6.0; Z*=kT;
    }
    i->V.set(X,Y,Z);
  }

  for (VAI i=begin; i<end; i++){
    sumvx+=i->V.x; sumvy+=i->V.y; sumvz+=i->V.z;
  }
  //remove net momentum
  double sumek=0;
  sumV.set(sumvx/(double)N, sumvy/(double)N, sumvz/(double)N);
  for (VAI i=begin; i<end; i++) 
    if(!i->is_fixed){
      i->V-=sumV;
      sumek+=i->Ek();
    }
  T=sumek/3/(double)(N-NFixed)/KB;
  
  //finally, initialize time
  t=0;

  Nmax=N;
}
//************************************************************
void config::Maxwell(float To, int seed){
  
  double Ta, Tscale;
  double sumvx=0; double sumvy=0; double sumvz=0;
  svector sumV;
  double X,Y,Z;
  int q=0;
  T=To;
 
  //velocity init
  if (seed==0)
    srand((unsigned)clock()+(unsigned)time(NULL));
  else srand(seed);

  for (VAI i=begin; i<end; i++){
    X=0; Y=0; Z=0;
    if (!i->is_fixed){
      double kT=sqrt(KB*To/i->m);
      for (q=0; q<12; q++){
	X+=rand01();
	Y+=rand01();
	Z+=rand01();
      }
      X-=6.0; X*=kT;
      Y-=6.0; Y*=kT;
      Z-=6.0; Z*=kT;
    }
    i->V.set(X,Y,Z);
  }
  //remove net momentum
  sumvx=0;
  for (VAI i=begin; i<end; i++){
    sumvx+=i->V.x; sumvy+=i->V.y; sumvz+=i->V.z;
  }
  double sumek=0;
  sumV.set(sumvx/(double)N, sumvy/(double)N, sumvz/(double)N);
  for (VAI i=begin; i<end; i++) 
    if(!i->is_fixed){
      i->V-=sumV;
      sumek+=i->Ek();
    }
  T=sumek/3/(double)(N-NFixed)/KB;
}
//**************************************************************************
void config::AddFixed(float rho){
  cerr<<"* Adding new bottom atoms.\n";
  if (rho==0) rho=0.049853;
  double a=pow(8.0/rho, 1.0/3.0);
  Lz+=a/2;
  atom dummy;
  float temp=T;
  set<int> fixed_list;
  int q; double X,Y,Z;
  for (VAI i=begin; i<end; i++)
    if (i->is_fixed){
      double kT=sqrt(KB*T/i->m);
      i->is_fixed=0;
      X=Y=Z=0;
      for (q=0; q<12; q++){
	X+=rand01();
	Y+=rand01();
	Z+=rand01();
      }
      X-=6.0; X*=kT;
      Y-=6.0; Y*=kT;
      Z-=6.0; Z*=kT;
      i->V.set(X,Y,Z);
      fixed_list.insert(i->ix);
    }
  for (set<int>::iterator si=fixed_list.begin(); si!=fixed_list.end(); si++){
    VAI i=atomix(*si);
    dummy.R=i->R;
    dummy.R.x+=a;
    dummy.R.y-=1.5*a;
    dummy.R.z-=0.5*a;
    dummy.R.minimg(Lx,Ly);
    dummy.id=i->id;
    dummy.ix=Nmax++;
    dummy.is_fixed=1;
    dummy.SetProps();
    atom_list->push_back(dummy);
    begin=atom_list->begin();
    end=atom_list->end();
  }
  N+=NFixed;
  Partition(); ReNeighbor();
  a=t;
  Run(0.05,0,0,0,0,-1,1,0,0);
  Thermo(temp, 0, 0, dt, 0, 0, 0, 0);
  t=a;
}
//*************************************************************************
void config::CheckDepth(){
  for (VAI s=begin; s<end; s++)
    if (s->is_fixed)
      for (VNI sp=s->nlist.begin(); sp!=s->nlist.end(); sp++){
	if (((atom*)sp->a2)->id != s->id){
	  s=end;
	  cerr<<"* Detected foreign species in proximity of fixed layer.\n";
	  AddFixed(0);
	  break;
	}
      }
}

  
