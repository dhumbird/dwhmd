#include "config.h"

//**********************************************************************
void config::ReNeighbor(){
  for (i=begin; i<end; i++) i->nlist.clear();
  for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
    for (i0=cells_it->abegin(); i0!=cells_it->aend(); i0++)  
      for (c=cells_it->nbegin(); c!=cells_it->nend(); c++)
	for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	  if (*i0 < *j0)
	    if (((*i0)->R - (*j0)->R).minsqmag(Lx, Ly) < 
		((*i0)->rc >? (*j0)->rc)){
	      (*i0)->nlist.insert(*j0);
	      (*j0)->nlist.insert(*i0);
	    }
}
//************************************************************************
void config::Thermo(float To, float beren_dt, int cfgevery, double del_t,
		    bool printe, bool dump_end, bool tso){
  
  //Integration with thermostatting.
  //Berendsen factor lambda is calculated on a layer-wise (z-subcell) basis.
  //This leaves an un-Maxwellian config, so the code tracks the time for 
  //T to cross To (setpoint), then runs for another 20% of that interval
  //to fill out the distribution.

  if (beren_dt!=0) bhb_dt=beren_dt;
  if (del_t!=0) Setdt(del_t);
  double lambda=0;
  bool cool=(To < T);
  double start=t;
  subcell* cx_it;
  while ((cool && T > To) || (!cool && T < To)){
    if (tso) dtOptimize();
    if (u!=0.0){
      CheckSput();
      FirstVV();
    }
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump();
    TimeStepInit();
    ReNeighbor();
    ForceEval();
    SecondVV();
    for (int z=0; z<Nz; z++){
      float Tz=0; int w=0;
      for (int x=0; x<Nx; x++)
	for (int y=0; y<Ny; y++){
	  cx_it=&cells[x+Nx*y+Nx*Ny*z];
	  for (j0=cx_it->abegin(); j0!=cx_it->aend(); j0++){
	    Tz+=(*j0)->Ek();
	    if (!(*j0)->is_fixed) w++;
	  }
	}
      Tz/=(3*w*KB);
      lambda=sqrt(1+dt/bhb_dt*(To/Tz-1));
      for (int x=0; x<Nx; x++)
	for (int y=0; y<Ny; y++){
	  cx_it=&cells[x+Nx*y+Nx*Ny*z];
	  for (j0=cx_it->abegin(); j0!=cx_it->aend(); j0++)
	    if (!(*j0)->is_fixed) (*j0)->V*=lambda;
	}
    }
    t+=dt;
  }  
  t-=dt;
  double runtime=(t-start)*0.2;
  Run(runtime, cfgevery, 0, printe, 0, -1, 1, tso,0);
  if (dump_end) Dump();
  cerr<<"* Cell thermalized. ("<<T<<" K). Time [ps]: "<<t<<endl;
}
//*************************************************************************
void config::dtOptimize(){
  double d=0;
  for (i=begin; i<end; i++) if (!i->is_fixed) d = d >? i->V.sqmag();
  if (mat==14) Setdt(0.01*SW_SI_SIGMA/sqrt(d));
}
//***************************************************************************
void config::Run(double run_time, int cfgevery, double del_t, bool printe, 
		 bool dump_end, int track, bool quiet, bool tso, bool nosput){
  
  //no-frills integration

  if (run_time==0) run_time=1;
  if (del_t!=0) Setdt(del_t);
  
  double start=t;
  double finish=start+run_time;

  if (!quiet) cerr<<"* Normal integration started. Cfg time [ps]: "<<t<<endl;
  
  while (t<finish+dt){
    if (tso) dtOptimize();
    if (u!=0.0) {
      if (!nosput) CheckSput();
      FirstVV();
    }
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump();
    TimeStepInit();
    ReNeighbor();
    ForceEval();
    SecondVV();
    if (printe) PrintE();
    if (track!=-1)
      cout<<t<<"\t"<<atomix(track)->R<<"\t"<<atomix(track)->V<<endl;
    t+=dt;
  }
  t-=dt;
  if (!quiet) cerr<<"* Normal integration ended. Cfg time [ps]: "<<t<<endl;
  if (dump_end) Dump();
}
//***********************************************************************
void config::RunQuenched(double run_time, int cfgevery, double del_t, 
			 bool printe, bool dump_end, int track, bool quiet, 
			 bool tso, bool nosput, float To){
  
  //no-frills integration

  if (run_time==0) run_time=1;
  if (del_t!=0) Setdt(del_t);
  
  double start=t;
  double finish=start+run_time;
  double lambda=0;
  if (!quiet) cerr<<"* Quenched integration started. Cfg time [ps]: "<<t<<endl;
  
  while (t<finish+dt){
    if (tso) dtOptimize();
    if (u!=0.0) {
      if (!nosput) CheckSput();
      FirstVV();
    }
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump();
    TimeStepInit();
    ReNeighbor();
    ForceEval();
    SecondVV();
    lambda=sqrt(1+dt/bhb_dt*(To/T-1));
    for (i=begin; i<end; i++)
      if (!i->is_fixed)
	if (i->ix!=inert) i->V*=lambda;
    if (printe) PrintE();
    if (track!=-1)
      cout<<t<<"\t"<<atomix(track)->R<<"\t"<<atomix(track)->V<<endl;
    t+=dt;
  }
  t-=dt;
  if (!quiet) cerr<<"* Quenched integration ended. Cfg time [ps]: "<<t<<endl;
  if (dump_end) Dump();
}




