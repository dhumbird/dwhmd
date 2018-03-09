#include "config.h"

//**********************************************************************
void config::ReNeighbor(){
  for (cells_it=cells.begin(); cells_it!=cells.end(); cells_it++)
    for (i0=cells_it->abegin(); i0!=cells_it->aend(); i0++)  
      for (c=cells_it->nbegin(); c!=cells_it->nend(); c++)
	for (j0=(*c)->abegin(); j0!=(*c)->aend(); j0++)
	  if (*i0 < *j0)
	    if (((*i0)->R - (*j0)->R).minsqmag(Lx, Ly, Lz) < 
		((*i0)->rc >? (*j0)->rc)){
	      (*i0)->nlist.insert(*j0);
	      (*j0)->nlist.insert(*i0);
	    }
}
//**************************************************************************
void config::Thermo(float To, float beren_dt, int cfgevery, double del_t,
		    bool printe, bool tso, bool dump_end){
  
  if (beren_dt!=0) bhb_dt=beren_dt;
  if (del_t!=0) Setdt(del_t);

  //integration with thermostatting-stops when T crosses To
 
  bool cool=(To < T);

  while ((cool && T > To) || (!cool && T < To)){
    if (tso) dtOptimize();
    if (u!=0.0) FirstVV();
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump();
    TimeStepInit();
    ReNeighbor();
    ForceEval(0);
    SecondVV();
    d1=sqrt(1.0+dt/bhb_dt*(To/T-1.0));
    for (i=begin; i<end; i++) if (!i->is_fixed) i->V*=d1;
    if (printe) PrintE();
    t+=dt;
  }  
  t-=dt;
  if (dump_end) Dump();
  cerr<<"* Cell thermalized. ("<<To<<" K). Time [ps]: "<<t<<endl;
}
//***************************************************************************
void config::dtOptimize(){
  d1=0;
  for (i=begin; i<end; i++) if (!i->is_fixed) d1 = d1 >? i->V.sqmag();
  if (mat==14) Setdt(0.01*SIGMA_SW/sqrt(d1));
}
//***************************************************************************
void config::Run(double run_time, int cfgevery, double del_t, bool printe, 
		 bool tso, bool dump_end){
  
  //no-frills integration

  if (run_time==0) run_time=1;
  if (del_t!=0) Setdt(del_t);
  
  double start=t;
  double finish=start+run_time;

  cerr<<"* Normal integration started. Cfg time [ps]: "<<t<<endl;
  
  while (t<finish+dt){
    if (tso) dtOptimize();
    if (u!=0.0) FirstVV();
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump();
    TimeStepInit();
    ReNeighbor();
    ForceEval(0);
    SecondVV();
    if (printe) PrintE();
    t+=dt;
  }
  t-=dt;
  cerr<<"* Normal integration ended. Cfg time [ps]: "<<t<<endl;
  if (dump_end) Dump();
}
//***********************************************************************




