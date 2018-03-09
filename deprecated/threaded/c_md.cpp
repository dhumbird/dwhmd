#include "config.h"

//*******************************************************************
void config::TimeStepInit(){
  u=0; ek=0;
  for (i=begin; i<end; i++) i->ts_init();
  for (a=0; a<Nc; a++) cells[a].u=0;
}
//**********************************************************************
void config::ReNeighbor(){
  Allocate();
  pthread_create(&threads[0], NULL, subcell::th_s_ReNeighbor, &cellset1);
  pthread_join(threads[0], (void**) &status);
  pthread_create(&threads[1], NULL, subcell::th_s_ReNeighbor, &cellset2);
  pthread_join(threads[1], (void**) &status);
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
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump(0);
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
  if (dump_end) Dump(0);
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
    if (cfgevery>0) if ((int)(t*1000)%cfgevery==0) Dump(0);
    TimeStepInit();
    ReNeighbor();
    ForceEval(0);
    SecondVV();
    if (printe) PrintE();
    t+=dt;
  }
  t-=dt;
  cerr<<"* Normal integration ended. Cfg time [ps]: "<<t<<endl;
  if (dump_end) Dump(0);
}
//***********************************************************************
bool config::CheckSput(VPI s){
 
  /*delete an atom if it:
    has zero force on it
    -AND has punched through the bottom in an open cell
    -AND is the same chemical as the crystal AND has z-coordinate > Lz/2 
    (original top surface)*/

  if (!s->is_fixed){
    if (s->F.is_zero()){
      if(s->id==mat && s->R.z < -1*Lz2){
	cerr<<"* Atom passed through bottom at "<<s->R<<" time: "<<t<<endl;
	erase(s);
	RePartition();
	return 0;
      }
      else if (s->id==mat && s->R.z > Lz2){
	cerr<<"* Atom sputtered from "<<s->R<<" time: "<<t<<endl;
	erase(s);
	RePartition();
	Nsput++;
	return 0;
      }
      else return 1;
    }
    else return 1;   
  }
  else return 1;
}




