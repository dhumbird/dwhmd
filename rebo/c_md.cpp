#include "config.h"
extern map<short,double> RC_SQ;

void config::dtOptimize(){
  double d=0;
  for (i=begin; i<end; i++) if (!i->is_fixed) d = mfmax(d,i->V.sqmag());
  if (d==0) d=25;
  Setdt(0.005/sqrt(d));
}
//***************************************************************************
void config::Run(double run_time, int cfgevery, double del_t, bool printe, 
		 bool dump_end, int track, bool quiet, bool tso, bool nosput){
  
  if (run_time==0) run_time=1;
  if (del_t!=0) Setdt(del_t);
  double start=t;
  double finish=start+run_time;
  bool dump=0;
  double next=0;
  if (cfgevery>0)
    next=(double)(((int)(t*1000)+cfgevery - 
		   ((int)(t*1000)+cfgevery)%cfgevery))/1000;
  while (t<finish){
    dump=0;
    if (tso) dtOptimize();
    if (cfgevery>0){
      if (t+dt > next){
	if (tso) Setdt(next-t);
	dump=1;
	next+=(double)cfgevery/1000;
      }
    }
    if (u!=0.0){
      if (!nosput) CheckSput(); 
      FirstVV(); 
    }
    TimeStepInit(); 
    ReNeighbor(); 
    ForceEval(); 
    SecondVV(); 
    t+=dt;
    if (printe) PrintE();
    if (track != -1){
      cout<<t;
      VAI sweet=atomix(track);
      cout<<" "<<sweet->R<<" "<<sweet->V<<endl;
    }
    if (dump) Dump();
  }
  t-=dt;
  if (dump_end) Dump();
}
//***********************************************************************
void config::Thermo(float To, float beren_dt, int cfgevery, double del_t,
		    bool printe, bool dump_end, bool tso, double runtime){
  
  if (del_t!=0) Setdt(del_t);
  double bhb_dt=dt*10;
  double lambda=0;
  bool cool=(To < T);
  double start=t;
  bool dump=0;
  double next=0;
  if (cfgevery>0)
    next=(double)(((int)(t*1000)+cfgevery - 
		   ((int)(t*1000)+cfgevery)%cfgevery))/1000;
  if (runtime!=0){
    double finish=start+runtime;
    while(t<finish){
      dump=0;
      if (tso) dtOptimize();
      if (cfgevery>0){
	if (t+dt > next){
	  if (tso) Setdt(next-t);
	  dump=1;
	  next+=(double)cfgevery/1000;
	}
      }
      if (u!=0.0){
	CheckSput();
	FirstVV();
      }
      TimeStepInit();
      ReNeighbor();
      ForceEval();
      SecondVV();
      for (int z=0; z<Nz; z++){
	float Tz=0; int w=0;
	for (int x=0; x<Nx; x++)
	  for (int y=0; y<Ny; y++){
	    cells_it=&cells[x+Nx*y+Nx*Ny*z];
	    for (sai_j=cells_it->abegin(); sai_j!=cells_it->aend(); sai_j++){
	      Tz+=(*sai_j)->Ek();
	      if (!(*sai_j)->is_fixed) w++;
	    }
	  }
	Tz/=(3*w*KB);
	lambda=sqrt(1+dt/bhb_dt*(To/Tz-1));
	for (int x=0; x<Nx; x++)
	  for (int y=0; y<Ny; y++){
	    cells_it=&cells[x+Nx*y+Nx*Ny*z];
	    for (sai_j=cells_it->abegin(); sai_j!=cells_it->aend(); sai_j++)
	      if (!(*sai_j)->is_fixed) (*sai_j)->V*=lambda;
	  }
      }  
      t+=dt;
      if (dump) Dump();
    }
    t-=dt;
  }
  else{ 
    while ((cool && T > To) || (!cool && T < To)){
      if (tso) dtOptimize();
      if (cfgevery>0){
	if (t+dt > next){
	  if (tso) Setdt(next-t);
	  dump=1;
	  next+=(double)cfgevery/1000;
	}
      }
      if (u!=0.0){
	CheckSput();
	FirstVV();
      }
      TimeStepInit();
      ReNeighbor();
      ForceEval();
      SecondVV();
      for (int z=0; z<Nz; z++){
	float Tz=0; int w=0;
	for (int x=0; x<Nx; x++)
	  for (int y=0; y<Ny; y++){
	    cells_it=&cells[x+Nx*y+Nx*Ny*z];
	    for (sai_j=cells_it->abegin(); sai_j!=cells_it->aend(); sai_j++){
	      Tz+=(*sai_j)->Ek();
	      if (!(*sai_j)->is_fixed) w++;
	    }
	  }
	Tz/=(3*w*KB);
	lambda=sqrt(1+dt/bhb_dt*(To/Tz-1));
	for (int x=0; x<Nx; x++)
	  for (int y=0; y<Ny; y++){
	    cells_it=&cells[x+Nx*y+Nx*Ny*z];
	    for (sai_j=cells_it->abegin(); sai_j!=cells_it->aend(); sai_j++)
	      if (!(*sai_j)->is_fixed) (*sai_j)->V*=lambda;
	  }
      }  
      t+=dt;
      if (dump) Dump();
    }  
    t-=dt;
  }
  Run(0.1, cfgevery, 0, printe, 0, -1, 1, tso,0);
  if (dump_end) Dump();
  cerr<<"* Cell thermalized. ("<<T<<" K). Time [ps]: "<<t<<endl;
}
//**********************************************************************
void config::RunQuenched(double run_time, int cfgevery, double del_t, 
                         bool printe, bool dump_end, int track, bool quiet, 
                         bool tso, bool nosput, float To){
  
  //no-frills integration

  if (run_time==0) run_time=1;
  if (del_t!=0) Setdt(del_t);
  double bhb_dt=dt*10;
  double start=t;
  double finish=start+run_time;
  double lambda=0;
  if (!quiet) cerr<<"* Quenched integration started. Cfg time [ps]: "<<t<<endl;
  bool dump=0;
  double next=0;
  if (cfgevery>0)
    next=(double)(((int)(t*1000)+cfgevery - 
		   ((int)(t*1000)+cfgevery)%cfgevery))/1000;
  while (t<finish+dt){
    if (tso) dtOptimize();
    if (cfgevery>0){
      if (t+dt > next){
	if (tso) Setdt(next-t);
	dump=1;
	next+=(double)cfgevery/1000;
      }
    }
    if (u!=0.0){
      if (!nosput) CheckSput();
      FirstVV();
    }
    TimeStepInit();
    ReNeighbor();
    ForceEval();
    SecondVV();
    lambda=sqrt(1+dt/bhb_dt*(To/T-1));
    for (i=begin; i<end; i++)
      if (!i->is_fixed)
        i->V*=lambda;
    t+=dt;
    if (printe) PrintE();
    if (track!= -1 )cout<<t<<"\t"<<atomix(track)->NC<<endl;
    if (dump) Dump();
  }
  t-=dt;
  if (!quiet) cerr<<"* Quenched integration ended. Cfg time [ps]: "<<t<<endl;
  if (dump_end) Dump();
}


//**********************************************************************
void config::ReNeighbor(){
  for (i=begin; i<end; i++) i->nlist.clear();
  if (inert!= -1){
    i=atomix(inert);
    for (vai_j=begin; vai_j<end; vai_j++){
      if (i!=vai_j){
	R = i->R-vai_j->R; R.minimg(Lx,Ly);
	if ((r=R.sqmag()) < 36){ //6 A cutoff for inerts
	  nbr n;
	  n.a1=(void*)&(*i);
	  n.a2=(void*)&(*vai_j);
	  r=sqrt(r);
	  n.r=r;
	  n.Rhat=R/r;
	  i->nlist.push_back(n);
	}
      }
    }
  }
  for (vscit=cells.begin(); vscit!=cells.end(); vscit++)
    for (sai_i=vscit->abegin(); sai_i!=vscit->aend(); sai_i++){
      ii=*sai_i;
      for (c=vscit->nbegin(); c!=vscit->nend(); c++)
	for (sai_j=(*c)->abegin(); sai_j!=(*c)->aend(); sai_j++){
	  j=*sai_j;
	  if (ii < j && ii->ix!=inert && j->ix!=inert){
	    switch (ii->id + j->id){
	    case 12: //C-C
	      type=0; break;
	    case 2:  //H-H
	      type=1; break;
	    case 7:  //C-H
	      type=2; break;
	    }
	    R = ii->R-j->R; R.minimg(Lx,Ly);
	    if ((r=R.sqmag()) < RC_SQ[type]){
	      nbr n;
	      n.type=type;
	      r=sqrt(r);
	      n.a1=(void*)ii;
	      n.a2=(void*)j;
	      n.r=r;
	      n.Rhat=R/r;
	      n.PreComp();
	      ii->nlist.push_back(n);
	      n.Invert();
	      j->nlist.push_back(n);
	    }
	  }
	}
    }
  for (i=begin; i<end; i++) {
    if (i->ix != inert) 
      i->PostComp();
  }
}
//**********************************************************************

