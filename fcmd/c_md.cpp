#include "config.h"
extern map<short,double> RC_SQ;

void config::dtOptimize(){
  double d=0;
  for (VAI i=begin; i<end; i++) if (!i->is_fixed) d = mfmax(d,i->V.sqmag());
  if (d==0) d=25;
  Setdt(0.05/sqrt(d));
}
//***************************************************************************
void config::Run(double run_time, int cfgevery, double del_t, bool printe, 
		 bool dump_end, int track, bool quiet, bool tso, 
		 short sputmode){
  
  //sputmode: 0=off 1=normal 2=sweep

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
    //cerr<<t<<" ";
    if (u!=0.0){
      if (sputmode){
	//cerr<<"chksp ";
	if (sputmode==1) CheckSput();
	else  Sweep();
      }
      //cerr<<"vv1 ";
      FirstVV();
    }
    //cerr<<"tsi ";
    TimeStepInit();
    //cerr<<"rn ";
    ReNeighbor();
    //cerr<<"fe ";
    ForceEval();
    //cerr<<u<<" vv2 ";
    SecondVV(); 
    //cerr<<T<<" done"<<endl;
    t+=dt;
    if (printe) PrintE();
    if (track != -1){
      cout<<t<<" ";
      VAI sweet=atomix(track);
      cout<<sweet->F;
      cout<<endl;
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
  Run(0.1, cfgevery, 0, printe, 0, -1, 1, tso,1);
  if (dump_end) Dump();
  cerr<<"* Cell thermalized. ("<<T<<" K). Time [ps]: "<<t<<endl;
}
//**********************************************************************
void config::RunQuenched(double run_time, int cfgevery, double del_t, 
                         bool printe, bool dump_end, int track, bool quiet, 
                         bool tso, short sputmode, float To){
  
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
      if (sputmode){
	if (sputmode==1) CheckSput();
	else Sweep();
      }
      FirstVV();
    }
    TimeStepInit();
    ReNeighbor();
    ForceEval();
    SecondVV();
    lambda=sqrt(1+dt/bhb_dt*(To/T-1));
    for (VAI i=begin; i<end; i++)
      if (!i->is_fixed)
        i->V*=lambda;
    t+=dt;
    if (printe) PrintE();
    if (track!= -1 )cout<<t<<"\t"<<atomix(track)->Nmap[6]<<endl;
    if (dump) Dump();
  }
  t-=dt;
  if (!quiet) cerr<<"* Quenched integration ended. Cfg time [ps]: "<<t<<endl;
  if (dump_end) Dump();
}


//**********************************************************************
void config::ReNeighbor()
{
  int type; svector R; double r;
  for (VAI i=begin; i<end; i++)
  {
    i->nlist.clear();
    i->Nmap[6]=0; i->Nmap[9]=0; i->Nmap[14]=0; i->Nmap[17]=0;
  }
  if (inert!= -1)
  {
    VAI ion=atomix(inert);
    ion_nbr.clear();
    for (VAI j=begin; j<end; j++)
    {
      if (ion!=j)
      {
        R = ion->R - j->R; R.minimg(Lx,Ly);
        if ((r=R.sqmag()) < 36)
        { //6 A cutoff for inerts
          ion_nbr.insert(&(*j));
        }
      }
    }
  }
  //cerr<<"inert ";
  for (vscit=cells.begin(); vscit!=cells.end(); vscit++)
  {
    for (sai_i=vscit->abegin(); sai_i!=vscit->aend(); sai_i++)
    {
      atom* i=*sai_i;
      if (i->ix != inert)
      {
      	for (c=vscit->nbegin(); c!=vscit->nend(); c++)
        {
   	      for (sai_j=(*c)->abegin(); sai_j!=(*c)->aend(); sai_j++)
          {
      	    atom* j=*sai_j;
      	    if (i < j && j->ix!=inert)
            {
      	      type = i->id+j->id;
      	      R = i->R-j->R; R.minimg(Lx,Ly);
      	      if ((r=R.sqmag()) < RC_SQ[type])
              {
            		nbr n;
            		n.type=type;
            		r=sqrt(r);
            		n.a1=(void*)i;
            		n.a2=(void*)j;
            		n.r=r;
            		n.Rhat=R/r;
            		n.PreComp();
            		i->nlist.push_back(n);
            		n.Invert();
            		j->nlist.push_back(n);
            		i->Nmap[j->id] += n.f;
            		j->Nmap[i->id] += n.f;
      	      }
            }
      	  }
        }
      }
      //cerr<<i->Nmap[9]<<" "<<i->Nmap[14]<<endl;
    }
  }
  for (VAI i=begin; i<end; i++) if (i->ix != inert)
  {
    i->AddNbrs();
  }
}
//**********************************************************************

