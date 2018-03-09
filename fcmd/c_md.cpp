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
void config::ReNeighbor(){
  int type; svector R; double r;
  for (VAI i=begin; i<end; i++){
    i->nlist.clear();
    i->Nmap[6]=0; i->Nmap[9]=0; i->Nmap[14]=0; i->Nmap[17]=0;
  }
  if (inert!= -1){
    VAI ion=atomix(inert);
    ion_nbr.clear();
    for (VAI j=begin; j<end; j++){
      if (ion!=j){
	R = ion->R - j->R; R.minimg(Lx,Ly);
	if ((r=R.sqmag()) < 36){ //6 A cutoff for inerts
	  ion_nbr.insert(&(*j));
	}
      }
    }
  }
  //cerr<<"inert ";
  for (vscit=cells.begin(); vscit!=cells.end(); vscit++)
    for (sai_i=vscit->abegin(); sai_i!=vscit->aend(); sai_i++){
      atom* i=*sai_i;
      if (i->ix != inert)
	for (c=vscit->nbegin(); c!=vscit->nend(); c++)
	  for (sai_j=(*c)->abegin(); sai_j!=(*c)->aend(); sai_j++){
	    atom* j=*sai_j;
	    if (i < j && j->ix!=inert){
	      type = i->id+j->id;
	      R = i->R-j->R; R.minimg(Lx,Ly);
	      if (i->ix==0) cerr<<R.mag()<<" ";
	      if ((r=R.sqmag()) < RC_SQ[type]){
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
      if (i->ix ==0) cerr<<endl;
    }
  for (VAI i=begin; i<end; i++) if (i->ix != inert){
    i->AddNbrs();
  }
}
//**********************************************************************
//***************************************************************************
void config::bonddist(double run_time,double del_t,bool tso,short mode){
  
  double rho=0.049853;
  double a=pow(8.0/rho, 1.0/3.0);
  svector k(-1,1,-1);
  k*=(2*PI/a);
  svector kz(0,0,1);
  float h=.5*a;
  double rij=0;
  int bin=0;
  int num_bins=(int)(Lz/h)+1;
  vector<float> opz(num_bins);
  vector<float> lbonds(num_bins);
  vector<float> meanang(num_bins);
  vector<float> xi(num_bins);
  vector<float> lambda(num_bins);
	
  for (int z=0; z<num_bins; z++){
    opz[z]=0;
    lbonds[z]=0;
    meanang[z]=0;
    xi[z]=0;
    lambda[z]=0;
  }

  if (run_time==0) run_time=1;
  if (del_t!=0) Setdt(del_t);
  double start=t;
  double finish=start+run_time;
  bool dump=0;
  double next=0;
  int steps=0;
  while (t<finish){
    dump=0;
    if (tso) dtOptimize();
    if (u!=0.0) FirstVV();
    
    TimeStepInit();
    ReNeighbor();
    ForceEval();
    SecondVV(); 
    t+=dt;
    int here=0;
    for (int z=0; z<num_bins; z++){
      float opcos=0; float opsin=0; int opatoms=0;
      int numbonds=0; float sumbonds=0;
      int numangs=0; float sumangs=0;
      float newxi=0;
      float lam=0; double lmin; double lmax;
      for (VAI i_=begin; i_<end; i_++){
	if (!i_->is_fixed && i_->id==14){
	  if (i_->R.z > the_bottom+h*(z) && i_->R.z < the_bottom+h*(z+1)){
	    //if (i_->R.z > the_bottom+h*(z-1) && i_->R.z < the_bottom+h*(z+2)){
	    //opcos+=cos(k^i_->R);
	    //opsin+=sin(k^i_->R);
	    lam+=cos(8*PI*i_->R.x/a)+cos(8*PI*i_->R.y/a)+cos(8*PI*i_->R.z/a);
	    //cerr<<l<<endl;
	    opatoms++;
	    for (VNI ij=i_->nlist.begin(); ij!=i_->nlist.end(); ij++){
	      if (ij->type==28){
		if (ij->bo >0.9){
		  //cerr<<ij->Rhat*PI/a<<endl;
		  //float xx=cos(ij->Rhat.x/
		  //cerr<<ij->Rhat<<" "<<" "<<(ij->Rhat^k)*a*sqrt(3)<<endl;
		  //newxi+=pow((ij->Rhat^k),2);
		  newxi+=cos((ij->Rhat^k)*a*sqrt(3));
		  numbonds++;
		  sumbonds+=ij->r;
		}
		atom* j = (atom*)ij->a2;
		for (VNI bond_ik=i_->nlist.begin(); bond_ik!=i_->nlist.end(); bond_ik++){
		  if (bond_ik->type==28){
		    if ((atom*)bond_ik->a2 !=j && bond_ik->bo>0.9){
		      sumangs+=acos(ij->Rhat ^ bond_ik->Rhat);
		      numangs++;
		    }	
		  }
		}
		for (VNI bond_jk=j->nlist.begin(); bond_jk!=j->nlist.end(); bond_jk++){
		  if (bond_jk->type==28){
		    if ((atom*)bond_jk->a2 !=&(*i_) && bond_jk->bo > 0.9){
		      sumangs+= acos(-1*ij->Rhat ^ bond_jk->Rhat);
		      numangs++;
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      //if(z==num_bins-1) cerr<<opatoms<<endl;
      if (opatoms!=0){
	opz[z]+=sqrt(opcos*opcos+opsin*opsin)/opatoms;
	lambda[z]+=lam/opatoms/3;
      }
      else opz[z]=0;
      if (numbonds!=0){
	lbonds[z]+=sumbonds/numbonds;
	//cerr<<newxi<<" "<<numbonds<<endl;
	//xi[z]=0.5*(3*newxi/numbonds-1);
	xi[z]+=newxi/numbonds;
      }
      else lbonds[z]=0;
      if (numangs!=0){
	meanang[z]+=sumangs/numangs;
      }
      else meanang[z]=0;
    }
    steps++;
  }
  t-=dt;
  float height=num_bins*h;
  //float height=AvgTop()-the_bottom;
  for (int z=num_bins-1; z>0; z--){
    opz[z]/=steps;
    float here=height-(z+1)*h;
    printf("%8.5f  %8.5f  %8.5f  %8.5f  %8.5f\n",here, lambda[z]/steps, xi[z]/steps, lbonds[z]/steps, meanang[z]/steps/PI*180);
  }
  if (mode==1){
    float adepth=0;
    for (int z=num_bins-1; z>0; z--){
      if (lambda[z]<0.15){
	float x1=height-(z+1)*h;
	float x2=height-z*h;
	adepth=x1+(0.15-lambda[z])/(lambda[z-1]-lambda[z])*(x2-x1);
      }
    }
    cout<<"crystal interface depth: "<<adepth<<endl;
  }
}
