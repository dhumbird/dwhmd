#include "config.h"
#include "errors.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <list>

void MainFStudy(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_000000.cfg";
  float runtime=0.2; float e=200; float Tion=300;
  float inc=0; float azi=0; float T=300;
  double dt=0; double bhb_dt=0;
  int runs=1; int id=18; int begin=1;
  bool tso=0;
  bool choose_inc=0; bool choose_azi=0;
  short cfgthermo=1; bool pickup=0;
  bool refresh=0; bool dump_each=0;
    
  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-runs") runs=atoi(argv[++i]);
    else if (arg=="-e") e=atof(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
    else if (arg=="-Tion") Tion=atof(argv[++i]);
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-pickup") pickup=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-inc"){
      arg=argv[++i];
      if (arg=="rand") choose_inc=1;
      else inc=atof(arg.c_str());
    }
    else if (arg=="-azi"){
      arg=argv[++i];
      if (arg=="rand") choose_azi=1;
      else azi=atof(arg.c_str());
    }
    else CmdError(arg.c_str());
  }

  //**************************************************************************
  
  for (int run=begin; run<=begin+runs-1; run++){
    config cfg(scfg);
    srand(TimeSeed());
    float top=cfg.AvgTop();
    VPI vpi_F, vpi_Si;
    int int_F, int_Si;
    bool goodatom=0;
    while(!goodatom){
      int_F=(int)floor(rand01()*cfg.N);
      vpi_F=cfg.atomix(int_F);
      if (vpi_F!=cfg.end){
	if (!vpi_F->is_fixed && vpi_F->R.z < top && vpi_F->id==cfg.mat){
	  vpi_F->id=9;
	  vpi_F->SetProps();
	  goodatom=1;
	}
      }
    }
    goodatom=0;
    cfg.Run(0.1, 0, 0, 0, 0, -1, 1, 0,0);
    vpi_F=cfg.atomix(int_F);
    if (vpi_F!=cfg.end){
      double d1=1000; double d2=1000;
      for (SPI s=vpi_F->nlist.begin(); s!=vpi_F->nlist.end(); s++){
	d2=d1;
	d1 = d1 <? ((*s)->R - vpi_F->R).minsqmag(cfg.Lx,cfg.Ly);
	if (d2>d1) int_Si=(*s)->ix;
      }
      vpi_Si=cfg.atomix(int_Si);
      if (vpi_Si!=cfg.end)
	if (!vpi_Si->is_fixed)
	  goodatom=1;
    }
    if (!goodatom) run--;
    else{
      cerr<<"------------"<<DateTime()<<"--------run "<<run
          <<"-----------------------"<<endl;
      cerr<<"* F: "<<int_F<<"  Partner: "<<int_Si<<endl;
      int maxbomb=50; int bomb=0;
      vpi_Si->O=vpi_Si->R;
      vpi_Si->xpass=0; vpi_Si->ypass=0;
      vpi_F->O=vpi_F->R;
      vpi_F->xpass=0; vpi_F->ypass=0;
      string outfile=cfg.name+"_"+itoa(run,4)+".fst";
      fstream fout(outfile.c_str(), ios::out);
      cfg.TimeStepInit();
      while(cfg.atomix(int_F)!=cfg.end 
	    && cfg.atomix(int_Si)!=cfg.end 
	    && bomb<maxbomb){
	bomb++;
	cfg.AddIon(e, id, Tion, choose_inc, inc, choose_azi, azi, -1, 1,0,0,0);
	double finish=cfg.t+runtime;
	for (VPI I=cfg.begin; I<cfg.end; I++) I->P=I->R;
	while (cfg.t<finish){
	  if (tso) cfg.dtOptimize();
	  if (cfg.u!=0){
	    cfg.CheckSput();
	    cfg.FirstVV();
	  }
	  cfg.TimeStepInit();    
	  cfg.ReNeighbor();    
	  cfg.ForceEval();
	  cfg.SecondVV();
	  cfg.t+=cfg.dt;
	}
	cfg.t-=cfg.dt;
	cfg.DelIon(1);
	cfg.Thermo(T, bhb_dt, 0, dt, 0, 0, tso);
	vpi_Si=cfg.atomix(int_Si);
	vpi_F=cfg.atomix(int_F);
	if (cfg.atomix(int_F)!=cfg.end  && cfg.atomix(int_Si)!=cfg.end){
	  int partner=0;
	  double d1=1000; double d2=1000;
	  for (SPI s=vpi_F->nlist.begin(); s!=vpi_F->nlist.end(); s++){
	    d2=d1;
	    d1 = d1 <? ((*s)->R - vpi_F->R).minsqmag(cfg.Lx,cfg.Ly);
	    if (d2>d1) partner=(*s)->ix;
	  }
	  if (partner==int_Si){
	    svector P(vpi_Si->xpass*cfg.Lx, vpi_Si->ypass*cfg.Ly, 0);
	    fout<<bomb<<"\t"<<(vpi_Si->R+P-vpi_Si->O).sqmag()<<endl;
	  }
	  else{
	    cerr<<"* F switched partners. Doing over.\n";
	    run--;
	    break; break; break;
	  }
	}
      }
    }
  }
  cerr<<"---------------------------------------------------------------------------"<<endl;
  cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
} 



