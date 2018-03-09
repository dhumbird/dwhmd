#include "config.h"
#include "errors.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <list>

void MainIonBombard(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_000000.cfg";
  float runtime=0.2; float e_ion=200; float Tion=300;
  float inc=0; float azi=0; float T=300; 
  double dt=0; double bhb_dt=0;
  int cfgevery=0; int runs=1; int ion_id=18; int begin=1;
  bool printe=0; bool dump_end=0; bool tso=0;
  bool choose_inc=0; bool choose_azi=0;
  short cfgthermo=1; bool pickup=0;
  bool refresh=0; bool dump_each=0;
  list<int> seeds;
  int pass=0;
  int neutral_id=0; int n2i=0;
  float e_neutral=0;

  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-runs") runs=atoi(argv[++i]);
    else if (arg=="-e"||arg=="-eion") e_ion=atof(argv[++i]);
    else if (arg=="-eneu") e_neutral=atof(argv[++i]);
    else if (arg=="-n2i") n2i=atoi(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
    else if (arg=="-Tion") Tion=atof(argv[++i]);
    else if (arg=="-id"||arg=="-ion") ion_id=atoi(argv[++i]);
    else if (arg=="-neut") neutral_id=atoi(argv[++i]);
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-pickup") pickup=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-pass") pass=atoi(argv[++i]);
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-refresh") refresh=1;
    else if (arg=="-dump"){
      arg=argv[++i];
      if (arg=="end") dump_end=1;
      else if (arg=="each") dump_each=1;
      else cfgevery=atoi(arg.c_str());
    }
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
    else if (arg=="-seeds")
      for (int j=0; j<runs; j++)
	seeds.push_back(atoi(argv[++i]));
    else CmdError(arg.c_str());
  }

  if (seeds.empty()) 
    for (int j=0; j<runs; j++)
      seeds.push_back(0);

  //**************************************************************************
  
  if (refresh){
    for (int run=begin; run<=begin+runs-1; run++){
      cerr<<"------------"<<DateTime()<<"--------run "<<run
	  <<"-----------------------"<<endl;
      {
	config cfg(scfg);
	cfg.AddIon(e_ion, ion_id, Tion, choose_inc, inc, choose_azi, azi,0,0 );
	VPI I;
	cfg.Run(runtime, 0, dt, printe, 0, -1, 0, tso,0);
	if (dump_each){
	  cfg.DelIon(0);
	  cfg.Dump();
	}
      }
    }
  }
  else{
    config cfg(scfg);
    list<int>::iterator seed=seeds.begin();
    for (int run=begin; run<=begin+runs-1; run++){
      cerr<<"------------"<<DateTime()<<"--------run "<<run
	  <<"-----------------------"<<endl;
      if (pass!=0){
	cfg.Passivate(pass);
	cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso);
      }
      for (int n=0; n<n2i; n++){
	cfg.AddIon(e_neutral, neutral_id, Tion, 0, 0, 0, 0, 0, 0);
	double start=cfg.t;
	double finish=start+runtime;
	for (VPI I=cfg.begin; I<cfg.end; I++) I->P=I->R;
	while (cfg.t<finish){
	  if (tso) cfg.dtOptimize();
	  if (cfg.u!=0){
	    cfg.CheckSput();
	    cfg.FirstVV();
	  }
	  if (cfgevery>0) if ((int)(cfg.t*1000)%cfgevery==0) cfg.Dump();
	  cfg.TimeStepInit();    
	  cfg.ReNeighbor();    
	  cfg.ForceEval();
	  cfg.SecondVV();
	  cfg.t+=cfg.dt;
	}
	cfg.t-=cfg.dt;
	cfg.Thermo(T, bhb_dt, cfgthermo*cfgevery, dt, printe, 0, tso);
      }
      cfg.AddIon(e_ion, ion_id, Tion, choose_inc, inc, choose_azi, azi, 
		 *seed, 0);
      double start=cfg.t;
      double finish=start+runtime;
      for (VPI I=cfg.begin; I<cfg.end; I++) I->P=I->R;
      while (cfg.t<finish){
	if (tso) cfg.dtOptimize();
	if (cfg.u!=0){
	  cfg.CheckSput();
	  cfg.FirstVV();
	}
	if (cfgevery>0) if ((int)(cfg.t*1000)%cfgevery==0) cfg.Dump();
	cfg.TimeStepInit();    
	cfg.ReNeighbor();    
	cfg.ForceEval();
	cfg.SecondVV();
	cfg.t+=cfg.dt;
      }
      cfg.t-=cfg.dt;
      cerr<<"* Normal integration ended. Time: "<<cfg.t<<endl;
      cfg.DelIon(0);
      //while(!cfg.Sweep());
      cfg.Thermo(T, bhb_dt, cfgthermo*cfgevery, dt, printe, 0, tso);
      cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
      if (dump_each) cfg.Dump();
      seed++;
    }
    if (dump_end) cfg.Dump();
    cerr<<"---------------------------------------------------------------------------"<<endl;
    cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
  } 
}

