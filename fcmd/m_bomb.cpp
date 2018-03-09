#include "config.h"
#include "errors.h"
#include "consts.h"
#include "sysfun.h"
#include <list>

void MainIonBombard(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_0000-000.cfg";
  float runtime=0.5; float e_ion=0; float Tion=300;
  float inc=0; float azi=0; float T=300; 
  double dt=0; double bhb_dt=0;
  int cfgevery=0; int runs=1; int ion_id=18; int begin=1;
  bool printe=0; bool dump_end=0; bool tso=0;
  bool choose_inc=0; bool choose_azi=0; bool dump_only=0;
  short cfgthermo=1;
  bool refresh=0; bool dump_each=0;
  list<int> seeds;
  string species="unknown";
  bool fancy=1;
  bool b_test=0; //re-desorption not a problem except for H
  float confine=0;

  for (i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-runs"||arg=="-rounds") runs=atoi(argv[++i]);
    else if (arg=="-e"||arg=="-eion") e_ion=atof(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
    else if (arg=="-Tion") Tion=atof(argv[++i]);
    else if (arg=="-type") species=argv[++i];
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-tso") tso=1;
    else if (arg=="-notest") b_test=0;
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-refresh") refresh=1;
    else if (arg=="-nonfancy") fancy=0;
    else if (arg=="-confine") confine=atof(argv[++i]);
    else if (arg=="-dump"){
      arg=argv[++i];
      if (arg=="end") dump_end=1;
      else if (arg=="each") dump_each=1;
      else if (arg=="only") dump_only=1;
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
 
  species=RectifySpecies(species);
  if (e_ion!=0) fancy=0;
  if (species=="F2") b_test=1;
  //*********************************************************
  if (fancy){
    if (refresh){
      int kept=0;
      float nbegin, nend;
      list<int>::iterator seed=seeds.begin();
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	{
	  config cfg(scfg);
	  int No=cfg.N;
	  VAI I=cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
			   azi,*seed,0,confine);
	  int Hix=I->ix;
	  nbegin=I->Nmap[6]+I->Nmap[14];
	  bool quit=0;
	  while (!quit){
	    //step 10 fs at a time.
	    cfg.Run(0.01, 0, dt, printe, 0, -1, 0, tso,1);
	    I=cfg.atomix(Hix);
	    if (I==cfg.end){
	      quit=1;
	      cerr<<"* Species left.\n";
	    }
	    else if (I->Nmap[6]+I->Nmap[14] > 0.3){
	      int bondix=-1;
	      float bondr=1000;
	      for (VNI sj=I->nlist.begin(); sj!=I->nlist.end(); sj++){
		atom* j=(atom*)sj->a2;
		if ((I->R - j->R).minsqmag(cfg.Lx,cfg.Ly) < bondr){
		  bondr=(I->R - j->R).minsqmag(cfg.Lx,cfg.Ly);
		  bondix=j->ix;
		}
	      }
	      cerr<<"* Species appears to be bound to # "<<bondix<<" "
		  <<I->R<<endl;
	      cerr<<"  Begin: "<<nbegin<<" End: "
		  <<I->Nmap[6]+I->Nmap[14]<<endl;
	      //make sure it's not going to come off
	      if (b_test) cfg.Run(2, 0, dt, printe, 0, -1, 0, tso,1);
	      if (cfg.atomix(Hix)!=cfg.end) kept++;
	      quit=1;
	    }
	  }
	}
	seed++;
      }
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
      cerr<<kept<<" "<<species<<"'s were deposited,";
      cerr<<" Sticking coefficient: "<<(float)kept/(float)runs<<endl;
    } //fancy-refresh
    else{ //fancy-cumulative
      config cfg(scfg);
      list<int>::iterator seed=seeds.begin();
      int kept=0;
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	int No=cfg.N;
	VAI I=cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
			 azi, *seed, 0, confine);
	for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;
	int Hix=I->ix;
	double nbegin=I->Nmap[6]+I->Nmap[14];
	bool quit=0;
	while (!quit){
	  //step 50 fs at a time.
	  cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,1);
	  I=cfg.atomix(Hix);
	  if (I==cfg.end){
	    quit=1;
	    cerr<<"* Species left.\n";
	    if (cfg.N!=No && dump_only) cfg.Dump();
	    if (fabs(1-T/cfg.T) > 0.05)
	      cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, 0, 0);
	  }
	  else if (I->Nmap[6]+I->Nmap[14] > 0.3){
	    int bondix=-1;
	    float bondr=1000;
	    for (VNI sj=I->nlist.begin(); sj!=I->nlist.end(); sj++){
	      atom* j=(atom*)sj->a2;
	      if ((I->R - j->R).minsqmag(cfg.Lx,cfg.Ly) < bondr){
		bondr=(I->R - j->R).minsqmag(cfg.Lx,cfg.Ly);
		bondix=j->ix;
	      }
	    }
	    cerr<<"* Species appears to be bound to # "<<bondix<<" "
	        <<I->R<<endl;
	    cerr<<"  Begin: "<<nbegin<<" End: "<<I->Nmap[6]+I->Nmap[14]<<endl;
	    //make sure it's not going to come off
	    if (b_test) cfg.Run(0.5, 0, dt, printe, 0, -1, 0, tso,0);
	    quit=1;
	    cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
	    if (cfg.atomix(Hix)!=cfg.end){
	      kept++;
	      //sweep
	      cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,2);
	    }
	    if (dump_only) cfg.Dump();
	  }
// 	  else if (I->id==17 && I->Nmap[17] > 0.5){
// 	    cerr<<"* Cl got stuck on Cl\n";
// 	    quit=1;
// 	    cfg.erase(I);
// 	  }
	}
	cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
	cfg.CheckDepth();
	if (dump_each) cfg.Dump();
	seed++;
      }
      if (dump_only || dump_end) cfg.Dump();
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
    }
  } //fancy 
  //********************************************************
  else{ //not fancy
    if (refresh){
      list<int>::iterator seed=seeds.begin();
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	{
	  config cfg(scfg);
	  cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
		     azi,*seed,0,confine);
	  for (VAI I=cfg.begin; I<cfg.end; I++) I->P=I->R;
	  cfg.Run(runtime, 0, dt, printe, 0, -1, 0, tso,1);
	  cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,2);
	  if (dump_each){
	    cfg.DelIon(0);
	    cfg.Dump();
	  }
	}
	seed++;
      }
    }
    else{
      config cfg(scfg);
      list<int>::iterator seed=seeds.begin();
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, azi, 
		   *seed, 0,confine);
	for (VAI I=cfg.begin; I<cfg.end; I++) I->P=I->R;
	cfg.Run(runtime, 0, dt, printe, 0, -1, 0, tso,1);
	cerr<<"* Normal integration ended. Time: "<<cfg.t<<endl;
	cfg.DelIon(0);
	cfg.Thermo(T, bhb_dt, cfgthermo*cfgevery, dt, printe, 0, tso, 0.5);
	//sweep
	cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,2);
	cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
	cfg.CheckDepth();
	if (dump_each || dump_only) cfg.Dump();
	seed++;
      }
      if (dump_end) cfg.Dump();
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
    } 
  }
}

