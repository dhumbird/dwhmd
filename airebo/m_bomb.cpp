#include "config.h"
#include "errors.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <list>

struct mixion{
  float ratio;
  float T;
  string type;
};

void MainIonBombard(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_0000-000.cfg";
  float runtime=0; float e_ion=0; float Tion=300;
  float inc=0; float azi=0; float T=300; 
  double dt=0; double bhb_dt=0;
  int cfgevery=0; int runs=1; int ion_id=18; int begin=1;
  bool printe=0; bool dump_end=0; bool tso=0;
  bool choose_inc=0; bool choose_azi=0; bool dump_only=0;
  short cfgthermo=1; bool pickup=0;
  bool refresh=0; bool dump_each=0;
  list<int> seeds;
  string species="unknown";
  string mix="";
  bool mixed=0;
  bool b_test=1;

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
    else if (arg=="-pickup") pickup=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-notest") b_test=0;
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-refresh") refresh=1;
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
    else if (arg=="-mix") {
      mix=argv[++i];
      mixed=1;
    }
    else CmdError(arg.c_str());
  }

  if (seeds.empty()) 
    for (int j=0; j<runs; j++)
      seeds.push_back(0);
 
  species=RectifySpecies(species);
 
  //****************HYDROGEN ONLY******************************
  if (species=="H"){
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
			   azi,*seed,0 );
	  int Hix=I->ix;
	  nbegin=I->NC;
	  bool quit=0;
	  while (!quit){
	    //step 10 fs at a time.
	    cfg.Run(0.01, 0, dt, printe, 0, -1, 0, tso,0);
	    I=cfg.atomix(Hix);
	    if (I==cfg.end){
	      quit=1;
	      if (cfg.N==No) cerr<<"* Species was scattered.\n";
	      else{ 
		cerr<<"* Species abstracted something.\n";
	      }
	    }
	    else if (I->NC > 0.5){
	      cerr<<"* Species appears to be bound. "<<I->R<<endl;
	      cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
	      //make sure it's not going to come off
	      if (b_test) cfg.Run(0.5, 0, dt, printe, 0, -1, 0, tso,0);
	      if (cfg.N!=No) kept++;
	      quit=1;
	    }
	  }
	}
	seed++;
      }
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
      cerr<<kept<<" hydrogens were deposited,";
      cerr<<" Sticking coefficient: "<<(float)kept/(float)runs<<endl;
    } // hydrogen - refresh
    //*******
    else{ //hydrogen - cumulative
      config cfg(scfg);
      list<int>::iterator seed=seeds.begin();
      int kept=0;
     
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	int No=cfg.N;
	VAI I=cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
			 azi, *seed, 0);
	for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;
	int Hix=I->ix;
	double nbegin=I->NC;
	bool quit=0;
	while (!quit){
	  //step 50 fs at a time.
	  cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,0);
	  I=cfg.atomix(Hix);
	  if (I==cfg.end){
	    if (cfg.N==No){
	      cerr<<"* Species was scattered.\n";
	      cfg.Nmax--;
	    }
	    else{
	      cerr<<"Species abstracted something.\n";
	      if (dump_only) cfg.Dump();
	    }
	    quit=1;
	  }
	  else if (I->NC > 0.5){
	    cerr<<"* Species appears to be bound. "<<I->R<<endl;
	    cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
	    //make sure it's not going to come off
	    if (b_test) cfg.Run(0.5, 0, dt, printe, 0, -1, 0, tso,0);
	    quit=1;
	    cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
	    if (cfg.atomix(Hix)!=cfg.end){
	      kept++;
	      if (dump_only) cfg.Dump();
	    }
	  }
	}
	cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
	if (dump_each) cfg.Dump();
	seed++;
      }
      if (dump_only || dump_end) cfg.Dump();
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
    }
  } //species is hydrogen
  
  //**********************HYDROCARBON ONLY************************************
  else if (species=="CH3" || species=="CH2" || species=="CH"){
    if (refresh){
      int kept=0;
      float nbegin;
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	{
	  config cfg(scfg);
	  VAI I=cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
			   azi,0,0 );
	  int Hix=I->ix;
	  nbegin=I->NC;
	  bool quit=0;
	  while (!quit){
	    //step 50 fs at a time.
	    cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,0);
	    I=cfg.atomix(Hix);
	    if (I==cfg.end){
	      cerr<<"* Species was scattered.\n";
	      quit=1;
	    }
	    else if (I->NC > 0.3){
	      cerr<<"* Species appears to be bound. "<<I->R<<endl;
	      cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
	      kept++;
	      quit=1;
	    }
	  }
	}
      }
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
      cerr<<kept<<" methyls were deposited,";
      cerr<<" Sticking coefficient: "<<(float)kept/(float)runs<<endl;
    } // hydrocarbon ions -- refresh 
    //*******
    else{ //hydrocarbon -- cumulative
      config cfg(scfg);
      list<int>::iterator seed=seeds.begin();
      int kept=0;
      for (int run=begin; run<=begin+runs-1; run++){
	cerr<<"------------"<<DateTime()<<"--------run "<<run
	    <<"-----------------------"<<endl;
	int No=cfg.N;
	VAI I=cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
			 azi, *seed, 0);
	for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;
	int Hix=I->ix;
	double nbegin=I->NC;
	bool quit=0;
	while (!quit){
	  //step 50 fs at a time.
	  cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,0);
	  I=cfg.atomix(Hix);
	  if (I==cfg.end){
	    if (cfg.N==No){
	      cerr<<"* Species was scattered.\n";
	    }
	    else{
	      cerr<<"* Species abstracted something.\n";
	      if (dump_only) cfg.Dump();
	    }
	    cfg.Nmax-=4;
	    quit=1;
	  }
	  else if (I->NC > 0.3){
	    cerr<<"* Species appears to be bound. "<<I->R<<endl;
	    cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
	    kept++;
	    quit=1;
	    cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
	    if (dump_only) cfg.Dump();
	  }
	}
	cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
	if (dump_each) cfg.Dump();
	seed++;
      }
      if (dump_end || dump_only) cfg.Dump();
      cerr<<"---------------------------------------------------------------------------"<<endl;
      cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
    }
  } //hydrocarbon ion 
  //**************************INERT ONLY*********************************
  else if (species=="Ar"){
    config cfg(scfg);
    list<int>::iterator seed=seeds.begin();
    if (runtime==0) runtime=0.5;
    for (int run=begin; run<=begin+runs-1; run++){
      cerr<<"------------"<<DateTime()<<"--------run "<<run
	  <<"-----------------------"<<endl;
     
      VAI I=cfg.AddIon(e_ion, species, Tion, choose_inc, inc, choose_azi, 
		       azi, *seed, 0);
      for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;
      cfg.Run(runtime, 0, dt, printe, 0, -1, 0, tso,0);
      cfg.DelIon(0);
      cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0);
      cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
      if (dump_each) cfg.Dump();
      seed++;
    }
    if (dump_end) cfg.Dump();
    cerr<<"---------------------------------------------------------------------------"<<endl;
    cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
  } //inert ion
  //***********************MIXED IONS***********************************
  else if (mixed){
    vector<mixion> ionbag;
    string::size_type p=0;
    string::size_type h=0;
    string m;
    vector<string> level1;
    for (p=0; p<mix.length(); p++){
      if(mix[p]==':'){
	level1.push_back(mix.substr(h,p-h));
	h=p+1;
      }
    }
    level1.push_back(mix.substr(h,1000));
    for (vector<string>::iterator vv=level1.begin(); vv!=level1.end(); vv++){
      m=*vv;
      mixion mi;
      vector<string> level2;
      h=0;
      for (p=0; p<m.length(); p++){
	if (m[p]==','){
	  level2.push_back(m.substr(h,p-h));
	  h=p+1;
	}
      }
      level2.push_back(m.substr(h,1000));
      mi.type=level2[0];
      mi.T=atof(level2[1].c_str());
      mi.ratio=atof(level2[2].c_str());
      ionbag.push_back(mi);
    }
    float total=0;
    vector<mixion>::iterator vmi;
    for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++)
      total+=vmi->ratio;
    float sum=0;
    for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++){
      vmi->ratio/=total;
      vmi->ratio+=sum;
      sum=vmi->ratio;
      vmi->type=RectifySpecies(vmi->type);
    }
    config cfg(scfg);
    list<int>::iterator seed=seeds.begin();
    int kept=0;
    //pick the first one. 
    srand(TimeSeed()-1); //don't want to reseed with the same number.
    float fff=rand01();
    for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++){
      if (fff < vmi->ratio) break;
    }
    for (int run=begin; run<=begin+runs-1; run++){
      cerr<<"------------"<<DateTime()<<"--------run "<<run
	  <<"-----------------------"<<endl;
      int No=cfg.N;
      VAI I;
      if (vmi->type=="Ar")
	I=cfg.AddIon(50, "Ar", vmi->T, 0, 0, 0, 0, *seed, 0);
      else
	I=cfg.AddIon(e_ion, vmi->type, vmi->T, choose_inc, inc, choose_azi, 
		       azi, *seed, 0);
      for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;
      int Hix=I->ix;
      double nbegin=I->NC;
      bool quit=0;
      if (vmi->type=="Ar"){
	cfg.Run(0.5, 0, dt, printe, 0, -1, 0, tso,0);
	cfg.DelIon(0);
	cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
	if (dump_only && cfg.N!=No) cfg.Dump();
      }
//       else if (vmi->type=="H"){
// 	while (!quit){
// 	  //step 10 fs at a time.
// 	  cfg.Run(0.01, 0, dt, printe, 0, -1, 0, tso,0);
// 	  if (cfg.atomix(Hix)==cfg.end){
// 	    quit=1;
// 	    cfg.Nmax--;
// 	    if (cfg.N==No) cerr<<"* Species was scattered.\n";
// 	    else{ 
// 	      cerr<<"* Species abstracted something.\n";
// 	    }
// 	  }
// 	  else if (I->NC > 0.5){
// 	    cerr<<"* Species appears to be bound. "<<I->R<<endl;
// 	    cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
// 	    //make sure it's not going to etch.
// 	    cfg.Run(0.5, 0, dt, printe, 0, -1, 0, tso,0);
// 	    if (cfg.N==No) kept++;
// 	    quit=1;
// 	  }
// 	}
//       }
      else{
	while (!quit){
	  //step 10 fs at a time.
	  cfg.Run(0.01, 0, dt, printe, 0, -1, 0, tso,0);
	  I=cfg.atomix(Hix);
	  if (I==cfg.end){
	    if (vmi->type=="H") cfg.Nmax--;
	    if (vmi->type=="CH") cfg.Nmax-=2;
	    if (vmi->type=="CH2") cfg.Nmax-=3;
	    if (vmi->type=="CH3") cfg.Nmax-=4;
	    quit=1;
	    if (cfg.N==No) cerr<<"* Species was scattered.\n";
	    else{ 
	      cerr<<"* Species abstracted something.\n";
	      if (dump_only) cfg.Dump();
	    }
	    if (fabs(1-T/cfg.T) > 0.05)
	      cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0);
	  }
	  else if (I->NC > 0.3){
	    cerr<<"* Species appears to be bound. "<<I->R<<endl;
	    cerr<<"  Begin: "<<nbegin<<" End: "<<I->NC<<endl;
	    kept++;
	    quit=1;
	    cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
	    if (dump_only) cfg.Dump();
	  }
	}
      }
      cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
      if (dump_each) cfg.Dump();
      seed++;
      //pick the next one. this avoids reseeding the RNG with the same number
      fff=rand01();
      for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++){
	if (fff < vmi->ratio) break;
      }
    }
    if (dump_end || dump_only) cfg.Dump();
    cerr<<"---------------------------------------------------------------------------"<<endl;
    cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
  }
  else cerr<<"Unknown species ("<<species<<").\n";
  exit(0);
}

