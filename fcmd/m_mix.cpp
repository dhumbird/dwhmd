#include "config.h"
#include "errors.h"
#include "consts.h"
#include "sysfun.h"
#include "strfun.h"
#include <list>

struct mixion{
  float ratio;
  float T;
  float e;
  string type;
};

void MainMixIon(int argc, char * argv[]){
  
  int i;
  string arg;
  
  //*******************parse & setup input*********************************
  //input identifiers
  string scfg="temp_0000-000.cfg";
  float runtime=0.5;
  float inc=0; float azi=0; float T=300; 
  double dt=0; double bhb_dt=0;
  int cfgevery=0; int runs=1; int ion_id=18; int begin=1;
  bool printe=0; bool dump_end=0; bool tso=0;
  bool choose_inc=0; bool choose_azi=0; bool dump_only=0;
  short cfgthermo=1; bool pickup=0; bool dump_before=0;
  bool dump_each=0;
  list<int> seeds;
  string mix="";
  bool b_test=1;
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
    else if (arg=="-T") T=atof(argv[++i]);
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-pickup") pickup=1;
    else if (arg=="-tso") tso=1;
    else if (arg=="-notest") b_test=0;
    else if (arg=="-nocfgthermo") cfgthermo=0;
    else if (arg=="-confine") confine=atof(argv[++i]);
    else if (arg=="-dump"){
      arg=argv[++i];
      if (arg=="end") dump_end=1;
      else if (arg=="each") dump_each=1;
      else if (arg=="only") dump_only=1;
      else if (arg=="before") dump_before=1;
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
    else if (arg=="-ions") {
      mix=argv[++i];
    }
    else CmdError(arg.c_str());
  }

  if (seeds.empty()) 
    for (int j=0; j<runs; j++)
      seeds.push_back(0);
  
  //***********************MIXED IONS***********************************
  //command line looks like -mix species1,T1,e1,ratio1:species2,T2,e2,ratio2
  
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
  level1.push_back(mix.substr(h,string::npos));
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
    level2.push_back(m.substr(h,string::npos));
    mi.type=RectifySpecies(level2[0]);
    mi.T=atof(level2[1].c_str());
    mi.e=atof(level2[2].c_str());
    mi.ratio=atof(level2[3].c_str());
    ionbag.push_back(mi);
  }
  float total=0;
  vector<mixion>::iterator vmi;
  for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++){
    total+=vmi->ratio;
  }
  float sum=0;
  for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++){
    vmi->ratio/=total;
    vmi->ratio+=sum;
    sum=vmi->ratio;
    vmi->type=RectifySpecies(vmi->type);
  }
  (ionbag.back()).ratio=1;
  config cfg(scfg);
  list<int>::iterator seed=seeds.begin();
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
    I=cfg.AddIon(vmi->e, vmi->type, vmi->T, choose_inc, inc, choose_azi, 
	       azi, *seed, 0, confine);
    for (VAI s=cfg.begin; s<cfg.end; s++) s->P=s->R;

    if (vmi->e != 0){
      cfg.Run(runtime, 0, dt, printe, 0, -1, 0, tso,1);
      cfg.DelIon(0);
      cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, tso, 0.5);
      //sweep
      cfg.Run(0.05, 0, dt, printe, 0, -1, 0, tso,2);
      if (dump_only) cfg.Dump();
    }
    else{
      cfg.Setdt(0.001);	
      int Hix=I->ix;
      double nbegin=I->Nmap[6]+I->Nmap[14];  
      bool quit=0;
      while (!quit){
	//step 10 fs at a time.
	cfg.Run(0.01, 0, dt, printe, 0, -1, 0, 0,1);
	I=cfg.atomix(Hix);
	if (I==cfg.end){
	  quit=1;
	  cerr<<"* Species left.\n";
	  if (cfg.N != No && dump_only) cfg.Dump();
	  if (fabs(1-T/cfg.T) > 0.05)
	    cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, 0, 0);
	}
	else if (I->Nmap[6]+I->Nmap[14] > 0.3){
	  cerr<<"* Species appears to be bound. "<<I->R<<endl;
	  cerr<<"  Begin: "<<nbegin<<" End: "<<I->Nmap[6]+I->Nmap[14]<<endl;
	  quit=1;
	  cfg.Thermo(T, bhb_dt, 0, dt, printe, 0, 0, 0.5);
	  //sweep
	  cfg.Run(0.05, 0, dt, printe, 0, -1, 0, 0 ,2);
	  if (dump_only) cfg.Dump();
	}
      }
    }
    cerr<<"* Impact complete. Time: "<<cfg.t<<endl;
    cfg.CheckDepth();
    seed++;
    //pick the next one. this avoids reseeding the RNG with the same number
    fff=rand01();
    for (vmi=ionbag.begin(); vmi!=ionbag.end(); vmi++){
      if (fff < vmi->ratio) break;
    }
    if ((vmi->type=="Ar" && dump_before) || dump_each) cfg.Dump();
  }
  if (dump_end || dump_only) cfg.Dump();
  cerr<<"---------------------------------------------------------------------------"<<endl;
  cerr<<DateTime()<<" Requested runs ("<<runs<<") completed."<<endl;
}
