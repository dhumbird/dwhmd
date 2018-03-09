//md_main.cpp
//An object-oriented molecular dynamics code in C++.
//By Dave Humbird

#include "config.h"
#include "errors.h"
#include "edep.h"
#include "bombard.h"
//#include "edist.h"

int main(int argc, char * argv[]){
  int i;

  if (argc==1){
    cerr<<"supply some more arguments.\n";
    exit(1);
  }
  
  string arg;

  for (i=0; i<argc; i++){
    arg=argv[i];
    cerr<<arg<<" ";
  }
  cerr<<endl;
  
  arg=argv[1];

  if (arg=="-init"){
    short Ix=4; short Iy=4; short Iz=4; short id=14;
    float T=300; float rho=0.049852723;
    int seed=0;
    bool open=0; bool isolated=0;
    string name="temp";

    for (i=2; i<argc; i++){
      arg=argv[i];
      if(arg=="-T") T=atof(argv[++i]);
      else if(arg=="-I") { Ix=atoi(argv[++i]); Iy=Ix; Iz=Ix;}
      else if(arg=="-Ix") Ix=atoi(argv[++i]);
      else if(arg=="-Iy") Iy=atoi(argv[++i]);
      else if(arg=="-Iz") Iz=atoi(argv[++i]);
      else if(arg=="-rho") rho=atof(argv[++i]);
      else if(arg=="-of" || arg=="-name") name=argv[++i];
      else if(arg=="-open") open=1;
      else if(arg=="-isolated") isolated=1;
      else if(arg=="-mat") id=atoi(argv[++i]);
      else if(arg=="-seed") seed=atoi(argv[++i]);
      else CmdError(arg.c_str());
    }
    config atom;
    atom.MakeCryst(T, Ix, Iy, Iz, rho, open, isolated, id, seed, name);
    atom.Dump(0);
    exit(0);
  }
  
  else if (arg=="-run"){
    string scfg="temp_000000.cfg";
    double runtime=0; double To=0; double dt=0; double bhb_dt=0;
    int cfgevery=0;
    bool printe=0; bool tso=0; bool dump_end=0;
    
    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
      else if (arg=="-time") runtime=atof(argv[++i]);
      else if (arg=="-T") To=atof(argv[++i]);
      else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
      else if (arg=="-dt") dt=atof(argv[++i]);
      else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
      else if (arg=="-printe") printe=1;
      else if (arg=="-tso") tso=1;
      else if (arg=="-dump_end") dump_end=1;
      else CmdError(arg.c_str());
    }
    config cfg(scfg);
    if (To==0.0) cfg.Run(runtime, cfgevery, dt, printe, tso, dump_end);
    else cfg.Thermo(To, bhb_dt, cfgevery, dt, printe, tso, dump_end);
    exit(0);
  }

  else if (arg=="-addion"){
    string scfg="temp_000000.cfg";
    float e=100; float T=300; float inc=0; float azi=0; bool confine=0;
    bool round=0; float conf=0;
    short id=18;
    int seed=0;
    bool choose_inc=0; bool choose_azi=0; bool no_noise=0;
    bool choose_pos=1;

    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
      else if (arg=="-e") e=atof(argv[++i]);
      else if (arg=="-id") id=atoi(argv[++i]);
      else if (arg=="-T") T=atof(argv[++i]);
      else if (arg=="-seed") seed=atoi(argv[++i]);
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
      else if (arg=="-norand") choose_pos=0;
      else if (arg=="-nonoise") no_noise=1;
      else if (arg=="-confine"){
	confine=1;
	arg=argv[++i];
	if (arg=="round") round=1;
	else if (arg=="square") round=0;
	else CmdError(arg.c_str());
	conf=atof(argv[++i]);
      }
      else CmdError(arg.c_str());
    }
    config atom(scfg);
    atom.AddIon(e, id, T, choose_inc, inc, choose_azi, azi, seed, 0, no_noise,
		choose_pos, confine, round, conf);
    atom.Dump(0);
    exit(0);
  }
  
  else if (arg=="-delion" || arg=="-rmion"){
    string scfg="temp_000000.cfg";

    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
    }
    config cfg(scfg);
    cfg.DelIon(0);
    cfg.Dump(0);
    exit(0);
  }
  
  else if (arg=="-reset"){
    string scfg="temp_000000.cfg";
    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
    }
    config cfg(scfg);
    cfg.t=0;
    cfg.Dump(0);
  }
  
  else if (arg=="-edep") EnergyDep(argc, argv);
  else if (arg=="-bomb") IonBombard(argc, argv);
  //  else if (arg=="-edist") EvolvEnergyDist(argc, argv);

  else if (arg=="-sortz"){
    string scfg="temp_000000.cfg";

    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
    }
    config cfg(scfg);
    cfg.Rz_sort_asc();
    cfg.Dump(0);
    exit(0);
  }

  else if (arg=="-textdump"){
    string scfg="temp_000000.cfg";
    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
    }
    config cfg(scfg);
    cfg.Dump(1);
    exit(0);
  }

  else if (arg=="-tellme"){
    string scfg="temp_000000.cfg";
    int ix=0;
    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
      else ix=atoi(arg.c_str());
    }
    config cfg(scfg);
    cout<<*cfg.atomix(ix)<<endl;
  }

  else if (arg=="-delete"){
    string scfg="temp_000000.cfg";
    int ix=0;
    for (i=2; i<argc; i++){
      arg=argv[i];
      if (sfind(arg, ".cfg")) scfg=arg;
      else ix=atoi(arg.c_str());
    }
    config cfg(scfg);
    if (cfg.erase(cfg.atomix(ix))) cfg.Dump(0);
    else cerr<<"Atom "<<ix<<" not found!\n";
  }

  else CmdError(arg.c_str());

};




