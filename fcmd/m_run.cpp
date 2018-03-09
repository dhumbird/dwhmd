#include "config.h"
#include <string>
#include "errors.h"

void MainRun(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  double runtime=0; double To=0; double dt=0; double bhb_dt=0;
  int cfgevery=0; int track=-1; bool tso=0;
  bool printe=0; bool dump_end=0; bool quiet=0; short sput=1;
  string arg; bool quench=0; 
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-T") To=atof(argv[++i]);
    else if (arg=="-tso") tso=1;
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-track") track=atoi(argv[++i]);
    else if (arg=="-dump_end") dump_end=1;
    else if (arg=="-quiet") quiet=1;
    else if (arg=="-nosput") sput=0;
    else if (arg=="-sweep") sput=2;
    else if (arg=="-quench") {quench=1; To=atof(argv[++i]);}
    else CmdError(arg.c_str());
  }
  config cfg(scfg);
  if (quench)
    cfg.RunQuenched(runtime, cfgevery, dt, printe, dump_end, track, quiet, 
		    tso, sput, To);
  else{
    if (To==0)
      cfg.Run(runtime, cfgevery, dt, printe, dump_end, track, quiet, tso, 
              sput);
    else
      cfg.Thermo(To, bhb_dt, cfgevery, dt, printe, dump_end, tso,runtime);
  }
  exit(0);
}


void MainBondDist(int argc, char * argv[]){
  string scfg="temp_0000-000.cfg";
  double runtime=1; double To=0; double dt=0; double bhb_dt=0;
  int cfgevery=0; int track=-1; bool tso=0;
  bool printe=0; bool dump_end=0; bool quiet=0; short mode=0;
  string arg; bool quench=0; 
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-T") To=atof(argv[++i]);
    else if (arg=="-tso") tso=1;
    else if (arg=="-cfgevery") cfgevery=atoi(argv[++i]);
    else if (arg=="-dt") dt=atof(argv[++i]);
    else if (arg=="-bhb_dt") bhb_dt=atof(argv[++i]);
    else if (arg=="-printe") printe=1;
    else if (arg=="-track") track=atoi(argv[++i]);
    else if (arg=="-dump_end") dump_end=1;
    else if (arg=="-quiet") quiet=1;
    else if (arg=="-a") mode=1;
    else if (arg=="-quench") {quench=1; To=atof(argv[++i]);}
    else CmdError(arg.c_str());
  }
  config cfg(scfg);
  cfg.bonddist(runtime, dt, tso, mode);
  exit(0);
}
