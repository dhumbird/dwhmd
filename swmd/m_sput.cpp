#include "config.h"

void MainSweep(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  bool dump=1;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-nodump") dump=0;
  }
  config cfg(scfg);
  cfg.Sweep();
  if (dump)  cfg.Dump();
  exit(0);
}

void MainChksput(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  bool dump=1;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-nodump") dump=0;
  }
  config cfg(scfg);
  cfg.CheckSput();
  if (dump) cfg.Dump();
  exit(0);
}
