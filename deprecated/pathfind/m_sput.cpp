#include "config.h"

void MainSweep(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
  }
  config cfg(scfg);
  while(!cfg.Sweep()){}
  cfg.Dump();
  exit(0);
}

void MainChksput(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
  }
  config cfg(scfg);
  cfg.CheckSput();
  cfg.Dump();
  exit(0);
}
