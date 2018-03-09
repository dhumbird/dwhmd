#include "config.h"
#include <string>

void MainSortek(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
  }
  config cfg(scfg);
  cfg.Ek_sort_asc();
  cfg.Dump();
  exit(0);
}

void MainSortix(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
  }
  config cfg(scfg);
  cfg.ix_sort_asc();
  cfg.Dump();
  exit(0);
}

void MainSortz(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  string arg;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
  }
  config cfg(scfg);
  cfg.Rz_sort_asc();
  cfg.Dump();
  exit(0);
}
