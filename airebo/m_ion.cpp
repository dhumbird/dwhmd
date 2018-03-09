#include "config.h"
#include "errors.h"
#include <string>
#include "strfun.h"

void MainAddion(int argc, char * argv[]){
  string scfg="temp_000000.cfg";
  float e=100; float T=300; float inc=0; float azi=0;
  string species="unknown";
  int seed=0;
  bool choose_inc=0; bool choose_azi=0;
  string arg; bool nosave=0;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-e") e=atof(argv[++i]);
    else if (arg=="-type") species=argv[++i];
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
    else if (arg=="-nosave") nosave=1;
    else CmdError(arg.c_str());
  }
  config atom(scfg);
  species=RectifySpecies(species);
  atom.AddIon(e, species, T, choose_inc, inc, choose_azi, azi, seed, 0);
  if (!nosave) atom.Dump();
  exit(0);
}

// void MainDelion(int argc, char * argv[]){
//   string scfg="temp_000000.cfg";
//   string arg;
//   for (int i=2; i<argc; i++){
//     arg=argv[i];
//     if (sfind(arg, ".cfg")) scfg=arg;
//   }
//   config cfg(scfg);
//   cfg.DelIon(0);
//   cfg.Dump();
//   exit(0);
// }
