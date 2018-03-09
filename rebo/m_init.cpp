
#include "config.h"
#include "errors.h"
#include <string>

void MainInit(int argc, char * argv[]){
  short Ix=4; short Iy=4; short Iz=4; short id=14;
  float T=300; float rho=0;
  int seed=0;
  string name="temp";
  string arg;

  for (int i=2; i<argc; i++){
    arg=argv[i];
    if(arg=="-T") T=atof(argv[++i]);
    else if(arg=="-I") { Ix=atoi(argv[++i]); Iy=Ix; Iz=Ix;}
    else if(arg=="-Ix") Ix=atoi(argv[++i]);
    else if(arg=="-Iy") Iy=atoi(argv[++i]);
    else if(arg=="-Iz") Iz=atoi(argv[++i]);
    else if(arg=="-rho") rho=atof(argv[++i]);
    else if(arg=="-name") name=argv[++i];
    else if(arg=="-mat") id=atoi(argv[++i]);
    else if(arg=="-seed") seed=atoi(argv[++i]);
    else CmdError(arg.c_str());
  }
  config atom;
  atom.MakeCryst(T, Ix, Iy, Iz, rho, id, seed, name);
  atom.Dump();
  exit(0);
}

//  void MainAddFixed(int argc, char * argv[]){
//    string arg; string scfg="temp_00000000.cfg";
//    float rho=0;
//    for (int i=2; i<argc; i++){
//      arg=argv[i];
//      if (sfind(arg,".cfg")) scfg=arg;
//      else if (arg=="-rho") rho=atof(argv[++i]);
//      else CmdError(arg.c_str());
//    }
//    config cfg(scfg);
//    cfg.AddFixed(rho);
//    cfg.Dump();
//  }
