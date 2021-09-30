//md_main.cpp
//An object-oriented molecular dynamics code in C++.
//By Dave Humbird

#include "md_main.h"

int main(int argc, char * argv[]){
  if (argc==1){
    cerr<<"supply some more arguments.\n";
    exit(1);
  }
  //echo the command line
  for (int i=0; i<argc; i++) cerr<<argv[i]<<" ";
  cerr<<endl;
  
  
  Pcc_genCoef(); Pcf_genCoef(); Psif_genCoef();
  Fcc_genCoef(); Psicl_genCoef();
  
  string arg=argv[1];

  if (arg=="-run") MainRun(argc, argv);
  else if (arg=="-addfix") MainAddFixed(argc,argv);
  else if (arg=="-reset") MainReset(argc, argv);
  else if (arg=="-delete") MainDelete(argc, argv);
  else if (arg=="-seto") MainShiftOrigin(argc,argv);
  else if (arg=="-rand") MainRandomize(argc,argv);
  else if (arg=="-init") MainInit(argc,argv);
  else if (arg=="-edep") MainEnergyDep(argc, argv);
  //else if (arg=="-sat") MainSat(argc,argv);
  //else if (arg=="-pass") MainPass(argc,argv);
  else if (arg=="-addion") MainAddion(argc, argv);
  //else if (arg=="-path") MainPath(argc,argv);
  else if (arg=="-bomb") MainIonBombard(argc,argv);
  else if (arg=="-mixion") MainMixIon(argc,argv);
  else if (arg=="-cgmin") MainCGmin(argc,argv);
  else if (arg=="-sweep") MainMisc(argc,argv);
  else if (arg=="-flist") MainForceList(argc,argv);

  else CmdError(arg.c_str());
};
