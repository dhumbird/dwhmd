#include "config.h"
#include "errors.h"

void MainSat(int argc, char * argv[]){
  string arg;
  string scfg="temp_00000000.cfg";
  short id=9;
  int rounds=10;
  int begin=1;
  float T=300;
  float runtime=0.5;
  for (int i=2; i<argc; i++){
    arg=argv[i];
    if (sfind(arg, ".cfg")) scfg=arg;
    else if (arg=="-id") id=atoi(argv[++i]);
    else if (arg=="-rounds") rounds=atoi(argv[++i]);
    else if (arg=="-time") runtime=atof(argv[++i]);
    else if (arg=="-begin") begin=atoi(argv[++i]);
    else if (arg=="-T") T=atof(argv[++i]);
  }
  config cfg(scfg);
  VPI i, target;
  for (int run=begin; run<=begin+rounds-1; run++){
    cerr<<"------------"<<DateTime()<<"--------round "<<run
	<<"-----------------------"<<endl;
    svector P;
    double u=cfg.FindWell(id,&P);
    VPI test=cfg.append(id);
    test->R=P;
    cerr<<"* Adding atom "<<test->ix<<" at "<<test->R<<endl;
    cerr<<"* Well depth was "<<u<<" eV."<<endl;
    cfg.RunQuenched(runtime,0,0,0,0,-1,0,0,0,T);
    cfg.Dump();
  }
  cerr<<"---------------------------------------------------------------------------"<<endl;
  cerr<<DateTime()<<" Requested rounds ("<<rounds<<") completed."<<endl;
}
  
